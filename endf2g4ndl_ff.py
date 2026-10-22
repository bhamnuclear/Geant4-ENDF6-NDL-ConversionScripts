#!/usr/bin/env python3
"""Convert ENDF-6 MF=8 fission product yields into a G4NDL Fission/FF/ file.

The G4NDL "FF" files are read by G4ENDFTapeRead::ReadInData() with plain
stream extraction, so the record layout is:

    [G4NDL]                       <- optional, consumed by GetDataStream()
    [<library name>]              <- optional, only read if the first token is G4NDL
    MT  MF  dummy  nEnergyBlocks
      E_inc(eV)  nProducts  interpolation
        ZZZAAA  metastate  yield      x nProducts
      ... x nEnergyBlocks

MT must be 454 (independent yields) for G4WendtFissionFragmentGenerator, which
is the class that consumes these files.

Usage:
    endf2g4ndl_ff.py <endf-tape> <output-file> [--mt 454] [--library NAME]
                     [--energies E1,E2,...]
"""

import argparse
import sys

FIELD_WIDTH = 11
FIELDS_PER_LINE = 6


def endf_float(s):
    """Parse an ENDF-6 real number.

    Handles both the explicit form (1.2345E+04) and the compact form in which
    the exponent marker is dropped (1.2345+04).
    """
    s = s.strip()
    if not s:
        return 0.0
    if "e" not in s and "E" not in s:
        # Skip a leading sign, then look for the exponent sign.
        for i in range(1, len(s)):
            if s[i] in "+-" and s[i - 1] not in "eE":
                s = s[:i] + "E" + s[i:]
                break
    return float(s)


def endf_int(s):
    s = s.strip()
    return int(s) if s else 0


class Tape:
    """Line-oriented reader over one MF/MT section of an ENDF-6 tape."""

    def __init__(self, path, mf, mt):
        self.lines = []
        with open(path, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if len(line) < 75:
                    continue
                if endf_int(line[70:72]) == mf and endf_int(line[72:75]) == mt:
                    self.lines.append(line[:66])
        self.pos = 0

    def record(self):
        """Return the six 11-character fields of the next line, as raw strings."""
        line = self.lines[self.pos].ljust(66)
        self.pos += 1
        return [line[i * FIELD_WIDTH:(i + 1) * FIELD_WIDTH] for i in range(FIELDS_PER_LINE)]

    def control(self):
        """[C1, C2, L1, L2, N1, N2] of the next record."""
        f = self.record()
        return (endf_float(f[0]), endf_float(f[1]),
                endf_int(f[2]), endf_int(f[3]), endf_int(f[4]), endf_int(f[5]))

    def values(self, count):
        """Read `count` floats from the data lines that follow a LIST head."""
        out = []
        while len(out) < count:
            for field in self.record():
                out.append(endf_float(field))
                if len(out) == count:
                    break
        return out


def read_yields(path, mf=8, mt=454):
    """Return (za, awr, [(energy_eV, [(zafp, state, yield), ...]), ...])."""
    tape = Tape(path, mf, mt)
    if not tape.lines:
        raise SystemExit(f"{path}: no MF={mf} MT={mt} section found")

    za, awr, nenergy, _, _, _ = tape.control()

    blocks = []
    for _ in range(nenergy):
        energy, _, _, _, nn, nfp = tape.control()
        data = tape.values(nn)
        products = []
        for j in range(nfp):
            zafp, fps, y, _dy = data[4 * j:4 * j + 4]
            products.append((int(round(zafp)), int(round(fps)), y))
        blocks.append((energy, products))

    return int(round(za)), awr, blocks


def write_g4ndl(path, mt, mf, awr, blocks, library, interpolation=2):
    with open(path, "w") as out:
        out.write("G4NDL\n")
        out.write(f"{library}\n")
        out.write(f"{mt:14d}\n")
        out.write(f"{mf:14d}\n")
        out.write(f"{awr:14.6e}\n")
        out.write(f"{len(blocks):14d}\n")
        for energy, products in blocks:
            out.write(f"{energy:14.6e}\n")
            out.write(f"{len(products):14d}\n")
            out.write(f"{interpolation:14d}\n")
            for i, (zafp, state, yld) in enumerate(products):
                out.write(f"{zafp:14d}{state:14d}{yld:14.6e}")
                out.write("\n" if i % 2 else "")
            if len(products) % 2:
                out.write("\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("endf", help="ENDF-6 tape containing MF=8 fission yields")
    p.add_argument("output", help="G4NDL Fission/FF file to write")
    p.add_argument("--mt", type=int, default=454,
                   help="454 = independent yields (default), 459 = cumulative")
    p.add_argument("--library", default="UKFY4.1",
                   help="library label written on line 2")
    p.add_argument("--energies",
                   help="comma-separated list of incident energies (eV) to keep; "
                        "the nearest available point is used for each")
    p.add_argument("--compress", action="store_true",
                   help="also write <output>.z, the zlib stream G4NDL ships and "
                        "G4ParticleHPManager::GetDataStream() prefers")
    args = p.parse_args()

    za, awr, blocks = read_yields(args.endf, mt=args.mt)

    if args.energies:
        wanted = [float(e) for e in args.energies.split(",")]
        keep = sorted({min(range(len(blocks)),
                           key=lambda i: abs(blocks[i][0] - e)) for e in wanted})
        blocks = [blocks[i] for i in keep]

    write_g4ndl(args.output, args.mt, 8, awr, blocks, args.library)

    if args.compress:
        import zlib
        with open(args.output, "rb") as fh:
            raw = fh.read()
        with open(args.output + ".z", "wb") as fh:
            fh.write(zlib.compress(raw, 9))

    z, a = divmod(za, 1000)
    print(f"ZA={za} (Z={z}, A={a})  AWR={awr:.5f}", file=sys.stderr)
    print(f"{len(blocks)} incident-energy points -> {args.output}", file=sys.stderr)
    for energy, products in blocks:
        total = sum(y for _, _, y in products)
        print(f"  E = {energy:12.5e} eV   {len(products):5d} products   "
              f"sum(Y) = {total:.6f}", file=sys.stderr)


if __name__ == "__main__":
    main()
