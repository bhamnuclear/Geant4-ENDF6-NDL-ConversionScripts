#!/usr/bin/env python3
"""Convert an ENDF-6 neutron evaluation into the G4NDL Fission/ final-state files.

Companion to endf2g4ndl_ff.py, which handles Fission/FF/. This script writes
every other directory under Fission/:

    CrossSection/  MF=3  MT=18            fission cross section
    FS/            MF=1  MT=452,455,456,458  nu-bar, delayed-neutron data, energy release
                   MF=4  MT=18            prompt-neutron angular distribution
                   MF=5  MT=18,455        prompt / delayed neutron spectra
                   MF=12,14,15 MT=18      fission photons
    FC/ SC/ TC/ LC/  MF=3,4,5  MT=19,20,21,38   first/second/third/last chance

Every G4NDL file is read with plain stream extraction, so the layout is just
the ENDF-6 records re-emitted in the token order the Geant4 readers expect.
The readers this mirrors (geant4-v11.4.2, particle_hp/):

    G4ParticleHPIsoData::Init            CrossSection/
    G4ParticleHPFSFissionFS::Init        FS/
    G4ParticleHPFissionBaseFS::Init      FC/ SC/ TC/ LC/
    G4ParticleHPVector::Init             every TAB1: NP, NR, (NBT INT)*NR, (x y)*NP

Cross sections must already be pointwise and lin-lin, i.e. the tape must have
been through PREPRO (LINEAR + RECENT), as described in INDC(NDS)-0758 section 2.
Everything else is copied from the tape as-is: Geant4 honours the ENDF
interpolation laws for distributions and multiplicities.

Usage:
    endf2g4ndl_fission.py <endf-tape>... <output-dir> [--only FS,FC,...] [--compress]
                          [--digits 7]
"""

import argparse
import os
import sys
import zlib

FIELD_WIDTH = 11

# G4ParticleHPNames::theString, spelling included ("Platinium", "Berylium").
ELEMENTS = [
    "Hydrogen", "Helium", "Lithium", "Berylium", "Boron", "Carbon", "Nitrogen",
    "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium", "Aluminum", "Silicon",
    "Phosphorous", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium",
    "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt",
    "Nickel", "Copper", "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium",
    "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium",
    "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium",
    "Silver", "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine",
    "Xenon", "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium",
    "Neodymium", "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium",
    "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium", "Lutetium",
    "Hafnium", "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinium",
    "Gold", "Mercury", "Thallium", "Lead", "Bismuth", "Polonium", "Astatine",
    "Radon", "Francium", "Radium", "Actinium", "Thorium", "Protactinium",
    "Uranium", "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium",
    "Californium", "Einsteinium", "Fermium",
]

CHANCE_DIRS = {"FC": 19, "SC": 20, "TC": 21, "LC": 38}
ALL_DIRS = ["CrossSection", "FS", "FC", "SC", "TC", "LC"]


class ConversionError(Exception):
    pass


# --------------------------------------------------------------------------
# ENDF-6 reading
# --------------------------------------------------------------------------

def endf_float(s):
    """Parse an ENDF-6 real, with or without the exponent marker (1.2345+04)."""
    s = s.strip()
    if not s:
        return 0.0
    if "e" not in s and "E" not in s:
        for i in range(1, len(s)):
            if s[i] in "+-" and s[i - 1] not in "eE":
                s = s[:i] + "E" + s[i:]
                break
    return float(s)


def endf_int(s):
    s = s.strip()
    return int(s) if s else 0


class Tab1:
    def __init__(self, c1, c2, l1, l2, nbt, interp, x, y):
        self.c1, self.c2, self.l1, self.l2 = c1, c2, l1, l2
        self.nbt, self.interp, self.x, self.y = nbt, interp, x, y


class Section:
    """Sequential record reader over one MF/MT section (lines without the SEND)."""

    def __init__(self, lines):
        self.lines = lines
        self.pos = 0

    def _fields(self):
        line = self.lines[self.pos].ljust(66)
        self.pos += 1
        return [line[i * FIELD_WIDTH:(i + 1) * FIELD_WIDTH] for i in range(6)]

    def cont(self):
        f = self._fields()
        return (endf_float(f[0]), endf_float(f[1]),
                endf_int(f[2]), endf_int(f[3]), endf_int(f[4]), endf_int(f[5]))

    def _floats(self, n):
        out = []
        while len(out) < n:
            for field in self._fields():
                out.append(endf_float(field))
                if len(out) == n:
                    break
        return out

    def _ints(self, n):
        out = []
        while len(out) < n:
            for field in self._fields():
                out.append(endf_int(field))
                if len(out) == n:
                    break
        return out

    def list(self):
        c1, c2, l1, l2, npl, n2 = self.cont()
        return (c1, c2, l1, l2, npl, n2), self._floats(npl)

    def _interp(self, nr):
        v = self._ints(2 * nr)
        return v[0::2], v[1::2]

    def tab1(self):
        c1, c2, l1, l2, nr, np_ = self.cont()
        nbt, interp = self._interp(nr)
        v = self._floats(2 * np_)
        return Tab1(c1, c2, l1, l2, nbt, interp, v[0::2], v[1::2])

    def tab2(self):
        c1, c2, l1, l2, nr, nz = self.cont()
        nbt, interp = self._interp(nr)
        return (c1, c2, l1, l2, nz), nbt, interp


def read_tape(path):
    """Return {(mf, mt): [66-column lines]} for the first material on the tape.

    Several tapes of the same material may be concatenated (e.g. a PREPRO
    output followed by the ENDFB.MF4/MF5/M12/M14/M15 files SIXPAK writes);
    a section appearing again replaces the earlier copy.
    """
    sections = {}
    mat = None
    current = None
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if len(line) < 75:
                continue
            m, mf, mt = endf_int(line[66:70]), endf_int(line[70:72]), endf_int(line[72:75])
            if m <= 0 or mf == 0 or mt == 0:
                current = None    # TPID, SEND, FEND, MEND, TEND
                continue
            if mat is None:
                mat = m
            elif m != mat:
                continue          # only the first material
            if (mf, mt) != current:
                current = (mf, mt)
                sections[current] = []
            sections[current].append(line[:66])
    if not sections:
        raise ConversionError(f"{path}: no ENDF-6 material found")
    return mat, sections


# --------------------------------------------------------------------------
# G4NDL writing
# --------------------------------------------------------------------------

# Significant digits for reals. 7 reproduces the IAEA files; above ~1 keV that
# cannot separate points inside narrow resolved resonances (--digits).
DIGITS = 7


def fe(v):
    return f"{v:.{DIGITS - 1}e}"


def fg(v):
    return f"{v:g}"


class G4Writer:
    """Accumulates G4NDL text in the layout of the IAEA-converted files:
    every token preceded by ' \\t', data pairs three to a line."""

    def __init__(self):
        self.lines = []

    def row(self, *tokens):
        self.lines.append("".join(f" \t{t}" for t in tokens))

    def floats(self, values, per_line=6):
        for i in range(0, len(values), per_line):
            self.row(*(fe(v) for v in values[i:i + per_line]))

    def pairs(self, x, y):
        flat = []
        for a, b in zip(x, y):
            flat += [a, b]
        self.floats(flat)

    def interp(self, nbt, interp):
        """G4InterpolationManager::Init: NR, then NBT INT per range."""
        self.row(len(nbt))
        for n, i in zip(nbt, interp):
            self.row(n, i)

    def vector(self, t):
        """G4ParticleHPVector::Init(istream): NP, interpolation, pairs."""
        self.row(len(t.x))
        self.interp(t.nbt, t.interp)
        self.pairs(t.x, t.y)

    def text(self):
        return "\n".join(self.lines) + "\n"


def is_linear(t):
    return all(i == 2 for i in t.interp) or len(t.x) < 2


def write_xs_pairs(w, t):
    """Cross sections are read without an interpolation table: lin-lin only."""
    w.pairs(t.x, t.y)


def write_angular(w, sec, awr_default):
    """MF=4 -> G4ParticleHPAngular::Init."""
    s = Section(sec)
    _za, awr_head, _l1, ltt, _n1, _n2 = s.cont()
    _c1, awr, li, lct, _nk, _nm = s.cont()
    # The AWR from the second record is passed through, as in the IAEA files, even
    # where an evaluation leaves it 0 (JEFF-3.3 236-238U, 240Pu). Geant4 only uses
    # it to leave the centre-of-mass frame, so substitute a real mass for LCT=2.
    if lct == 2 and awr == 0:
        awr = awr_head or awr_default
    if li == 1 or ltt == 0:
        w.row(0, fg(awr), lct)
        return
    if ltt == 1:
        (_, _, _, _, ne), nbt, interp = s.tab2()
        w.row(1, fg(awr), lct)
        w.row(ne)
        w.interp(nbt, interp)
        for _ in range(ne):
            (t, e, lt, _l2, _nl, _n2), coeff = s.list()
            w.row(fe(t), fe(e), lt, len(coeff))
            w.floats(coeff)
        return
    if ltt == 2:
        (_, _, _, _, ne), nbt, interp = s.tab2()
        w.row(2, fg(awr), lct)
        w.row(ne)
        w.interp(nbt, interp)
        for _ in range(ne):
            f = s.tab1()
            w.row(fe(f.c1), fe(f.c2), f.l1)
            w.vector(f)
        return
    raise ConversionError(f"MF=4 LTT={ltt} has no Geant4 representation "
                          "(G4ParticleHPAngular knows isotropic, Legendre, tabulated)")


def write_energy(w, sec):
    """MF=5 -> G4ParticleHPEnergyDistribution::Init."""
    s = Section(sec)
    _za, awr, _l1, _l2, nk, _n2 = s.cont()
    w.row(fg(awr), nk)
    for _ in range(nk):
        p = s.tab1()
        lf = p.l2
        w.row(lf)
        w.vector(p)
        if lf == 1:
            (_, _, _, _, ne), nbt, interp = s.tab2()
            w.row(ne)
            w.interp(nbt, interp)
            for _ in range(ne):
                g = s.tab1()
                w.row(fe(g.c2), len(g.x))
                w.interp(g.nbt, g.interp)
                w.pairs(g.x, g.y)
        elif lf == 5:
            w.vector(s.tab1())          # theta(E)
            w.vector(s.tab1())          # g(x)
        elif lf in (7, 9):
            w.vector(s.tab1())          # theta(E)
        elif lf == 11:
            w.vector(s.tab1())          # a(E)
            w.vector(s.tab1())          # b(E)
        elif lf == 12:
            tm = s.tab1()
            w.row(fe(tm.c1), fe(tm.c2))  # EFL, EFH
            w.vector(tm)
            print("  warning: MF=5 LF=12 (Madland-Nix) written as-is; INDC(NDS)-0758 "
                  "converts it to LF=1 with PREPRO SPECTRA because of Geant4 problems",
                  file=sys.stderr)
        else:
            raise ConversionError(f"MF=5 LF={lf} has no Geant4 representation")


# ---- CrossSection/ --------------------------------------------------------

def build_cross_section(sections):
    if (3, 18) not in sections:
        return None
    t = Section(sections[(3, 18)])
    t.cont()
    xs = t.tab1()
    if not is_linear(xs):
        raise ConversionError("MF=3 MT=18 is not lin-lin; run PREPRO LINEAR first")
    w = G4Writer()
    w.row(0, 0, len(xs.x))
    write_xs_pairs(w, xs)
    return w.text()


# ---- FC/ SC/ TC/ LC/ ------------------------------------------------------

def build_chance(sections, mt, awr):
    if (3, mt) not in sections:
        return None
    if (6, mt) in sections:
        raise ConversionError(f"MT={mt} uses MF=6; convert to MF=4/5 with PREPRO SIXPAK first")
    t = Section(sections[(3, mt)])
    t.cont()
    xs = t.tab1()
    if not is_linear(xs):
        raise ConversionError(f"MF=3 MT={mt} is not lin-lin; run PREPRO LINEAR first")
    w = G4Writer()
    w.row(int(round(xs.c2)), xs.l2, len(xs.x))      # QM [eV], LR
    write_xs_pairs(w, xs)
    # G4ParticleHPFissionBaseFS::Init stops here when there is no further token.
    if (5, mt) in sections:
        w.row(1)
        w.row(4)
        if (4, mt) in sections:
            write_angular(w, sections[(4, mt)], awr)
        else:
            w.row(0, fg(awr), 1)                   # isotropic, LAB
        w.row(3)
        w.row(5)
        write_energy(w, sections[(5, mt)])
    elif (4, mt) in sections:
        print(f"  warning: MT={mt} has MF=4 but no MF=5; final state omitted",
              file=sys.stderr)
    return w.text()


# ---- FS/ ------------------------------------------------------------------

def write_nubar_total(w, sec):
    s = Section(sec)
    _za, awr, _l1, lnu, _n1, _n2 = s.cont()
    w.row(fg(awr), lnu)
    if lnu == 1:
        _, coeff = s.list()
        w.row(len(coeff))
        w.floats(coeff)
    else:
        w.vector(s.tab1())


def write_nubar_delayed(w, sec):
    s = Section(sec)
    _za, awr, ldg, lnu, _n1, _n2 = s.cont()
    if ldg != 0:
        raise ConversionError("MF=1 MT=455 LDG=1 (energy-dependent decay constants) "
                              "has no Geant4 representation")
    _, lambdas = s.list()
    w.row(fg(awr), lnu)
    w.row(len(lambdas))
    w.floats(lambdas)
    if lnu == 2:
        w.vector(s.tab1())
    else:
        _, coeff = s.list()
        if len(coeff) > 1:
            print("  warning: MT=455 polynomial nu-bar reduced to its constant term",
                  file=sys.stderr)
        w.row(fe(coeff[0]))


def write_nubar_prompt(w, sec):
    s = Section(sec)
    _za, awr, _l1, lnu, _n1, _n2 = s.cont()
    w.row(fg(awr), lnu)
    if lnu == 2:
        w.vector(s.tab1())
    else:
        _, coeff = s.list()
        if len(coeff) > 1:
            print("  warning: MT=456 polynomial nu-bar reduced to its constant term",
                  file=sys.stderr)
        w.row(fe(coeff[0]))


def write_energy_release(w, sec):
    """MF=1 MT=458 -> G4ParticleHPFissionERelease::Init: a dummy, then the
    constant terms of EFR ENP END EGP EGD EB ENU ER ET (uncertainties dropped)."""
    s = Section(sec)
    s.cont()
    _, v = s.list()
    e = v[0:18:2]
    w.row(0)
    w.row(*(fg(x) for x in e[0:4]))
    w.row(*(fg(x) for x in e[4:8]))
    w.row(fg(e[8]))


def write_photon_mult(w, sec):
    """MF=12 -> G4ParticleHPPhotonDist::InitMean. Returns the list of LF flags."""
    s = Section(sec)
    _za, awr, lo, _lg, nk, _n2 = s.cont()
    if lo != 1:
        raise ConversionError(f"MF=12 MT=18 LO={lo} not supported (fission uses LO=1)")
    if nk > 1:
        s.tab1()                    # total yield: Geant4 does not read it
    w.row(lo, fg(awr))
    w.row(nk)
    flags = []
    for _ in range(nk):
        y = s.tab1()
        flags.append(y.l2)
        w.row(y.l2, fg(y.c1))      # LF, Eg
        w.vector(y)
    return flags


def write_photon_angular(w, sec):
    """MF=14 -> G4ParticleHPPhotonDist::InitAngular."""
    s = Section(sec)
    _za, _awr, li, ltt, nk, ni = s.cont()
    w.row(li)
    if li == 1:
        return
    w.row(ltt, nk, ni)
    for _ in range(ni):
        eg, es, *_ = s.cont()
        w.row(fe(eg), fe(es))
    for _ in range(nk - ni):
        (eg, es, _, _, ne), nbt, interp = s.tab2()
        w.row(fe(eg), fe(es), ne)
        if ltt == 1:
            w.interp(nbt, interp)
            for _ in range(ne):
                (_c1, e, _l1, _l2, _nl, _n2), coeff = s.list()
                w.row(fe(e), len(coeff))
                w.floats(coeff)
        elif ltt == 2:
            # G4ParticleHPAngularP carries no incident-energy interpolation table.
            for _ in range(ne):
                f = s.tab1()
                w.row(fe(f.c2), len(f.x))
                w.interp(f.nbt, f.interp)
                w.pairs(f.x, f.y)
        else:
            raise ConversionError(f"MF=14 LTT={ltt} not supported")


def write_photon_energy(w, sec):
    """MF=15 -> G4ParticleHPPhotonDist::InitEnergies."""
    s = Section(sec)
    _za, _awr, _l1, _l2, nc, _n2 = s.cont()
    w.row(nc)
    for _ in range(nc):
        p = s.tab1()
        if p.l2 != 1:
            raise ConversionError(f"MF=15 LF={p.l2} not supported")
        w.row(p.l2)
        w.vector(p)
        (_, _, _, _, ne), nbt, interp = s.tab2()
        w.row(ne)
        w.interp(nbt, interp)
        for _ in range(ne):
            g = s.tab1()
            w.row(fe(g.c2), len(g.x))
            w.interp(g.nbt, g.interp)
            w.pairs(g.x, g.y)


def build_fs(sections):
    if (3, 18) not in sections:
        return None
    if (6, 18) in sections and (5, 18) not in sections:
        raise ConversionError("MT=18 uses MF=6; convert to MF=4/5/12/14/15 with "
                              "PREPRO SIXPAK first (INDC(NDS)-0758, 232Th/231Pa/233Pa)")
    awr = Section(sections[(3, 18)]).cont()[1]
    flags = []

    def photon_mult(w, sec):
        flags.extend(write_photon_mult(w, sec))

    def photon_energy(w, sec):
        # InitEnergies reads nothing unless some photon is continuous (LF=1);
        # writing the block anyway would desynchronise the stream.
        if 1 in flags:
            write_photon_energy(w, sec)
            return True
        return False

    # (infoType, dataType) tags as read by G4ParticleHPFSFissionFS::Init.
    blocks = {
        (1, 452): ((2, 1), write_nubar_total),
        (1, 455): ((3, 1), write_nubar_delayed),
        (1, 456): ((4, 1), write_nubar_prompt),
        (1, 458): ((5, 1), write_energy_release),
        (4, 18): ((1, 4), lambda w, sec: write_angular(w, sec, awr)),
        (5, 18): ((1, 5), write_energy),
        (5, 455): ((3, 5), write_energy),
        (12, 18): ((1, 12), photon_mult),
        (14, 18): ((1, 14), write_photon_angular),
        (15, 18): ((1, 15), photon_energy),
    }
    # Blocks follow tape order, like the IAEA files (sections SIXPAK added at the
    # end of the tape come last); Geant4 itself only needs MF=12 before 14 and 15.
    w = G4Writer()
    for key in sections:
        if key not in blocks:
            continue
        if key in ((14, 18), (15, 18)) and (12, 18) not in sections:
            continue
        (info, dtype), fn = blocks[key]
        sub = G4Writer()
        sub.row(info)
        sub.row(dtype)
        if fn(sub, sections[key]) is not False:
            w.lines += sub.lines
    return w.text()


# --------------------------------------------------------------------------

def target_name(sections):
    head = Section(sections[(1, 451)])
    za, _awr, lrp, _lfi, _nlib, _nmod = head.cont()
    _elis, _sta, _lis, liso, _n1, _nfor = head.cont()
    z, a = divmod(int(round(za)), 1000)
    m = f"m{liso}" if liso > 0 else ""
    return z, a, liso, lrp, f"{z}_{a}{m}_{ELEMENTS[z - 1]}"


def convert(tape, outdir, only=None, compress=False):
    _mat, sections = read_tape(tape)
    z, a, liso, lrp, name = target_name(sections)
    awr = Section(sections[(1, 451)]).cont()[1]
    dirs = only or ALL_DIRS

    if lrp == 1 and "CrossSection" in dirs:
        print(f"  warning: LRP=1 - resonance parameters in MF=2 not yet added to MF=3; "
              "run PREPRO RECENT or the cross section will miss the resonance region",
              file=sys.stderr)

    written = []
    for d in dirs:
        if d == "CrossSection":
            text = build_cross_section(sections)
        elif d == "FS":
            text = build_fs(sections)
        elif d in CHANCE_DIRS:
            text = build_chance(sections, CHANCE_DIRS[d], awr)
        else:
            raise ConversionError(f"unknown directory {d}")
        if text is None:
            continue
        os.makedirs(os.path.join(outdir, d), exist_ok=True)
        path = os.path.join(outdir, d, name)
        with open(path, "w") as fh:
            fh.write(text)
        if compress:
            with open(path + ".z", "wb") as fh:
                fh.write(zlib.compress(text.encode(), 9))
        written.append(path)
    return name, written


def main():
    global DIGITS
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("endf", nargs="+", help="ENDF-6 neutron tape(s), PREPRO-processed")
    p.add_argument("outdir", help="output Fission/ directory")
    p.add_argument("--only", help="comma-separated subset of " + ",".join(ALL_DIRS))
    p.add_argument("--compress", action="store_true",
                   help="also write <file>.z (raw zlib, as G4NDL ships)")
    p.add_argument("--digits", type=int, default=7,
                   help="significant digits for reals (default 7, as the IAEA files)")
    args = p.parse_args()

    DIGITS = args.digits

    only = args.only.split(",") if args.only else None
    failed = 0
    for tape in args.endf:
        try:
            name, written = convert(tape, args.outdir, only, args.compress)
            dirs = " ".join(os.path.basename(os.path.dirname(w)) for w in written)
            print(f"{tape}: {name} -> {dirs or '(no fission data)'}", file=sys.stderr)
        except (ConversionError, IndexError, ValueError) as e:
            failed += 1
            print(f"{tape}: ERROR {e}", file=sys.stderr)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
