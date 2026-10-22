#!/usr/bin/env python3
"""Compare G4NDL files token by token, or two whole Fission/ trees.

G4NDL files are read with operator>>, so whitespace is meaningless and two
files are equivalent when their token streams agree. Integers must match
exactly; reals to a relative tolerance (default 1e-5, i.e. the 7 significant
digits G4NDL carries). Accepts .z (raw zlib) files.

A cross-section table (CrossSection/, and the head of FC/ SC/ TC/ LC/) is
also compared as a function: lin-lin interpolation of A on the energy grid
of B. Pointwise grids from different PREPRO versions or tolerances never
match token for token, but the functions should agree to within the
reconstruction tolerance.

Usage:
    g4ndl_compare.py A B                  two files
    g4ndl_compare.py --tree DIR_A DIR_B   every file under DIR_B/<sub>/ that DIR_A also has
"""

import argparse
import bisect
import math
import os
import sys
import zlib

XS_DIRS = ("CrossSection", "FC", "SC", "TC", "LC")


def tokens(path):
    if not os.path.exists(path) and os.path.exists(path + ".z"):
        path += ".z"
    with open(path, "rb") as fh:
        raw = fh.read()
    if path.endswith(".z"):
        raw = zlib.decompress(raw)
    return raw.decode().split()


def same(a, b, rtol):
    if a == b:
        return True
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False
    return abs(x - y) <= rtol * max(abs(x), abs(y)) + 1e-300


def token_diff(ta, tb, rtol):
    return [i for i in range(min(len(ta), len(tb))) if not same(ta[i], tb[i], rtol)]


def xs_table(tok):
    """(header, energies, values, index of first token after the table)."""
    n = int(tok[2])
    v = [float(t) for t in tok[3:3 + 2 * n]]
    return tok[:2], v[0::2], v[1::2], 3 + 2 * n


def interp(xs, ys, x):
    i = bisect.bisect_left(xs, x)
    if i < len(xs) and xs[i] == x:
        return ys[i]
    if i == 0 or i == len(xs):
        return 0.0
    x0, x1, y0, y1 = xs[i - 1], xs[i], ys[i - 1], ys[i]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def integral(xs, ys, lo, hi):
    """Integral of the lin-lin function over [lo, hi]."""
    pts = [lo] + [x for x in xs if lo < x < hi] + [hi]
    vals = [interp(xs, ys, x) for x in pts]
    return sum(0.5 * (vals[i] + vals[i + 1]) * (pts[i + 1] - pts[i]) for i in range(len(pts) - 1))


def xs_diff(ea, sa, eb, sb):
    """Largest relative difference between the integrals of A and B over any
    energy decade, and the fraction of B's points where A differs by more than 1%.

    The pointwise figure is only indicative: G4NDL prints energies to 7 digits,
    so points in narrow resonances collapse onto the same printed energy in both
    files. Decade integrals are insensitive to that and to the choice of grid.
    """
    lo, hi = max(ea[0], eb[0]), min(ea[-1], eb[-1])
    worst, e = 0.0, 10 ** math.floor(math.log10(lo))
    while e < hi:
        a0, a1 = max(e, lo), min(10 * e, hi)
        if a1 > a0:
            ia, ib = integral(ea, sa, a0, a1), integral(eb, sb, a0, a1)
            if max(ia, ib) > 0:
                worst = max(worst, abs(ia - ib) / max(abs(ia), abs(ib)))
        e *= 10
    over = sum(abs(interp(ea, sa, x) - y) > 0.01 * abs(y) for x, y in zip(eb, sb))
    return worst, over / len(eb)


def compare_file(pa, pb, sub, rtol):
    ta, tb = tokens(pa), tokens(pb)
    if sub in XS_DIRS:
        ha, ea, sa, ra = xs_table(ta)
        hb, eb, sb, rb = xs_table(tb)
        worst, frac = xs_diff(ea, sa, eb, sb)
        rest = token_diff(ta[ra:], tb[rb:], rtol)
        head_ok = all(same(a, b, rtol) for a, b in zip(ha, hb))
        ok = head_ok and not rest and len(ta) - ra == len(tb) - rb
        return (f"xs {len(ea)}/{len(eb)} pts, decade-integral diff {worst:.2e}, points >1% off: {frac:.2%}; "
                f"header {'ok' if head_ok else ' '.join(ha) + ' vs ' + ' '.join(hb)}; "
                f"final state {len(ta) - ra}/{len(tb) - rb} tokens, {len(rest)} mismatches"), ok
    bad = token_diff(ta, tb, rtol)
    ok = not bad and len(ta) == len(tb)
    msg = f"{len(ta)}/{len(tb)} tokens, {len(bad)} mismatches"
    if bad:
        i = bad[0]
        msg += f"; first @{i}: [{' '.join(ta[max(0, i - 3):i + 4])}] vs [{' '.join(tb[max(0, i - 3):i + 4])}]"
    return msg, ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--tree", action="store_true")
    p.add_argument("--rtol", type=float, default=1e-5)
    p.add_argument("--quiet", action="store_true", help="only list files that differ")
    args = p.parse_args()

    if not args.tree:
        sub = os.path.basename(os.path.dirname(os.path.abspath(args.b)))
        msg, ok = compare_file(args.a, args.b, sub, args.rtol)
        print(msg)
        sys.exit(0 if ok else 1)

    total = failed = missing = 0
    for sub in sorted(os.listdir(args.b)):
        db = os.path.join(args.b, sub)
        if not os.path.isdir(db) or sub == "FF":
            continue
        for name in sorted(os.listdir(db)):
            base = name[:-2] if name.endswith(".z") else name
            pa = os.path.join(args.a, sub, base)
            if not (os.path.exists(pa) or os.path.exists(pa + ".z")):
                missing += 1
                print(f"MISSING {sub}/{base}")
                continue
            total += 1
            msg, ok = compare_file(pa, os.path.join(db, name), sub, args.rtol)
            failed += not ok
            if not ok or not args.quiet:
                print(f"{'ok  ' if ok else 'DIFF'} {sub}/{base}: {msg}")
    print(f"compared {total}, differing {failed}, missing {missing}")


if __name__ == "__main__":
    main()
