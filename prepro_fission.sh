#!/bin/sh
# Pre-process one ENDF-6 neutron tape with PREPRO the way INDC(NDS)-0758 §2
# describes for the IAEA G4NDL libraries, ready for endf2g4ndl_fission.py:
#
#   LINEAR  0.1%   MF=3 to lin-lin
#   RECENT  0.1%   add resonance contribution (MF=2) to MF=3
#   LINEAR  0.1%   thin the reconstructed MF=3 (INDC says 1%, but the IAEA
#                  JEFF-3.3 files are reproduced exactly only with 0.1%)
#   FIXUP          tabulate polynomial nu-bar, insert a 0.0253 eV point in every table
#   DICTIN         rebuild the MF=1/MT=451 index
#   SIXPAK         only if fission is given in MF=6: writes MF=4/5/12/14/15,
#                  appended to the output tape
#
# Usage: prepro_fission.sh <endf-tape> <output-tape>
# Tolerances can be overridden with LIN1_ERR, RECENT_ERR, LIN2_ERR (fractions).

set -e
[ $# -eq 2 ] || { echo "usage: $0 <endf-tape> <output-tape>" >&2; exit 2; }

here=$(cd "$(dirname "$0")" && pwd)
bin=${PREPRO_BIN:-$here/../bin}
in=$(realpath "$1")
out=$(realpath -m "$2")
LIN1_ERR=${LIN1_ERR:-1.0000E-03}
RECENT_ERR=${RECENT_ERR:-1.0000E-03}
LIN2_ERR=${LIN2_ERR:-1.0000E-03}

work=$(mktemp -d "${TMPDIR:-/tmp}/prepro_fission.XXXXXX")
trap 'rm -rf "$work"' EXIT
cd "$work"
cp "$in" ENDFB.IN

printf '          0          0 1.0000E-10          0\nENDFB.IN\nLINEAR1.OUT\n\n %s\n\n' \
    "0.0000E+00 $LIN1_ERR" > LINEAR.INP
"$bin/linear" > /dev/null

printf '          0 1.0000E-10          1          0          0          0\nLINEAR1.OUT\nRECENT.OUT\n\n %s\n\n' \
    "0.0000E+00 $RECENT_ERR" > RECENT.INP
"$bin/recent" > /dev/null

printf '          0          0 1.0000E-10          0\nRECENT.OUT\nLINEAR2.OUT\n\n %s\n\n' \
    "0.0000E+00 $LIN2_ERR" > LINEAR.INP
"$bin/linear" > /dev/null

printf '10002111110101\nLINEAR2.OUT\nFIXUP.OUT\n 2.53000- 2   0  0\n\n' > FIXUP.INP
"$bin/fixup" > /dev/null

printf 'FIXUP.OUT\nDICTIN.OUT\n' > DICTIN.INP
"$bin/dictin" > /dev/null

cp DICTIN.OUT "$out"

# MF=6 for any fission channel -> SIXPAK
if awk 'substr($0,71,2)+0==6 { mt=substr($0,73,3)+0; if (mt==18||mt==19||mt==20||mt==21||mt==38) f=1 }
        END { exit !f }' DICTIN.OUT; then
    printf 'DICTIN.OUT\n\n' > SIXPAK.INP
    "$bin/sixpak" > /dev/null
    for f in ENDFB.MF4 ENDFB.MF5 ENDFB.M12 ENDFB.M14 ENDFB.M15; do
        [ -s "$f" ] && cat "$f" >> "$out"
    done
    echo "$1: fission given in MF=6, SIXPAK output appended" >&2
fi
