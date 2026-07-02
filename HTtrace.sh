#!/bin/sh

# Copyright 2026 Seiji Kuga, Andrei Seymour-Howell and Satoshi Wakatsuki

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Trace stage (per weight).  Assembles tr T_p for every totally positive prime element p with
# N(p) <= A in weight (k1,k2), reusing the weight-independent class-number
# table produced by HTclassno.sh.
#
# Usage:  sh HTtrace.sh <D> <k1> <k2> <A> [ncores] [hf-table]
#   requires the class-number table; default hf_D<D>_A<A>.txt (from HTclassno.sh).
#   A 6th arg overrides that table path.
# Handles any weight: parallel (k1=k2) and non-parallel (k1!=k2) alike.
# Output: traces_<D>_<k1><k2>_A<A>.txt   ("Np u v e1 e2" per prime, by norm;
#   e1,e2 are the two Galois embeddings of a_p, with e1=e2 when k1=k2)
cd "$(dirname "$0")"

D=${1:-5}; k1=${2:-6}; k2=${3:-6}; A=${4:-10000}; p=${5:-16}
PARISIZE=256000000
PARISIZEMAX=2000000000
HF=${6:-"hf_D${D}_A${A}.txt"}
OUT="traces_${D}_${k1}${k2}_A${A}.txt"

if [ ! -f "$HF" ]; then
    echo "ERROR: $HF not found -- run the class-number stage first:" >&2
    echo "       sh HTclassno.sh $D $A" >&2
    exit 1
fi

NLINES=`wc -l < "$HF"`           # table size, so each worker pre-sizes its arrays

run() { echo "read(\"HTprime.gp\"); $1" | gp -q -s "$PARISIZE" --default parisizemax="$PARISIZEMAX"; }

echo "=== [D=$D, weight ($k1,$k2), A=$A] assemble tr T_p from $HF ($NLINES class numbers, all $p cores) ==="
for z in `seq 1 $p`; do run "assemble($D,$k1,$k2,$A,\"$HF\",$NLINES,$z,$p,\"tr.tmp.$z\");" & done
wait
sort -n -k1,1 -k2,2 --parallel=$p tr.tmp.* > "$OUT"; rm -f tr.tmp.*
echo "    traces written to $OUT: `wc -l < "$OUT"`"
