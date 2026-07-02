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


# Computes the generalised Hurwitz class number. For F = Q(sqrt D) it computes
# every Hurwitz class number H_F(4 p - t^2) that the trace formula needs for
# the totally positive prime element p with N(p) <= A.
#
# The H_F values do not depend on the weight (k1,k2) -- only the final trace
# assembly does -- so run this once per (D,A), then run HTtrace.sh for as many
# weights as you like, all reusing the table produced here.
#
# Usage:  sh HTclassno.sh <D> <A> [ncores]
# Phase 1 streams arguments to disk with NO per-core dedup map (emit_stream),
# so each worker uses ~constant RAM regardless of A and all cores run; the
# dedup is done by an external `sort -u` with a bounded buffer (SORTBUF) and
# its temp files on this filesystem (SORTTMP), not /tmp.  Trades RAM for
# transient disk (~tens of GB at A=10^7).  Phase 2 (the class numbers) uses
# all cores; its only memory is hurwitz's per-worker field cache (modest).
# Output: hf_D<D>_A<A>.txt   ("wa wb num den" per needed argument; UNCONDITIONAL)
cd "$(dirname "$0")"

D=${1:-5}; A=${2:-10000000}; p=${3:-16}
PARISIZE=256000000
PARISIZEMAX=2000000000
SORTBUF=8G                 # external-sort memory cap (one shared sort, NOT per-core)
SORTTMP=.                  # sort scratch on this (big) filesystem, not /tmp
ARGS="args_D${D}_A${A}.txt"
HF="hf_D${D}_A${A}.txt"

run() { echo "read(\"HTprime.gp\"); $1" | gp -q -s "$PARISIZE" --default parisizemax="$PARISIZEMAX"; }

echo "=== [D=$D, A=$A] phase 1: stream H_F arguments to disk ($p workers, no per-core map) ==="
for z in `seq 1 $p`; do run "emit_stream($D,$A,$z,$p,\"args.tmp.$z\");" & done
wait
echo "    deduping on disk (sort -u, $SORTBUF buffer, temp in $SORTTMP) ..."
sort -u --parallel=$p -S "$SORTBUF" -T "$SORTTMP" args.tmp.* > "$ARGS"; rm -f args.tmp.*
echo "    distinct arguments: `wc -l < "$ARGS"`"

echo "=== [D=$D, A=$A] phase 2: compute H_F (unconditional; weight-independent) ==="
for z in `seq 1 $p`; do run "compute_HF($D,\"$ARGS\",$z,$p,\"hf.tmp.$z\");" & done
wait
sort -n -k1,1 -k2,2 --parallel=$p hf.tmp.* > "$HF"; rm -f hf.tmp.*
echo "    class numbers written to $HF: `wc -l < "$HF"`"
echo "    reuse for any weight:   sh HTtrace.sh $D <k1> <k2> $A"
