#!/bin/bash
#
# Make a universal-depth bigWig (local, single process)
# =====================================================
# Builds the value source for the `taffy lift --bigwig` query shim: a depth
# bigWig keyed on the INTEGER universal column (chrom uni<chunk>), from a
# `cactus-hal2maf --universal` MAF/TAF.  Three steps:
#
#   1. taffy index -u        build <input>.tui if missing (universal column index)
#   2. taffy gerp --universal --depthOnly --bin N
#                            mean leaf depth per N-column bin -> bedGraph
#                            (uni<chunk> axis, monotone / sorted by construction)
#   3. wigToBigWig           -> <out>  (+ a <out>.uni.sizes chrom-sizes file)
#
# Query the result with:
#   taffy lift -i <input> --bigwig <out> -g GENOME -r GENOME.chrom:A-B -B BINBP
#
# Single process: fine up to ~a fish-subtree alignment.  For the full 577-way
# the gerp pass won't finish in one process -- shard it over SLURM instead
# (uni_depth_bigwig_slurm.sh; same column-range model as gerp_shard_slurm.sh).
#
# Requirements: taffy (with `gerp --bin` -- commit 8f9c5d7+), wigToBigWig (UCSC).
# Env overrides: TAFFY, WIGTOBIGWIG.

set -euo pipefail

INPUT=""; OUT=""; BIN=1000; THREADS=8; TUI_THREADS=""; TMPDIR_OVR=""; KEEP_BG=0
TAFFY="${TAFFY:-taffy}"
WIGTOBIGWIG="${WIGTOBIGWIG:-wigToBigWig}"

usage() {
  cat >&2 <<'U'
Usage: uni_depth_bigwig_local.sh -i UNI.maf.gz -o OUT.bw [options]
  -i  universal MAF/TAF (cactus-hal2maf --universal); its .tui is built if missing
  -o  output bigWig
  -b  bin width in universal columns; must divide 4e9 (default 1000)
  -T  threads for the gerp pass (CPU-bound).  Default 8.
  -j  threads for the .tui build (MEMORY-bound: each worker holds one genome's
      runs[], up to tens of GB at vertebrate scale).  Set BELOW -T on small-RAM
      boxes to avoid OOM.  Default: same as -T.
  -d  tmp/spill dir for the .tui build + the intermediate bedGraph (put it on a
      ROOMY fs for big inputs).  Default: the output directory.
  -k  keep the intermediate bedGraph (default: delete after wigToBigWig)
Env: TAFFY, WIGTOBIGWIG override the binaries.
U
  exit "${1:-1}"
}

while getopts "i:o:b:T:j:d:kh" opt; do
  case "$opt" in
    i) INPUT=$OPTARG ;;
    o) OUT=$OPTARG ;;
    b) BIN=$OPTARG ;;
    T) THREADS=$OPTARG ;;
    j) TUI_THREADS=$OPTARG ;;
    d) TMPDIR_OVR=$OPTARG ;;
    k) KEEP_BG=1 ;;
    h) usage 0 ;;
    *) usage ;;
  esac
done

[ -n "$INPUT" ] && [ -n "$OUT" ] || usage
[ -e "$INPUT" ] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
command -v "$TAFFY" >/dev/null 2>&1 || { echo "ERROR: taffy not found (set \$TAFFY)" >&2; exit 1; }
command -v "$WIGTOBIGWIG" >/dev/null 2>&1 || { echo "ERROR: wigToBigWig not found (set \$WIGTOBIGWIG)" >&2; exit 1; }
# gerp --bin requires N to divide the 4e9 universal chunk (so bins never straddle
# a uni<chunk> boundary); the same value drives the uni<chunk> sizes below.
if (( 4000000000 % BIN != 0 )); then
  echo "ERROR: -b BIN ($BIN) must divide 4000000000 (e.g. 1000, 10000)" >&2; exit 1
fi

OUTDIR=$(dirname "$OUT"); mkdir -p "$OUTDIR"
WORK="${TMPDIR_OVR:-$OUTDIR}"; mkdir -p "$WORK"
BASE=$(basename "$OUT"); BASE=${BASE%.bw}; BASE=${BASE%.bigWig}
BG="$WORK/$BASE.bg"
SIZES="$OUTDIR/$BASE.uni.sizes"

# --- 1. universal column index (.tui) ---------------------------------------
[ -n "$TUI_THREADS" ] || TUI_THREADS=$THREADS
if [ -e "$INPUT.tui" ]; then
  echo "[1/3] using existing $INPUT.tui" >&2
else
  echo "[1/3] building $INPUT.tui (taffy index -u -T $TUI_THREADS) ..." >&2
  "$TAFFY" index -i "$INPUT" -u -T "$TUI_THREADS" ${TMPDIR_OVR:+-d "$TMPDIR_OVR"}
fi

T=$("$TAFFY" stats -i "$INPUT" -u)
echo "[*]   T=$T universal columns; bin=$BIN  (~$((T/BIN)) bins)" >&2

# --- 2. binned universal depth bedGraph (monotone -> needs no external sort) -
echo "[2/3] taffy gerp --universal --depthOnly --bin $BIN -T $THREADS ..." >&2
"$TAFFY" gerp -i "$INPUT" --universal --depthOnly --bin "$BIN" --minLeaves 1 \
  -T "$THREADS" -o /dev/null -D "$BG"

# --- 3. uni<chunk> chrom sizes + bigWig -------------------------------------
awk -v T="$T" -v C=4000000000 \
  'BEGIN{for(c=0;c*C<T;c++){s=((c+1)*C<=T)?C:T-c*C; print "uni"c"\t"s}}' > "$SIZES"
echo "[3/3] wigToBigWig -> $OUT ..." >&2
"$WIGTOBIGWIG" "$BG" "$SIZES" "$OUT"
[ "$KEEP_BG" = 1 ] || rm -f "$BG"

echo "" >&2
echo "done: $OUT ($(du -h "$OUT" 2>/dev/null | cut -f1)); sizes: $SIZES" >&2
echo "query: $TAFFY lift -i $INPUT --bigwig $OUT -g GENOME -r GENOME.chrom:A-B -B BINBP" >&2
