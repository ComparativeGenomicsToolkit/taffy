#!/bin/bash
#
# Make a universal per-species (or total) depth bigWig (local, single process)
# ===========================================================================
# Builds a 64-bit per-species coverage bigWig on the UNIVERSAL-COLUMN axis from a
# `cactus-hal2maf --universal` MAF/TAF, using taffy's in-tree 64-bit vector bigWig
# writer -- one component per leaf genome, a single chrom uni0 of length T, plus
# an <out>.names sidecar (leaf names in component order).  No more wigToBigWig and
# no more 2e9 uni0..uniK chunking (the 64-bit axis holds [0,T) directly).
#
#   1. taffy index -u                       build <input>.tui if missing (T / coords)
#   2. taffy depth --bin N --perSpecies --bigwig OUT.bw
#                                           per-species covered-count vector bigWig
#
# Query per leaf with:
#   taffy lift --bigwig OUT.bw -g <leaf> -r <leaf.seq:start-end> -B <bin>
# (maps the leaf window to columns via the .tui; emits per-species coverage).
#
# Single process: fine up to ~a fish-subtree alignment.  For the full 577-way,
# shard it over SLURM instead (uni_depth_bigwig_slurm.sh).
#
# Requirements: taffy (`depth` + `lift`).  Env override: TAFFY.

set -euo pipefail

INPUT=""; OUT=""; BIN=1000; THREADS=8; TUI_THREADS=""; TMPDIR_OVR=""; SCALAR=0
TAFFY="${TAFFY:-taffy}"

usage() {
  cat >&2 <<'U'
Usage: uni_depth_bigwig_local.sh -i UNI.maf.gz -o OUT.bw [options]
  -i  universal MAF/TAF (cactus-hal2maf --universal); its .tui is built if missing
  -o  output bigWig (universal-column uni0 axis; also writes OUT.bw.names)
  -b  bin width in bp (default 1000)
  -T  threads for the depth pass (CPU-bound).  Default 8.
  -j  threads for the .tui build (MEMORY-bound: each worker holds one genome's
      runs[], up to tens of GB at vertebrate scale).  Set BELOW -T on small-RAM
      boxes to avoid OOM.  Default: same as -T.
  -d  tmp/spill dir for the .tui build (put it on a ROOMY fs for big inputs).
      Default: the output directory.
  -S  scalar total-depth bigWig (mean leaves/bin) instead of the per-species
      vector.  Default: per-species (one component per leaf).
Env: TAFFY overrides the binary.
U
  exit "${1:-1}"
}

while getopts "i:o:b:T:j:d:Sh" opt; do
  case "$opt" in
    i) INPUT=$OPTARG ;;
    o) OUT=$OPTARG ;;
    b) BIN=$OPTARG ;;
    T) THREADS=$OPTARG ;;
    j) TUI_THREADS=$OPTARG ;;
    d) TMPDIR_OVR=$OPTARG ;;
    S) SCALAR=1 ;;
    h) usage 0 ;;
    *) usage ;;
  esac
done

[ -n "$INPUT" ] && [ -n "$OUT" ] || usage
[ -e "$INPUT" ] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
command -v "$TAFFY" >/dev/null 2>&1 || { echo "ERROR: taffy not found (set \$TAFFY)" >&2; exit 1; }
(( BIN > 0 )) || { echo "ERROR: -b BIN ($BIN) must be > 0" >&2; exit 1; }

OUTDIR=$(dirname "$OUT"); mkdir -p "$OUTDIR"

# --- 1. universal column index (.tui; needed for T / column coords) ----------
[ -n "$TUI_THREADS" ] || TUI_THREADS=$THREADS
if [ -e "$INPUT.tui" ]; then
  echo "[1/2] using existing $INPUT.tui" >&2
else
  echo "[1/2] building $INPUT.tui (taffy index -u -T $TUI_THREADS) ..." >&2
  "$TAFFY" index -i "$INPUT" -u -T "$TUI_THREADS" ${TMPDIR_OVR:+-d "$TMPDIR_OVR"}
fi

# --- 2. binned depth -> 64-bit bigWig directly (no wigToBigWig) ---------------
if (( SCALAR )); then
  echo "[2/2] taffy depth --bin $BIN --bigwig (scalar total-depth) -T $THREADS ..." >&2
  "$TAFFY" depth -i "$INPUT" --bin "$BIN" --minLeaves 1 -T "$THREADS" --bigwig "$OUT"
else
  echo "[2/2] taffy depth --bin $BIN --perSpecies --bigwig -T $THREADS ..." >&2
  "$TAFFY" depth -i "$INPUT" --bin "$BIN" --minLeaves 1 -T "$THREADS" --perSpecies --bigwig "$OUT"
fi

echo "" >&2
echo "done: $OUT ($(du -h "$OUT" 2>/dev/null | cut -f1))" >&2
(( SCALAR )) || echo "  per-species: $(wc -l < "$OUT.names" 2>/dev/null) leaf components in $OUT.names" >&2
echo "  on the universal-column uni0 axis (64-bit)." >&2
echo "  query per leaf: taffy lift --bigwig $OUT -g <leaf> -r <leaf.seq:start-end> -B <bin>" >&2
