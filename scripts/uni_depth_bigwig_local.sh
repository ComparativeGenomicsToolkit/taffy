#!/bin/bash
#
# Make a universal-depth bigWig (local, single process)
# =====================================================
# Builds a mean-leaf-depth bigWig on the UNIVERSAL-COLUMN axis from a
# `cactus-hal2maf --universal` MAF/TAF.  ONE `taffy depth` command produces both
# the bedGraph and the chrom-sizes (the axis is the raw column [0,T) 2e9-chunked
# into chroms uni0,uni1,...); then wigToBigWig.  Three steps:
#
#   1. taffy index -u   build <input>.tui if missing (needed for T / column coords)
#   2. taffy depth --bin N --depth BG --sizes SIZES
#                       mean leaf depth per N-bp bin on the chunked uni axis,
#                       already monotone/sorted (no downstream sort/merge), plus
#                       the uni0..uniK chrom-sizes file
#   3. wigToBigWig      BG SIZES -> <out>
#
# The bigWig is on the universal-column axis, so query it per leaf with
# `taffy lift --bigwig BW -g <leaf> -r <leaf.seq:start-end>` (which maps the
# leaf window to columns via the .tui and reads the chunked uni axis).
#
# Single process: fine up to ~a fish-subtree alignment.  For the full 577-way
# shard it over SLURM instead (uni_depth_bigwig_slurm.sh).
#
# Requirements: taffy (`depth` subcommand), wigToBigWig (UCSC).
# Env overrides: TAFFY, WIGTOBIGWIG.

set -euo pipefail

INPUT=""; OUT=""; BIN=1000; THREADS=8; TUI_THREADS=""; TMPDIR_OVR=""; KEEP_BG=0
TAFFY="${TAFFY:-taffy}"
WIGTOBIGWIG="${WIGTOBIGWIG:-wigToBigWig}"

usage() {
  cat >&2 <<'U'
Usage: uni_depth_bigwig_local.sh -i UNI.maf.gz -o OUT.bw [options]
  -i  universal MAF/TAF (cactus-hal2maf --universal); its .tui is built if missing
  -o  output bigWig (on the universal-column uni0,uni1,... axis)
  -b  bin width in bp (default 1000; must divide 2e9, the uni-axis chunk size)
  -T  threads for the depth pass (CPU-bound).  Default 8.
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
# --bin must divide 2e9 (the uni-axis chunk size) so no bin straddles a chunk.
(( BIN > 0 )) || { echo "ERROR: -b BIN ($BIN) must be > 0" >&2; exit 1; }
(( 2000000000 % BIN == 0 )) || { echo "ERROR: -b BIN ($BIN) must divide 2000000000 (the uni-axis chunk size)" >&2; exit 1; }

OUTDIR=$(dirname "$OUT"); mkdir -p "$OUTDIR"
WORK="${TMPDIR_OVR:-$OUTDIR}"; mkdir -p "$WORK"
BASE=$(basename "$OUT"); BASE=${BASE%.bw}; BASE=${BASE%.bigWig}
BG="$WORK/$BASE.bg"                  # chunked uni-axis mean-depth bedGraph (sorted)
SIZES="$OUTDIR/$BASE.sizes"          # uni0..uniK chrom sizes (from taffy depth)

# --- 1. universal column index (.tui; needed for T / column coords) ----------
[ -n "$TUI_THREADS" ] || TUI_THREADS=$THREADS
if [ -e "$INPUT.tui" ]; then
  echo "[1/3] using existing $INPUT.tui" >&2
else
  echo "[1/3] building $INPUT.tui (taffy index -u -T $TUI_THREADS) ..." >&2
  "$TAFFY" index -i "$INPUT" -u -T "$TUI_THREADS" ${TMPDIR_OVR:+-d "$TMPDIR_OVR"}
fi

# --- 2. binned depth on the chunked uni axis + chrom sizes (one command) ------
echo "[2/3] taffy depth --bin $BIN --depth --sizes -T $THREADS ..." >&2
"$TAFFY" depth -i "$INPUT" --bin "$BIN" --minLeaves 1 \
  -T "$THREADS" --depth "$BG" --sizes "$SIZES"

# --- 3. bigWig (bedGraph already monotone/sorted; no sort/merge) --------------
echo "[3/3] wigToBigWig -> $OUT ..." >&2
"$WIGTOBIGWIG" "$BG" "$SIZES" "$OUT"
[ "$KEEP_BG" = 1 ] || rm -f "$BG"

echo "" >&2
echo "done: $OUT ($(du -h "$OUT" 2>/dev/null | cut -f1)); sizes: $SIZES" >&2
echo "  on the universal-column uni0,uni1,... axis." >&2
echo "  query per leaf: taffy lift --bigwig $OUT -g <leaf> -r <leaf.seq:start-end>" >&2
