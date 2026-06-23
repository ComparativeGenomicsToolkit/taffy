#!/bin/bash
#
# Make a universal-depth bigWig (local, single process)
# =====================================================
# Builds a mean-leaf-depth bigWig on NAMED row-0 (ancestor) coordinates from a
# `cactus-hal2maf --universal` MAF/TAF.  Four steps:
#
#   1. taffy index -u   build <input>.tui if missing (needed for the chrom sizes)
#   2. taffy depth --coords named --bin N --depth
#                       mean leaf depth per N-bp bin, keyed on the row-0 seq, as
#                       per-(name,bin) records `name start end sum cnt` (UNSORTED:
#                       the row-0 ancestor seq recurs through the column order)
#   3. sort + merge     sort -k1,1 -k2,2n then sum (sum,cnt) per (name,bin) ->
#                       a clean 4-col bedGraph of mean depth
#   4. wigToBigWig      -> <out>  (+ a <out>.sizes from `taffy stats -s`)
#
# NOTE: `taffy lift --bigwig` currently queries the older integer-column axis;
# the named-coord query is a pending follow-up.  This builds the named bigWig
# (browsable directly on ancestor coords) in the meantime.
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
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat >&2 <<'U'
Usage: uni_depth_bigwig_local.sh -i UNI.maf.gz -o OUT.bw [options]
  -i  universal MAF/TAF (cactus-hal2maf --universal); its .tui is built if missing
  -o  output bigWig (keyed on row-0 ancestor coords)
  -b  bin width in bp (default 1000)
  -T  threads for the depth pass (CPU-bound).  Default 8.
  -j  threads for the .tui build (MEMORY-bound: each worker holds one genome's
      runs[], up to tens of GB at vertebrate scale).  Set BELOW -T on small-RAM
      boxes to avoid OOM.  Default: same as -T.
  -d  tmp/spill dir for the .tui build + the intermediate bedGraphs (put it on a
      ROOMY fs for big inputs).  Default: the output directory.
  -k  keep the intermediate bedGraphs (default: delete after wigToBigWig)
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
# Named --bin has no chunk constraint -- any positive bin width works.
(( BIN > 0 )) || { echo "ERROR: -b BIN ($BIN) must be > 0" >&2; exit 1; }

OUTDIR=$(dirname "$OUT"); mkdir -p "$OUTDIR"
WORK="${TMPDIR_OVR:-$OUTDIR}"; mkdir -p "$WORK"
BASE=$(basename "$OUT"); BASE=${BASE%.bw}; BASE=${BASE%.bigWig}
PARTIALS="$WORK/$BASE.partials.bg"   # raw per-(name,bin) records from `depth`
BG="$WORK/$BASE.bg"                  # merged 4-col mean-depth bedGraph
SIZES="$OUTDIR/$BASE.sizes"

# --- 1. universal column index (.tui; needed by `taffy stats -s` below) ------
[ -n "$TUI_THREADS" ] || TUI_THREADS=$THREADS
if [ -e "$INPUT.tui" ]; then
  echo "[1/4] using existing $INPUT.tui" >&2
else
  echo "[1/4] building $INPUT.tui (taffy index -u -T $TUI_THREADS) ..." >&2
  "$TAFFY" index -i "$INPUT" -u -T "$TUI_THREADS" ${TMPDIR_OVR:+-d "$TMPDIR_OVR"}
fi

# --- 2. named binned depth (per-(name,bin) records, UNSORTED) ----------------
echo "[2/4] taffy depth --coords named --bin $BIN --depth -T $THREADS ..." >&2
"$TAFFY" depth -i "$INPUT" --coords named --bin "$BIN" --minLeaves 1 \
  -T "$THREADS" --depth "$PARTIALS"

# --- 3. sort + merge per (name,bin) -> clean mean-depth bedGraph -------------
echo "[3/4] sort + merge -> $BG ..." >&2
LC_ALL=C sort -k1,1 -k2,2n -T "$WORK" "$PARTIALS" \
  | awk -v N="$BIN" -f "$SCRIPTDIR/uni_depth_merge.awk" > "$BG"

# --- 4. chrom sizes (ancestor refChrs + leaves) + bigWig --------------------
echo "[4/4] taffy stats -s -> sizes; wigToBigWig -> $OUT ..." >&2
"$TAFFY" stats -i "$INPUT" -s > "$SIZES"
"$WIGTOBIGWIG" "$BG" "$SIZES" "$OUT"
[ "$KEEP_BG" = 1 ] || rm -f "$PARTIALS" "$BG"

echo "" >&2
echo "done: $OUT ($(du -h "$OUT" 2>/dev/null | cut -f1)); sizes: $SIZES" >&2
echo "  keyed on row-0 ancestor coords (e.g. AncNrefChrM); browse directly." >&2
echo "  (taffy lift --bigwig named query is a pending follow-up.)" >&2
