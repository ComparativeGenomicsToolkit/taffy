#!/usr/bin/env bash
# Local (non-SLURM) version of taffy_view_bench_slurm.sh, scoped to the
# specific question: how does `taffy view -r` on a chained .tui sidecar
# compare against the base .tui and against UCSC bigBedToBed?
#
# Strips out: SLURM, parallel-within-wave, --no-stage-local, multi-input
# tui formats.  Keeps: log-decade size ladder, /usr/bin/time -f '%e %M'
# timing, per-cell timeout, output-byte capture, TSV emit.
#
# Tools per size (each runs sequentially -- no within-wave parallelism
# so the timings aren't I/O-contended):
#   base.tui_maf      taffy view -U query -m   on the base .uni.maf.gz
#   chained.tui_maf   taffy view -U query -m   on the chained .uni.maf.gz
#                     (via <input>.chained.maf.gz symlink + .tui/.tai siblings)
#   bb                bigBedToBed              on the hg38-anchored bigMAF
#
# Notes (DO NOT remove):
#   * The chained .tui sidecar needs three siblings next to whatever
#     path `taffy view -i` is given: <X>.maf.gz, <X>.maf.gz.tui,
#     <X>.maf.gz.tai.  We use symlinks; see the apes setup in
#     /home/hickey/dev/work/unitaf for the pattern.
#   * taffy view -U query reorients output onto the queried genome so it's
#     hg38-anchored, comparable to bigBedToBed which is also hg38-anchored.
#     Without -U query the default `-U ancestor` emits ancestor-anchored
#     blocks (often ~12x more data) and the comparison is meaningless.

set -uo pipefail

BASE=""
CHAINED_BASE=""
REF_MAF=""           # optional: reference-anchored MAF (with .tai sidecar)
BB=""
OUTDIR=""
REF_CHROM="hg38.chr20"
BB_CHROM="chr20"
START=1000000
MAX_SIZE=""
TIME_BUDGET=600
TAFFY="${TAFFY:-$(command -v taffy || true)}"
BIGBED2BED="${BIGBED2BED:-$(command -v bigBedToBed || echo /home/hickey/local/bin/bigBedToBed)}"

usage() {
    cat >&2 <<EOF
taffy_view_bench_local.sh -- local view bench, 3 tools across a size ladder

Required:
  -u FILE       Base .uni.maf.gz (with .tui sibling).  No .tai needed --
                .uni.maf.gz files are queried via universal lift (.tui).
  -c FILE       Chained-sidecar base path.  Symlinks next to <c>.maf.gz
                and <c>.maf.gz.tui must already exist; the script does
                NOT create them.  E.g. /path/foo.chained where
                /path/foo.chained.maf.gz -> source uni.maf.gz, .tui ->
                chained .tui.
  -b FILE       Reference-anchored bigMAF (.bb)
  -o DIR        Output directory (created if missing)

Optional fourth tool:
  --refMaf FILE Reference-anchored MAF (with .tai sibling).  This is the
                "normal" tai-indexed MAF (e.g. mm39.maf.gz or
                mexican_tetra.maf.gz) -- NOT the universal MAF.  When set,
                a 4th tool tai_maf is added per size, running
                  taffy view -i <FILE> -r <bbChrom>:<start>-<end> -m

Optional:
  --refChrom NAME      view -r prefix (default $REF_CHROM)
  --bbChrom NAME       bigBedToBed chrom (default $BB_CHROM)
  --start INT          start of the queried region (default $START)
  --maxSize INT        cap of the size ladder (default = chrom_length - start
                       from taffy stats -s on the base input)
  --timeBudget SEC     per-cell timeout (default $TIME_BUDGET)
  -h --help            this help
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)            BASE="$2"; shift 2;;
        -c)            CHAINED_BASE="$2"; shift 2;;
        -b)            BB="$2"; shift 2;;
        --refMaf)      REF_MAF="$2"; shift 2;;
        -o)            OUTDIR="$2"; shift 2;;
        --refChrom)    REF_CHROM="$2"; shift 2;;
        --bbChrom)     BB_CHROM="$2"; shift 2;;
        --start)       START="$2"; shift 2;;
        --maxSize)     MAX_SIZE="$2"; shift 2;;
        --timeBudget)  TIME_BUDGET="$2"; shift 2;;
        -h|--help)     usage 0;;
        *)             echo "unknown arg: $1" >&2; usage 1;;
    esac
done

for v in BASE CHAINED_BASE BB OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: missing required" >&2; usage 1; }
done
[[ -n "$TAFFY" ]] || { echo "ERROR: taffy not on PATH" >&2; exit 1; }
[[ -x "$BIGBED2BED" ]] || { echo "ERROR: bigBedToBed not at $BIGBED2BED" >&2; exit 1; }
[[ -f "$BASE" && -f "$BASE.tui" ]] || \
    { echo "ERROR: base needs $BASE + .tui" >&2; exit 1; }
[[ -f "$CHAINED_BASE.maf.gz" && -f "$CHAINED_BASE.maf.gz.tui" ]] || \
    { echo "ERROR: chained base needs $CHAINED_BASE.maf.gz + .tui siblings" >&2; exit 1; }
# .tai is NOT required for view -r on a universal MAF: with .tui present
# view uses the universal-lift path (any genome.seq query).  .tai is only
# needed for the legacy row-0-anchored path (which doesn't apply here).
[[ -f "$BB" ]] || { echo "ERROR: $BB not found" >&2; exit 1; }
if [[ -n "$REF_MAF" ]]; then
    [[ -f "$REF_MAF" && -f "$REF_MAF.tai" ]] || \
        { echo "ERROR: --refMaf needs $REF_MAF + .tai" >&2; exit 1; }
fi

mkdir -p "$OUTDIR" "$OUTDIR/logs"

# Resolve max size if not given.
if [[ -z "$MAX_SIZE" ]]; then
    CHROM_LEN=$("$TAFFY" stats -i "$BASE" -s 2>/dev/null | awk -v c="$REF_CHROM" '$1==c {print $2}')
    if [[ -z "$CHROM_LEN" || ! "$CHROM_LEN" =~ ^[0-9]+$ ]]; then
        echo "ERROR: could not get $REF_CHROM length; pass --maxSize" >&2; exit 1
    fi
    MAX_SIZE=$(( CHROM_LEN - START ))
fi

# Log-decade ladder + chrom end.
SIZES=()
n=1
while (( n <= MAX_SIZE )); do
    SIZES+=("$n")
    n=$(( n * 10 ))
done
(( ${SIZES[-1]:-0} != MAX_SIZE )) && SIZES+=("$MAX_SIZE")

BENCH="$OUTDIR/bench.tsv"
printf "tool\tsize_bp\twall_s\tpeak_rss_kb\texit\ttimed_out\tout_bytes\n" > "$BENCH"

echo ">> base:         $BASE"
echo ">> chained base: $CHAINED_BASE.maf.gz"
echo ">> bb:           $BB"
echo ">> ref chrom:    $REF_CHROM   (bb chrom: $BB_CHROM)"
echo ">> start:        $START"
echo ">> sizes:        ${SIZES[*]}"
echo ">> outdir:       $OUTDIR"
echo

# One cell.  Args: tool N cmd...
run_cell() {
    local tool="$1" N="$2"; shift 2
    local stem="${tool}_${N}"
    local time_file="$OUTDIR/logs/time_${stem}.txt"
    local out_file="$OUTDIR/logs/out_${stem}"
    local err_file="$OUTDIR/logs/err_${stem}.log"
    /usr/bin/time -q -f '%e %M' -o "$time_file" \
        timeout --signal=KILL "$TIME_BUDGET" "$@" \
        > "$out_file" 2> "$err_file"
    local rc=$?
    local wall rss
    if [[ -s "$time_file" ]]; then
        read -r wall rss < "$time_file"
        [[ "$wall" =~ ^[0-9.]+$ && "$rss" =~ ^[0-9]+$ ]] || { wall="NA"; rss="NA"; }
    else
        wall="NA"; rss="NA"
    fi
    local timed_out=0
    (( rc == 137 || rc == 124 )) && timed_out=1
    local out_bytes
    out_bytes=$(stat -c %s "$out_file" 2>/dev/null || echo 0)
    rm -f "$out_file"
    printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\n" \
        "$tool" "$N" "$wall" "$rss" "$rc" "$timed_out" "$out_bytes" >> "$BENCH"
    printf "  %-18s N=%-12s wall=%-7s rss=%-10s bytes=%-12s exit=%d\n" \
        "$tool" "$N" "$wall" "$rss" "$out_bytes" "$rc"
}

for N in "${SIZES[@]}"; do
    END=$((START + N))
    echo "=== N=$N  ${REF_CHROM}:${START}-${END} ==="
    run_cell base.tui_maf "$N" \
        "$TAFFY" view -i "$BASE" -r "${REF_CHROM}:${START}-${END}" -U query -m
    run_cell chained.tui_maf "$N" \
        "$TAFFY" view -i "$CHAINED_BASE.maf.gz" -r "${REF_CHROM}:${START}-${END}" -U query -m
    if [[ -n "$REF_MAF" ]]; then
        # tai_maf: query the reference-anchored MAF by its row-0 seq name.
        # The .tai uses the SAME genome.seq format the .tui uses (e.g.
        # GCF_023375975.1.NC_064427.1, NOT just NC_064427.1).  Only
        # bigBedToBed uses the bare chrom name.
        run_cell tai_maf "$N" \
            "$TAFFY" view -i "$REF_MAF" -r "${REF_CHROM}:${START}-${END}" -m
    fi
    run_cell bb "$N" \
        "$BIGBED2BED" -chrom="$BB_CHROM" -start="$START" -end="$END" "$BB" stdout
done

echo
echo "done.  bench.tsv: $BENCH"
column -t -s $'\t' "$BENCH"
