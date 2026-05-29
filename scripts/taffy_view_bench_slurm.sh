#!/bin/bash
#
# taffy view benchmark driver -- SLURM
#
# Benchmarks four region-extract tools on a common (chrom, start=0, length=N)
# query at a log-decade ladder of N values: 1, 10, 100, 1k, ..., chr10_len.
#
# Tools (each gets one cell per N, run in parallel within a size wave):
#   tui      taffy view -i UNI.maf.gz -r REF_CHROM:0-N           (universal MAF + .tui)
#   tai      taffy view -i HG38.maf.gz -r REF_CHROM:0-N          (hg38-anchored + .tai)
#   hal      hal2maf --refGenome HG --refSequence CHR ... HAL    (HAL baseline)
#   bb       bigBedToBed BB CHR 0 N stdout                       (bigbed floor)
#
# Per cell we record wall seconds + max RSS (KB) + exit + timed-out flag +
# output byte count to bench.tsv.
#
# Run model: one SLURM job; within it, each of the 10 sizes runs all four
# cells concurrently (so timings are comparable -- same FS load on every
# tool at a given N) and waves are sequential.  Sub-tool threading: each
# taffy gets T_TOTAL / 4 bgzf threads; the other two are single-threaded.
#
# All four tools READ FROM THE INPUT PATHS AS GIVEN -- there is no local
# stage for now (intentionally; we want to measure realistic
# network-FS-backed behavior).  See --no-stage-local in
# gerp_shard_slurm.sh for the staging pattern when we want to add it.
#
# Usage:
#   taffy_view_bench_slurm.sh \
#       -u UNI.maf.gz   -m HG38.maf.gz  -H HAL.hal  -b BB.bb \
#       -o OUTDIR  -T 48  [options]

set -euo pipefail

UNI=""
HG38=""
HAL=""
BB=""
OUTDIR=""
T_TOTAL=48
REF_CHROM="hg38.chr10"        # taffy view -r prefix (genome.seq for universal MAF)
HAL_GENOME="hg38"             # hal2maf --refGenome
HAL_SEQ="chr10"               # hal2maf --refSequence
BB_CHROM="chr10"              # bigBedToBed positional chrom name
MAX_SIZE=""                   # cap of the size ladder; default = chr10's known length
TIME_BUDGET=1800              # per-cell wall seconds (timeout sends SIGKILL)
SBATCH_TIME=24
SBATCH_MEM=64
PARTITION=""
ACCOUNT=""
DRY_RUN=0
TAFFY="${TAFFY:-$(command -v taffy || true)}"
HAL2MAF="${HAL2MAF:-$(command -v hal2maf || true)}"
BIGBED2BED="${BIGBED2BED:-$(command -v bigBedToBed || true)}"

usage() {
    cat >&2 <<EOF
taffy_view_bench_slurm.sh -- bench 4 region-extract tools on a log-decade
                             ladder, 1 SLURM job, intra-size parallelism

Required:
  -u FILE       Universal MAF (.uni.maf.gz with .tui sibling)
  -m FILE       hg38-anchored MAF (.maf.gz with .tai sibling)
  -H FILE       HAL file (for hal2maf)
  -b FILE       BigBed (for bigBedToBed)
  -o DIR        Output directory

Optional:
  -T INT        Total CPUs (cpus-per-task).  Each taffy gets T/4 bgzf
                threads (default 48 -> 12 per taffy)
  --refChrom NAME   taffy view -r prefix (default $REF_CHROM)
  --halGenome NAME  hal2maf --refGenome (default $HAL_GENOME)
  --halSeq NAME     hal2maf --refSequence (default $HAL_SEQ)
  --bbChrom NAME    bigBedToBed chrom (default $BB_CHROM)
  --maxSize INT     Cap on the size ladder.  Default = chr10's length
                    pulled from \`taffy stats -s\` on the .tai input.
  --timeBudget SEC  Per-cell wall cap (timeout) (default $TIME_BUDGET)
  --time HRS    sbatch wall (default $SBATCH_TIME)
  --mem GB      sbatch mem (default $SBATCH_MEM)
  --partition X --account X
  --dry-run     Print sbatch; do not submit
  -h            Help

Override binary paths via env: TAFFY, HAL2MAF, BIGBED2BED
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)             UNI="$2"; shift 2;;
        -m)             HG38="$2"; shift 2;;
        -H)             HAL="$2"; shift 2;;
        -b)             BB="$2"; shift 2;;
        -o)             OUTDIR="$2"; shift 2;;
        -T)             T_TOTAL="$2"; shift 2;;
        --refChrom)     REF_CHROM="$2"; shift 2;;
        --halGenome)    HAL_GENOME="$2"; shift 2;;
        --halSeq)       HAL_SEQ="$2"; shift 2;;
        --bbChrom)      BB_CHROM="$2"; shift 2;;
        --maxSize)      MAX_SIZE="$2"; shift 2;;
        --timeBudget)   TIME_BUDGET="$2"; shift 2;;
        --time)         SBATCH_TIME="$2"; shift 2;;
        --mem)          SBATCH_MEM="$2"; shift 2;;
        --partition)    PARTITION="$2"; shift 2;;
        --account)      ACCOUNT="$2"; shift 2;;
        --dry-run)      DRY_RUN=1; shift;;
        -h|--help)      usage 0;;
        *)              echo "unknown arg: $1" >&2; usage 1;;
    esac
done

for v in UNI HG38 HAL BB OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: -$(echo $v | cut -c1) required" >&2; usage 1; }
done
[[ -n "$TAFFY"      ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$HAL2MAF"    ]] || { echo "WARN: hal2maf not found; that tool will be skipped" >&2; }
[[ -n "$BIGBED2BED" ]] || { echo "WARN: bigBedToBed not found; that tool will be skipped" >&2; }
[[ -f "$UNI"        ]] || { echo "ERROR: $UNI not found" >&2; exit 1; }
[[ -f "${UNI}.tui"  ]] || { echo "ERROR: $UNI has no .tui sibling" >&2; exit 1; }
[[ -f "$HG38"       ]] || { echo "ERROR: $HG38 not found" >&2; exit 1; }
[[ -f "${HG38}.tai" ]] || { echo "ERROR: $HG38 has no .tai sibling" >&2; exit 1; }
[[ -f "$HAL"        ]] || { echo "ERROR: $HAL not found" >&2; exit 1; }
[[ -f "$BB"         ]] || { echo "ERROR: $BB not found" >&2; exit 1; }

# Default --maxSize: ask the .tai for $REF_CHROM's length so we cap the
# ladder at the actual chrom end (no point bench-ing 134M+1 byte queries).
if [[ -z "$MAX_SIZE" ]]; then
    # No `awk ... exit`: closing the pipe early SIGPIPEs `taffy stats` which
    # propagates as exit 141 under pipefail.  Scan the whole list instead --
    # it's a few hundred ref chroms, microseconds.
    MAX_SIZE=$("$TAFFY" stats -i "$HG38" -s | awk -v c="$REF_CHROM" '$1==c {print $2}')
    if [[ -z "$MAX_SIZE" || ! "$MAX_SIZE" =~ ^[0-9]+$ ]]; then
        echo "ERROR: could not get $REF_CHROM length from $HG38; pass --maxSize explicitly" >&2
        exit 1
    fi
fi

mkdir -p "$OUTDIR" "$OUTDIR/logs"
echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> hal2maf:       ${HAL2MAF:-(skip)}"
echo ">> bigBedToBed:   ${BIGBED2BED:-(skip)}"
echo ">> uni MAF:       $UNI  (+.tui)"
echo ">> hg38 MAF:      $HG38 (+.tai)"
echo ">> HAL:           $HAL"
echo ">> BigBed:        $BB"
echo ">> ref chrom:     $REF_CHROM (hal $HAL_GENOME / $HAL_SEQ, bb $BB_CHROM)"
echo ">> max size:      $MAX_SIZE bp"
echo ">> cpus/task:     $T_TOTAL (each taffy gets $((T_TOTAL / 4)) bgzf threads)"
echo ">> time budget:   $TIME_BUDGET s per cell"

# --- Build the size ladder: 1, 10, ..., chr10_len.  Capped to MAX_SIZE
SIZES=()
n=1
while (( n <= MAX_SIZE )); do
    SIZES+=("$n")
    n=$(( n * 10 ))
done
# Final entry: the chrom end itself (so we always bench the full query).
if (( ${SIZES[-1]:-0} != MAX_SIZE )); then SIZES+=("$MAX_SIZE"); fi
echo ">> sizes:         ${SIZES[*]}"

# --- Generate the runner script (the thing sbatch executes). -----------
RUNNER="$OUTDIR/bench.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -uo pipefail
# Don't 'set -e': we want per-cell exits captured, not job aborted.

UNI="$UNI"
HG38="$HG38"
HAL="$HAL"
BB="$BB"
OUTDIR="$OUTDIR"
T_TOTAL=$T_TOTAL
TAFFY_T=\$(( T_TOTAL / 4 ))
REF_CHROM="$REF_CHROM"
HAL_GENOME="$HAL_GENOME"
HAL_SEQ="$HAL_SEQ"
BB_CHROM="$BB_CHROM"
TIME_BUDGET=$TIME_BUDGET
TAFFY="$TAFFY"
HAL2MAF="$HAL2MAF"
BIGBED2BED="$BIGBED2BED"
SIZES=( ${SIZES[*]} )

BENCH_TSV="\$OUTDIR/bench.tsv"
LOGDIR="\$OUTDIR/logs"
mkdir -p "\$LOGDIR"

# Write header if file is empty / fresh.
if [[ ! -s "\$BENCH_TSV" ]]; then
    printf "tool\tsize_bp\twall_s\tpeak_rss_kb\texit\ttimed_out\tout_bytes\n" > "\$BENCH_TSV"
fi

# Run one cell.  Args: tool_name N cmd...
# Writes a single TSV row to stdout.
run_cell() {
    local tool="\$1" N="\$2"
    shift 2
    local stem="\${tool}_\${N}"
    local time_file="\$LOGDIR/time_\${stem}.txt"
    local out_file="\$LOGDIR/out_\${stem}"   # discarded after wc
    local err_file="\$LOGDIR/err_\${stem}.log"

    # -q suppresses /usr/bin/time's "Command exited with non-zero status N"
    # line on failure (which would otherwise occupy the time_file's first
    # row and break our read).
    /usr/bin/time -q -f '%e %M' -o "\$time_file" \\
        timeout --signal=KILL "\$TIME_BUDGET" "\$@" \\
        > "\$out_file" 2> "\$err_file"
    local rc=\$?

    local wall rss out_bytes timed_out=0
    if [[ -s "\$time_file" ]]; then
        read -r wall rss < "\$time_file"
        # Belt-and-suspenders: if the first line isn't a "<float> <int>"
        # pair, search the file for the last matching line instead.
        if ! [[ "\$wall" =~ ^[0-9.]+$ && "\$rss" =~ ^[0-9]+$ ]]; then
            read -r wall rss < <(awk '/^[0-9.]+[ \t][0-9]+\$/ {l=\$0} END{print l}' "\$time_file")
            [[ -z "\$wall" ]] && { wall="NA"; rss="NA"; }
        fi
    else
        wall="NA"; rss="NA"
    fi
    if (( rc == 137 || rc == 124 )); then timed_out=1; fi
    out_bytes=\$(stat -c %s "\$out_file" 2>/dev/null || echo 0)
    # Drop the output dump; we have its size + we don't want to flood
    # OUTDIR.  Keep the err log for failure debugging.
    rm -f "\$out_file"

    printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\n" \\
        "\$tool" "\$N" "\$wall" "\$rss" "\$rc" "\$timed_out" "\$out_bytes"
}

# Per-size wave: fire 4 tools in parallel, wait, append rows in cell order.
for N in "\${SIZES[@]}"; do
    echo "=== wave N=\$N ==="
    t_wave=\$SECONDS

    declare -A pids
    declare -A rowfiles
    for tool in tui tai hal bb; do
        rowfiles[\$tool]="\$LOGDIR/row_\${tool}_\${N}.tsv"
    done

    # tui: taffy view universal MAF
    ( run_cell tui "\$N" \\
        "\$TAFFY" view -i "\$UNI" -r "\${REF_CHROM}:0-\${N}" -T "\$TAFFY_T" \\
      ) > "\${rowfiles[tui]}" &
    pids[tui]=\$!

    # tai: taffy view hg38-anchored MAF
    ( run_cell tai "\$N" \\
        "\$TAFFY" view -i "\$HG38" -r "\${REF_CHROM}:0-\${N}" -T "\$TAFFY_T" \\
      ) > "\${rowfiles[tai]}" &
    pids[tai]=\$!

    # hal: hal2maf (if installed)
    if [[ -n "\$HAL2MAF" ]]; then
        ( run_cell hal "\$N" \\
            "\$HAL2MAF" --refGenome "\$HAL_GENOME" --refSequence "\$HAL_SEQ" \\
                       --start 0 --length "\$N" "\$HAL" /dev/stdout \\
          ) > "\${rowfiles[hal]}" &
        pids[hal]=\$!
    fi

    # bb: bigBedToBed (if installed)
    if [[ -n "\$BIGBED2BED" ]]; then
        ( run_cell bb "\$N" \\
            "\$BIGBED2BED" "\$BB" "\$BB_CHROM" 0 "\$N" stdout \\
          ) > "\${rowfiles[bb]}" &
        pids[bb]=\$!
    fi

    # Wait for each cell and append its row in canonical order.
    for tool in tui tai hal bb; do
        if [[ -n "\${pids[\$tool]:-}" ]]; then
            wait "\${pids[\$tool]}" || true
            if [[ -s "\${rowfiles[\$tool]}" ]]; then
                cat "\${rowfiles[\$tool]}" >> "\$BENCH_TSV"
            fi
        fi
    done

    echo "=== wave N=\$N took \$((SECONDS - t_wave)) s ==="
done

echo "bench done.  TSV: \$BENCH_TSV"
EOF
chmod +x "$RUNNER"

# --- Companion plot script for the user to run after. -----------------
PLOT="$OUTDIR/plot.py"
cat > "$PLOT" <<'PY'
#!/usr/bin/env python3
"""Plot wall + peak RSS vs query size for each tool (log-log)."""
import csv, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

bench_dir = os.path.dirname(os.path.abspath(__file__))
rows = []
with open(os.path.join(bench_dir, "bench.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            r["size_bp"]    = int(r["size_bp"])
            r["wall_s"]     = float(r["wall_s"]) if r["wall_s"] != "NA" else None
            r["peak_rss_kb"]= float(r["peak_rss_kb"]) if r["peak_rss_kb"] != "NA" else None
            r["timed_out"]  = int(r["timed_out"])
            rows.append(r)
        except ValueError:
            continue

tools = sorted({r["tool"] for r in rows})
colors = {"tui":"#1f77b4", "tai":"#2ca02c", "hal":"#d62728", "bb":"#ff7f0e"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.13, wspace=0.25)

for tool in tools:
    xs, ys_t, ys_m = [], [], []
    for r in rows:
        if r["tool"] != tool: continue
        if r["timed_out"]: continue
        if r["wall_s"] is None or r["peak_rss_kb"] is None: continue
        xs.append(r["size_bp"])
        ys_t.append(r["wall_s"])
        ys_m.append(r["peak_rss_kb"] / 1024.0)  # KB -> MB
    if xs:
        ax1.plot(xs, ys_t, "o-", label=tool, color=colors.get(tool))
        ax2.plot(xs, ys_m, "o-", label=tool, color=colors.get(tool))

for ax, title, ylab in [(ax1, "wall time", "seconds"),
                         (ax2, "peak RSS", "MB")]:
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("query size (bp)")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

out = os.path.join(bench_dir, "bench.png")
fig.savefig(out, dpi=140)
print(f"wrote {out}")
PY
chmod +x "$PLOT"

# --- Submit the SLURM job. --------------------------------------------
SBATCH_ARGS=(
    --cpus-per-task="$T_TOTAL"
    --mem="${SBATCH_MEM}G"
    --time="${SBATCH_TIME}:00:00"
    --output="$OUTDIR/slurm_%j.log"
    -J taffy_view_bench
)
[[ -n "$PARTITION" ]] && SBATCH_ARGS+=( --partition="$PARTITION" )
[[ -n "$ACCOUNT"   ]] && SBATCH_ARGS+=( --account="$ACCOUNT" )

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo ">> DRY RUN -- would submit:"
    echo "sbatch ${SBATCH_ARGS[*]} --parsable $RUNNER"
else
    echo ">> submitting..."
    JOB=$(sbatch "${SBATCH_ARGS[@]}" --parsable "$RUNNER")
    echo ">> job id: $JOB"
fi

echo ">> done."
echo ">> after the job completes:"
echo "     cat $OUTDIR/bench.tsv"
echo "     python3 $OUTDIR/plot.py   # writes bench.png"
