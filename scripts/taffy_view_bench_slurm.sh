#!/bin/bash
#
# taffy view benchmark driver -- SLURM
#
# Benchmarks four region-extract tools on a common (chrom, start=0, length=N)
# query at a log-decade ladder of N values: 1, 10, 100, 1k, ..., chr10_len.
#
# Tools (each gets one cell per N, run in parallel within a size wave):
#   tui_taf       view -U query                                  (universal -> TAF)
#   tui_maf       view -U query -m                               (universal -> MAF)
#   tui_taf_norm  view -U query        | taffy norm              (universal -> TAF + merge)
#   tui_maf_norm  view -U query        | taffy norm -k           (universal -> MAF + merge)
#   tai_taf       view                                           (hg38-anchored -> TAF)
#   tai_maf       view -m                                        (hg38-anchored -> MAF)
#   hal           hal2maf                                        (HAL baseline)
#   bb            bigBedToBed                                    (bigbed floor)
#
# -U query reorients .tui-extracted blocks onto the queried genome so
# the output is hg38-anchored (comparable to the .tai path).  Without
# it the universal lift's default `-U ancestor` keeps blocks in their
# native ancestor anchor and emits ~12x more data.
#
# The _norm variants pipe through `taffy norm`, which merges adjacent
# blocks with compatible row sets.  The raw -U query output is highly
# fragmented (one block per original universal-anchored block that
# overlapped), which defeats TAF's per-block coord-amortization.
# Normalization recovers contiguity and brings TAF size down toward
# the .tai path's ratio.
#
# Per cell we record wall seconds + max RSS (KB) + exit + timed-out flag +
# output byte count to bench.tsv.
#
# Run model: one SLURM job; within it, each of the 10 sizes runs all four
# cells concurrently (so timings are comparable -- same FS load on every
# tool at a given N) and waves are sequential.  Sub-tool threading: each
# taffy cell gets TAFFY_T bgzf threads (fixed default 4); the other two are
# single-threaded.  A per-wave concurrency throttle bounds the SUM of live
# thread-slots by T_TOTAL so a wave of taffy cells doesn't oversubscribe.
#
# Local-stage mode (default ON)
# -----------------------------
# The four inputs (UNI + .tui, HG38 + .tai, HAL, BB) get copied to $TMPDIR
# at the start of the job; all tools then read from the local copies.
# Without staging the four parallel cells per wave fight each other for
# network-FS bandwidth, which inflates the timings and conflates the
# numbers we actually want to measure.  --no-stage-local skips the copy
# and reads from the original paths (useful only for small inputs or
# local-FS deployments).
#
# Usage:
#   taffy_view_bench_slurm.sh \
#       -u UNI.maf.gz   -m HG38.maf.gz  -H HAL.hal  -b BB.bb \
#       -o OUTDIR  -T 48  [options]

set -euo pipefail

UNI=""
UNI_TAF=""
HG38=""
HAL=""
BB=""
OUTDIR=""
T_TOTAL=48
TAFFY_T=""                    # bgzf threads per taffy cell; default 4 (floored at 1).
                              # A wave fires several taffy cells concurrently, so a
                              # fixed small count + the concurrency throttle (<= T_TOTAL
                              # thread-slots) avoids oversubscription.  Override with -T-ish
                              # tuning via --taffyT.
REF_CHROM="hg38.chr10"        # taffy view -r prefix (genome.seq for universal MAF)
HAL_GENOME="hg38"             # hal2maf --refGenome
HAL_SEQ="chr10"               # hal2maf --refSequence
BB_CHROM="chr10"              # bigBedToBed positional chrom name
START=1000000                 # start of the queried region; default 1 Mb to skip the
                              # telomere-/repeat-dominated 0..1Mb prefix of chr10 that
                              # has near-zero universal-MAF coverage
MAX_SIZE=""                   # cap of the size ladder; default = chr10's length - START
HAL_MAX_SIZE=""               # --halMaxSize: skip (DO NOT LAUNCH) the hal2maf cell for
                              # any ladder size > this (bp).  Empty = no cap (hal runs the
                              # whole ladder).  The HAL has no LOD, so a whole-chrom hal2maf
                              # would burn hours; capping skips it rather than run-then-timeout.
MAF_ONLY=0                    # --mafOnly: emit MAF from every tool -- skip the TAF-output
                              # cells (_taf, _taf_norm, tai_taf), keep tui _maf / _maf_norm,
                              # tai_maf, hal2maf, bb.  Apples-to-apples MAF-extraction compare.
TIME_BUDGET=1800              # per-cell wall seconds (timeout sends SIGKILL)
HAL_TIME_BUDGET=""            # per-cell cap just for hal; defaults to TIME_BUDGET
SBATCH_TIME=24
SBATCH_MEM=64
TMP_GB=""                     # passed to sbatch --tmp= when set; default unset
STAGE_LOCAL=1
PARTITION=""
ACCOUNT=""
DRY_RUN=0
WAIT=1                        # block driver until SLURM finishes; --no-wait to detach
TAFFY="${TAFFY:-$(command -v taffy || true)}"
HAL2MAF="${HAL2MAF:-$(command -v hal2maf || true)}"
BIGBED2BED="${BIGBED2BED:-$(command -v bigBedToBed || true)}"

usage() {
    cat >&2 <<EOF
taffy_view_bench_slurm.sh -- bench 4 region-extract tools on a log-decade
                             ladder, 1 SLURM job, intra-size parallelism

Required (at least one of -u / --uniTaf must be set):
  -u FILE       Universal MAF (.uni.maf.gz with .tui sibling).  If set,
                generates four maf.tui_* cells per size (taffy view -U query
                against the MAF-anchored .tui, with each output format).
  --uniTaf FILE Universal TAF (.uni.taf.gz with .tui sibling).  If set,
                generates four taf.tui_* cells per size.  Both flags can
                be combined to bench both .tui formats side by side.
  -m FILE       hg38-anchored MAF (.maf.gz with .tai sibling)
  -H FILE       HAL file (for hal2maf)
  -b FILE       BigBed (for bigBedToBed)
  -o DIR        Output directory

Optional:
  -T INT        Total CPUs (cpus-per-task) = the concurrency thread-slot
                budget (default 48).  A wave fires several taffy cells in
                parallel; the throttle bounds the SUM of live thread-slots
                by this.
  --taffyT INT  bgzf threads per taffy cell (default 4, floored at 1).
                Fixed -- NOT T/4 -- because a wave runs several taffy cells
                concurrently; the throttle keeps the total within -T.
  --refChrom NAME   taffy view -r prefix (default $REF_CHROM)
  --halGenome NAME  hal2maf --refGenome (default $HAL_GENOME)
  --halSeq NAME     hal2maf --refSequence (default $HAL_SEQ)
  --bbChrom NAME    bigBedToBed chrom (default $BB_CHROM)
  --start INT       Start coord of the queried region (default $START).
                    Each ladder cell N becomes \`REF_CHROM:START-(START+N)\`.
                    Default 1 Mb skips the 0..1Mb chr10 prefix which is
                    telomeric / repeat-dominated and has near-zero
                    universal-MAF coverage; queries there are not
                    representative of typical "what does taffy view do".
  --maxSize INT     Cap on the size ladder (= max length N, not end).
                    Default = chrom_length - START pulled from
                    \`taffy stats -s\` on the .tai input.
  --halMaxSize INT  Skip (DO NOT LAUNCH) the hal2maf cell for any ladder
                    size > this (bp).  Default unset = hal runs the whole
                    ladder.  Because the HAL has no LOD, a whole-chrom
                    hal2maf can burn hours; capping skips those cells
                    outright instead of run-then-timeout.  taffy / tai /
                    bb cells always run the full ladder.
  --mafOnly         Emit MAF from every tool: skip the TAF-output cells
                    (_taf, _taf_norm, tai_taf); keep tui _maf / _maf_norm,
                    tai_maf, hal2maf and bb.  Apples-to-apples MAF compare.
  --timeBudget SEC  Per-cell wall cap (timeout) (default $TIME_BUDGET)
  --halTimeBudget SEC  Tighter cap just for the hal2maf cell.  Useful
                    because hal2maf scales much worse than the other
                    tools and you may want to short-circuit it at,
                    say, 60 s while leaving --timeBudget generous for
                    taffy on the big N values.  Default = --timeBudget.
  --time HRS    sbatch wall (default $SBATCH_TIME)
  --mem GB      sbatch mem (default $SBATCH_MEM)
  --tmp GB      Per-task local scratch requirement (sbatch --tmp=N).
                Default unset.  Required on clusters that enforce
                \`--tmp\`; size to ~(sum of input sizes + 10%).
  --no-stage-local
                Disable copying UNI/.tui, HG38/.tai, HAL, BB to
                \$TMPDIR.  All cells then read from the network paths
                and the 6 parallel cells per wave fight each other
                for FS bandwidth.  Only sensible for small inputs.
  --partition X --account X
  --no-wait     Submit and detach (default: driver blocks until SLURM
                completes).  stdout / stderr split: <OUTDIR>/slurm_<jobid>.log
                and .err.log.
  --dry-run     Print sbatch; do not submit
  -h            Help

Override binary paths via env: TAFFY, HAL2MAF, BIGBED2BED
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)             UNI="$2"; shift 2;;
        --uniTaf)       UNI_TAF="$2"; shift 2;;
        -m)             HG38="$2"; shift 2;;
        -H)             HAL="$2"; shift 2;;
        -b)             BB="$2"; shift 2;;
        -o)             OUTDIR="$2"; shift 2;;
        -T)             T_TOTAL="$2"; shift 2;;
        --taffyT)       TAFFY_T="$2"; shift 2;;
        --refChrom)     REF_CHROM="$2"; shift 2;;
        --halGenome)    HAL_GENOME="$2"; shift 2;;
        --halSeq)       HAL_SEQ="$2"; shift 2;;
        --bbChrom)      BB_CHROM="$2"; shift 2;;
        --start)        START="$2"; shift 2;;
        --maxSize)      MAX_SIZE="$2"; shift 2;;
        --halMaxSize)   HAL_MAX_SIZE="$2"; shift 2;;
        --mafOnly)      MAF_ONLY=1; shift;;
        --timeBudget)    TIME_BUDGET="$2"; shift 2;;
        --halTimeBudget) HAL_TIME_BUDGET="$2"; shift 2;;
        --time)             SBATCH_TIME="$2"; shift 2;;
        --mem)              SBATCH_MEM="$2"; shift 2;;
        --tmp)              TMP_GB="$2"; shift 2;;
        --no-stage-local)   STAGE_LOCAL=0; shift;;
        --partition)        PARTITION="$2"; shift 2;;
        --account)          ACCOUNT="$2"; shift 2;;
        --no-wait)          WAIT=0; shift;;
        --dry-run)          DRY_RUN=1; shift;;
        -h|--help)      usage 0;;
        *)              echo "unknown arg: $1" >&2; usage 1;;
    esac
done

for v in HG38 HAL BB OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: -$(echo $v | cut -c1) required" >&2; usage 1; }
done
[[ -n "$UNI" || -n "$UNI_TAF" ]] || {
    echo "ERROR: at least one of -u / --uniTaf must be set" >&2; usage 1;
}
[[ -n "$TAFFY"      ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$HAL2MAF"    ]] || { echo "WARN: hal2maf not found; that tool will be skipped" >&2; }
[[ -n "$BIGBED2BED" ]] || { echo "WARN: bigBedToBed not found; that tool will be skipped" >&2; }
if [[ -n "$UNI" ]]; then
    [[ -f "$UNI"        ]] || { echo "ERROR: $UNI not found" >&2; exit 1; }
    [[ -f "${UNI}.tui"  ]] || { echo "ERROR: $UNI has no .tui sibling" >&2; exit 1; }
fi
if [[ -n "$UNI_TAF" ]]; then
    [[ -f "$UNI_TAF"        ]] || { echo "ERROR: $UNI_TAF not found" >&2; exit 1; }
    [[ -f "${UNI_TAF}.tui"  ]] || { echo "ERROR: $UNI_TAF has no .tui sibling" >&2; exit 1; }
fi
[[ -f "$HG38"       ]] || { echo "ERROR: $HG38 not found" >&2; exit 1; }
[[ -f "${HG38}.tai" ]] || { echo "ERROR: $HG38 has no .tai sibling" >&2; exit 1; }
[[ -f "$HAL"        ]] || { echo "ERROR: $HAL not found" >&2; exit 1; }
[[ -f "$BB"         ]] || { echo "ERROR: $BB not found" >&2; exit 1; }

# Default --maxSize: ask the .tai for $REF_CHROM's length and subtract
# START so the ladder caps at the actual remaining bp (no point bench-ing
# past the chrom end).
if [[ -z "$MAX_SIZE" ]]; then
    # No `awk ... exit`: closing the pipe early SIGPIPEs `taffy stats` which
    # propagates as exit 141 under pipefail.  Scan the whole list instead --
    # it's a few hundred ref chroms, microseconds.
    CHROM_LEN=$("$TAFFY" stats -i "$HG38" -s | awk -v c="$REF_CHROM" '$1==c {print $2}')
    if [[ -z "$CHROM_LEN" || ! "$CHROM_LEN" =~ ^[0-9]+$ ]]; then
        echo "ERROR: could not get $REF_CHROM length from $HG38; pass --maxSize explicitly" >&2
        exit 1
    fi
    if (( START >= CHROM_LEN )); then
        echo "ERROR: --start $START is at/past $REF_CHROM end ($CHROM_LEN)" >&2
        exit 1
    fi
    MAX_SIZE=$(( CHROM_LEN - START ))
fi

if [[ -n "$HAL_MAX_SIZE" ]]; then
    [[ "$HAL_MAX_SIZE" =~ ^[0-9]+$ ]] || { echo "ERROR: --halMaxSize must be a non-negative integer (got '$HAL_MAX_SIZE')" >&2; exit 1; }
fi

# Per-cell taffy bgzf thread count: fixed default 4, floored at 1 (NOT
# T_TOTAL/4 -- a wave runs several taffy cells in parallel, and the
# concurrency throttle bounds the SUM of thread-slots by T_TOTAL).
if [[ -z "$TAFFY_T" ]]; then
    TAFFY_T=4
fi
[[ "$TAFFY_T" =~ ^[0-9]+$ ]] || { echo "ERROR: --taffyT must be a non-negative integer (got '$TAFFY_T')" >&2; exit 1; }
(( TAFFY_T >= 1 )) || TAFFY_T=1

mkdir -p "$OUTDIR" "$OUTDIR/logs"
echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> hal2maf:       ${HAL2MAF:-(skip)}"
echo ">> bigBedToBed:   ${BIGBED2BED:-(skip)}"
[[ -n "$UNI"     ]] && echo ">> uni MAF:       $UNI  (+.tui)"
[[ -n "$UNI_TAF" ]] && echo ">> uni TAF:       $UNI_TAF  (+.tui)"
echo ">> hg38 MAF:      $HG38 (+.tai)"
echo ">> HAL:           $HAL"
echo ">> BigBed:        $BB"
echo ">> ref chrom:     $REF_CHROM (hal $HAL_GENOME / $HAL_SEQ, bb $BB_CHROM)"
echo ">> start:         $START bp"
echo ">> max size:      $MAX_SIZE bp  (query range = [start, start+N))"
echo ">> hal max size:  ${HAL_MAX_SIZE:-(none)} bp  (hal2maf cells skipped above this; taffy/tai/bb run full ladder)"
echo ">> cpus/task:     $T_TOTAL (each taffy cell gets $TAFFY_T bgzf threads; concurrency throttle <= $T_TOTAL thread-slots)"
echo ">> time budget:   $TIME_BUDGET s per cell${HAL_TIME_BUDGET:+ (hal: ${HAL_TIME_BUDGET}s)}"
echo ">> local-stage:   $([[ $STAGE_LOCAL -eq 1 ]] && echo "ON (copies to \$TMPDIR)" || echo "OFF (reads from network)")"
[[ -n "$TMP_GB" ]] && echo ">> --tmp request: ${TMP_GB} GB per task"

# Belt-and-suspenders sizing hint when local-stage is on.
if [[ "$STAGE_LOCAL" -eq 1 ]]; then
    STAGE_BYTES=0
    STAGE_LIST=()
    [[ -n "$UNI"     ]] && STAGE_LIST+=( "$UNI" "${UNI}.tui" )
    [[ -n "$UNI_TAF" ]] && STAGE_LIST+=( "$UNI_TAF" "${UNI_TAF}.tui" )
    STAGE_LIST+=( "$HG38" "${HG38}.tai" "$HAL" "$BB" )
    for f in "${STAGE_LIST[@]}"; do
        if [[ -f "$f" ]]; then
            STAGE_BYTES=$(( STAGE_BYTES + $(stat -c %s "$f" 2>/dev/null || stat -f %z "$f") ))
        fi
    done
    STAGE_GB=$(( STAGE_BYTES / (1024**3) ))
    echo ">> stage-in size: ~${STAGE_GB} GB total"
    # Promote the stage hint to the actual --tmp default: request local
    # scratch sized to what we stage (+50 GB headroom for the cells' own
    # temp + FS slack).  An explicit --tmp still overrides.
    TMP_GB=${TMP_GB:-$(( STAGE_GB + 50 ))}
    echo ">> --tmp default: ${TMP_GB} GB per task (stage ~${STAGE_GB} GB + 50 GB headroom; override with --tmp)"
fi

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
UNI_TAF="$UNI_TAF"
HG38="$HG38"
HAL="$HAL"
BB="$BB"
OUTDIR="$OUTDIR"
T_TOTAL=$T_TOTAL
# Fixed default bgzf-thread count per taffy cell (NOT T_TOTAL/4): a wave
# fires several taffy cells in parallel, so sizing each at T/4 oversubscribed
# the alloc.  The concurrency throttle below bounds the SUM of live thread
# slots by T_TOTAL.  Floor at 1.
TAFFY_T=$TAFFY_T
REF_CHROM="$REF_CHROM"
HAL_GENOME="$HAL_GENOME"
HAL_SEQ="$HAL_SEQ"
BB_CHROM="$BB_CHROM"
START=$START
TIME_BUDGET=$TIME_BUDGET
HAL_TIME_BUDGET=${HAL_TIME_BUDGET:-$TIME_BUDGET}
HAL_MAX_SIZE="$HAL_MAX_SIZE"
MAF_ONLY=$MAF_ONLY
STAGE_LOCAL=$STAGE_LOCAL
TAFFY="$TAFFY"
HAL2MAF="$HAL2MAF"
BIGBED2BED="$BIGBED2BED"
SIZES=( ${SIZES[*]} )

BENCH_TSV="\$OUTDIR/bench.tsv"
LOGDIR="\$OUTDIR/logs"
mkdir -p "\$LOGDIR"

# --- Stage inputs to local scratch (\$TMPDIR or /tmp fallback). -----
# Copy UNI + .tui, HG38 + .tai, HAL, and BB up front so the 6 parallel
# cells per wave don't fight each other for network-FS bandwidth.
# Sequential cp: friendlier to the source FS than parallel pulls.
# Trap-cleanup so an aborted job doesn't leave TB of leftover scratch.
if [[ "\$STAGE_LOCAL" -eq 1 ]]; then
    SCRATCH="\${TMPDIR:-/tmp/taffy_view_bench_\${SLURM_JOB_ID:-\$\$}}"
    # Per-job subdir under SCRATCH so multiple jobs sharing a node-local
    # \$TMPDIR (e.g. /data/tmp without a per-task suffix) don't race on
    # the stage path.
    STAGE_DIR="\$SCRATCH/taffy_view_bench_stage_\${SLURM_JOB_ID:-\$\$}"
    mkdir -p "\$STAGE_DIR"
    trap 'rm -rf "\$STAGE_DIR" 2>/dev/null || true' EXIT
    # stage_one prints the destination path on stdout (for command
    # substitution capture) and status to stderr so the dst doesn't get
    # polluted by the progress messages.
    stage_one() {
        local src="\$1"
        local dst="\$STAGE_DIR/\$(basename "\$src")"
        echo "stage: \$src -> \$dst (\$(stat -c %s "\$src" 2>/dev/null || echo ?) bytes)" >&2
        local t0=\$SECONDS
        cp "\$src" "\$dst"
        echo "       done in \$((SECONDS - t0)) s" >&2
        echo "\$dst"
    }
    # Universal MAF/.tui: each provided source is staged (source file +
    # its .tui sidecar).  taffy view DOES open the source file (unlike
    # taffy lift) to extract bases, so for view-bench we need the source
    # itself staged.
    if [[ -n "\$UNI" ]]; then
        LOCAL_UNI=\$(stage_one "\$UNI");        stage_one "\$UNI.tui"  > /dev/null
        UNI="\$LOCAL_UNI"
    fi
    if [[ -n "\$UNI_TAF" ]]; then
        LOCAL_UNI_TAF=\$(stage_one "\$UNI_TAF");  stage_one "\$UNI_TAF.tui" > /dev/null
        UNI_TAF="\$LOCAL_UNI_TAF"
    fi
    LOCAL_HG38=\$(stage_one "\$HG38");      stage_one "\$HG38.tai" > /dev/null
    LOCAL_HAL=\$(stage_one "\$HAL")
    LOCAL_BB=\$(stage_one "\$BB")
    HG38="\$LOCAL_HG38"
    HAL="\$LOCAL_HAL"
    BB="\$LOCAL_BB"
    echo "stage: all inputs staged to \$STAGE_DIR" >&2
fi

# Write header if file is empty / fresh.
if [[ ! -s "\$BENCH_TSV" ]]; then
    printf "tool\tsize_bp\twall_s\tpeak_rss_kb\texit\ttimed_out\tout_bytes\n" > "\$BENCH_TSV"
fi

# Run one cell.  Args: tool_name N budget_secs cmd...
# Writes a single TSV row to stdout.
run_cell() {
    local tool="\$1" N="\$2" budget="\$3"
    shift 3
    local stem="\${tool}_\${N}"
    local time_file="\$LOGDIR/time_\${stem}.txt"
    local out_file="\$LOGDIR/out_\${stem}"   # discarded after wc
    local err_file="\$LOGDIR/err_\${stem}.log"

    # -q suppresses /usr/bin/time's "Command exited with non-zero status N"
    # line on failure (which would otherwise occupy the time_file's first
    # row and break our read).
    /usr/bin/time -q -f '%e %M' -o "\$time_file" \\
        timeout --signal=KILL "\$budget" "\$@" \\
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

# --- Concurrency throttle ----------------------------------------------
# Bound the SUM of running-cell thread-slots by T_TOTAL so a wave of taffy
# cells doesn't oversubscribe the SLURM alloc.  Before this, a wave fired
# ~5 \`taffy view\` cells at TAFFY_T threads each on a 32-core alloc, so the
# cells fought for cores and the marquee timings were corrupted.
#
# Each cell charges its thread-slot cost: taffy view = TAFFY_T, a _norm
# pipe = TAFFY_T+1 (view | norm), tai = TAFFY_T, hal2maf = 1, bb = 1.
#
# acquire_slot N : block until \`launched + N <= THREAD_BUDGET\`, then
#                  increment \`launched\`.  Uses \`wait -n\` to sleep on any
#                  child completion, then reaps every dead tracked pid.
# register_pid PID N : record a backgrounded cell's pid + thread count so
#                  the next acquire_slot can decrement \`launched\` when it
#                  finishes.
# The view runner runs per-size WAVES (it \`wait\`s for all cells at the end
# of each wave), so \`launched\` + \`pid_threads\` are reset at the START of
# each wave -- no live children carry across a wave boundary.
THREAD_BUDGET=\$T_TOTAL
launched=0
declare -A pid_threads

acquire_slot() {
    local threads=\$1
    # A single cell that asks for more than the whole budget would never
    # fit; clamp so we just wait for the alloc to drain to empty instead
    # of spinning forever.
    (( threads > THREAD_BUDGET )) && threads=\$THREAD_BUDGET
    while (( launched > 0 && launched + threads > THREAD_BUDGET )); do
        wait -n 2>/dev/null || true   # block until any child finishes
        local p
        for p in "\${!pid_threads[@]}"; do
            if ! kill -0 \$p 2>/dev/null; then
                launched=\$(( launched - pid_threads[\$p] ))
                unset pid_threads[\$p]
            fi
        done
    done
    launched=\$(( launched + threads ))
}
register_pid() {
    pid_threads[\$1]=\$2
}

# Per-size wave: fire all tools in parallel (throttled), wait, append rows.
for N in "\${SIZES[@]}"; do
    END=\$((START + N))
    echo "=== wave N=\$N  range=\${REF_CHROM}:\${START}-\${END} ==="
    t_wave=\$SECONDS

    declare -A pids
    declare -A rowfiles
    # Reset the throttle at the start of each wave: the previous wave's
    # \`wait\` reaped every child, so no live pid_threads entries carry over.
    launched=0
    unset pid_threads; declare -A pid_threads
    # Thread-slot charges for acquire_slot / register_pid (see throttle above).
    TAFFY_NORM_T=\$(( TAFFY_T + 1 ))
    # Cells per wave: 4 per provided .tui source (1 or 2 sources) +
    # 2 tai (taf / maf) + hal + bb.  With both .tui sources set this is
    # up to 12 cells; with one source it's 8 (same as before this flag
    # was added).
    #
    # tui-side cells re-orient onto the queried genome (-U query) so
    # blocks come out hg38-anchored, comparable to the .tai path / bigBed.
    # Without -U query the default is \`-U ancestor\` which emits every
    # overlapping universal-column block (~12x larger, conflates bench
    # question).  -m forces MAF, absence => TAF.  The _norm variants
    # pipe through \`taffy norm\` to merge fragmented universal output;
    # \`sh -c\` captures the pipe's RSS via wait4.
    #
    # prefix = "maf.tui" or "taf.tui" (matches the .tui source format).
    launch_tui_cells() {
        local prefix="\$1" input="\$2"
        local region="\${REF_CHROM}:\${START}-\${END}"
        # MAF-output cells (always run): _maf + _maf_norm.
        rowfiles["\${prefix}_maf"]="\$LOGDIR/row_\${prefix}_maf_\${N}.tsv"
        rowfiles["\${prefix}_maf_norm"]="\$LOGDIR/row_\${prefix}_maf_norm_\${N}.tsv"
        acquire_slot \$TAFFY_T
        ( run_cell "\${prefix}_maf"      "\$N" "\$TIME_BUDGET" \\
            "\$TAFFY" view -i "\$input" -r "\$region" -U query -m -T "\$TAFFY_T" \\
          ) > "\${rowfiles[\${prefix}_maf]}" &
        pids[\${prefix}_maf]=\$!; register_pid \$! \$TAFFY_T
        acquire_slot \$TAFFY_NORM_T
        ( run_cell "\${prefix}_maf_norm" "\$N" "\$TIME_BUDGET" \\
            sh -c '"\$1" view -i "\$2" -r "\$3" -U query -T "\$4" | "\$1" norm -k' \\
            _ "\$TAFFY" "\$input" "\$region" "\$TAFFY_T" \\
          ) > "\${rowfiles[\${prefix}_maf_norm]}" &
        pids[\${prefix}_maf_norm]=\$!; register_pid \$! \$TAFFY_NORM_T
        # TAF-output cells: skipped under --mafOnly.
        if [[ "\$MAF_ONLY" -eq 0 ]]; then
            rowfiles["\${prefix}_taf"]="\$LOGDIR/row_\${prefix}_taf_\${N}.tsv"
            rowfiles["\${prefix}_taf_norm"]="\$LOGDIR/row_\${prefix}_taf_norm_\${N}.tsv"
            acquire_slot \$TAFFY_T
            ( run_cell "\${prefix}_taf"      "\$N" "\$TIME_BUDGET" \\
                "\$TAFFY" view -i "\$input" -r "\$region" -U query -T "\$TAFFY_T" \\
              ) > "\${rowfiles[\${prefix}_taf]}" &
            pids[\${prefix}_taf]=\$!; register_pid \$! \$TAFFY_T
            acquire_slot \$TAFFY_NORM_T
            ( run_cell "\${prefix}_taf_norm" "\$N" "\$TIME_BUDGET" \\
                sh -c '"\$1" view -i "\$2" -r "\$3" -U query -T "\$4" | "\$1" norm' \\
                _ "\$TAFFY" "\$input" "\$region" "\$TAFFY_T" \\
              ) > "\${rowfiles[\${prefix}_taf_norm]}" &
            pids[\${prefix}_taf_norm]=\$!; register_pid \$! \$TAFFY_NORM_T
        fi
    }

    # Pre-allocate non-tui cell rowfiles too.
    for tool in tai_taf tai_maf hal bb; do
        rowfiles[\$tool]="\$LOGDIR/row_\${tool}_\${N}.tsv"
    done

    [[ -n "\$UNI"     ]] && launch_tui_cells "maf.tui" "\$UNI"
    [[ -n "\$UNI_TAF" ]] && launch_tui_cells "taf.tui" "\$UNI_TAF"

    # tai_maf: hg38-anchored MAF via .tai.  tai_taf (TAF output) is skipped
    # under --mafOnly.
    if [[ "\$MAF_ONLY" -eq 0 ]]; then
        acquire_slot \$TAFFY_T
        ( run_cell tai_taf "\$N" "\$TIME_BUDGET" \\
            "\$TAFFY" view -i "\$HG38" -r "\${REF_CHROM}:\${START}-\${END}"    -T "\$TAFFY_T" \\
          ) > "\${rowfiles[tai_taf]}" &
        pids[tai_taf]=\$!; register_pid \$! \$TAFFY_T
    fi

    acquire_slot \$TAFFY_T
    ( run_cell tai_maf "\$N" "\$TIME_BUDGET" \\
        "\$TAFFY" view -i "\$HG38" -r "\${REF_CHROM}:\${START}-\${END}" -m -T "\$TAFFY_T" \\
      ) > "\${rowfiles[tai_maf]}" &
    pids[tai_maf]=\$!; register_pid \$! \$TAFFY_T

    # hal: hal2maf (if installed); uses HAL_TIME_BUDGET (defaults to
    # TIME_BUDGET) so the user can short-circuit it independently.
    # --halMaxSize cap: SKIP (do not launch) the hal cell when N exceeds it.
    # The HAL has no LOD, so a whole-chrom hal2maf would burn hours; we skip
    # it outright rather than run-then-timeout.  No hal row is emitted for
    # capped sizes (the plot simply shows the hal line stopping at the cap).
    if [[ -n "\$HAL2MAF" ]] && { [[ -z "\$HAL_MAX_SIZE" ]] || (( N <= HAL_MAX_SIZE )); }; then
        acquire_slot 1
        ( run_cell hal "\$N" "\$HAL_TIME_BUDGET" \\
            "\$HAL2MAF" --refGenome "\$HAL_GENOME" --refSequence "\$HAL_SEQ" \\
                       --start "\$START" --length "\$N" "\$HAL" /dev/stdout \\
          ) > "\${rowfiles[hal]}" &
        pids[hal]=\$!; register_pid \$! 1
    elif [[ -n "\$HAL2MAF" ]]; then
        echo "   (hal skipped at N=\$N > halMaxSize=\$HAL_MAX_SIZE)" >&2
    fi

    # bb: bigBedToBed (if installed).  Uses the flag-based CLI
    # (-chrom=/-start=/-end=) which is what every UCSC build since ~2015
    # has shipped; the positional "chrom start end" form some older
    # builds had isn't there in v1 (kent/src/utils/bigBedToBed.c).
    if [[ -n "\$BIGBED2BED" ]]; then
        acquire_slot 1
        ( run_cell bb "\$N" "\$TIME_BUDGET" \\
            "\$BIGBED2BED" -chrom="\$BB_CHROM" -start="\$START" -end="\$END" "\$BB" stdout \\
          ) > "\${rowfiles[bb]}" &
        pids[bb]=\$!; register_pid \$! 1
    fi

    # Wait for each cell and append its row.  Iterate the pids set we
    # actually launched -- which varies by which .tui sources were given.
    for tool in "\${!pids[@]}"; do
        wait "\${pids[\$tool]}" || true
        if [[ -s "\${rowfiles[\$tool]}" ]]; then
            cat "\${rowfiles[\$tool]}" >> "\$BENCH_TSV"
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
"""Plot wall + peak RSS vs query size for each tool (log-log).

Timed-out cells get a hollow 'X' at the budget value -- communicates
"killed at budget" instead of silently dropping the point.  Cells with
non-zero exit and zero output (bigBedToBed range failures at large N)
are dropped: they are not real measurements.
"""
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
            r["exit"]       = int(r["exit"])
            r["out_bytes"]  = int(r["out_bytes"])
            rows.append(r)
        except (ValueError, KeyError):
            continue

# maf.tui_* uses blues; taf.tui_* uses purples (visually distinct but
# easy to read pairs across colour-coded output formats).  tai_*, hal,
# bb keep their existing colours.
colors = {
    "maf.tui_taf":      "#1f77b4", "maf.tui_maf":      "#aec7e8",
    "maf.tui_taf_norm": "#17becf", "maf.tui_maf_norm": "#9edae5",
    "taf.tui_taf":      "#9467bd", "taf.tui_maf":      "#c5b0d5",
    "taf.tui_taf_norm": "#bd5fbb", "taf.tui_maf_norm": "#e6abe6",
    "tai_taf":          "#2ca02c", "tai_maf":          "#98df8a",
    "hal":              "#d62728", "bb":               "#ff7f0e",
}
# Show the tools present in the data (the driver may have launched only
# the maf.tui set, or only the taf.tui set, or both).  Preserve a stable
# render order so legend reads cleanly.
order = ["maf.tui_taf", "maf.tui_maf", "maf.tui_taf_norm", "maf.tui_maf_norm",
         "taf.tui_taf", "taf.tui_maf", "taf.tui_taf_norm", "taf.tui_maf_norm",
         "tai_taf", "tai_maf", "hal", "bb"]
present = {r["tool"] for r in rows}
tools = [t for t in order if t in present]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.subplots_adjust(left=0.06, right=0.97, top=0.92, bottom=0.12, wspace=0.22)

for tool in tools:
    rs = [r for r in rows if r["tool"] == tool
          and r["wall_s"] is not None and r["peak_rss_kb"] is not None]
    rs.sort(key=lambda r: r["size_bp"])
    ok_x, ok_w, ok_r = [], [], []
    to_x, to_w        = [], []
    for r in rs:
        x = r["size_bp"]; w = r["wall_s"]; m = r["peak_rss_kb"] / 1024.0
        if r["timed_out"]:
            to_x.append(x); to_w.append(w)
        elif r["exit"] != 0 or (r["out_bytes"] == 0 and x > 0):
            # bigBedToBed range fail / other no-output -- not a measurement.
            continue
        else:
            ok_x.append(x); ok_w.append(w); ok_r.append(m)
    color = colors.get(tool)
    if ok_x:
        ax1.plot(ok_x, ok_w, "o-", label=tool, color=color)
        ax2.plot(ok_x, ok_r, "o-", label=tool, color=color)
    if to_x:
        ax1.plot(to_x, to_w, "X", color=color, markerfacecolor="none",
                 markeredgewidth=2, markersize=11, label=f"{tool} (timed out)")

for ax, title, ylab in [(ax1, "wall time", "seconds"),
                         (ax2, "peak RSS", "MB")]:
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("query size (bp)")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

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
    --error="$OUTDIR/slurm_%j.err.log"
    -J taffy_view_bench
)
[[ -n "$TMP_GB"    ]] && SBATCH_ARGS+=( --tmp="${TMP_GB}G" )
[[ -n "$PARTITION" ]] && SBATCH_ARGS+=( --partition="$PARTITION" )
[[ -n "$ACCOUNT"   ]] && SBATCH_ARGS+=( --account="$ACCOUNT" )

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo ">> DRY RUN -- would submit:"
    echo "sbatch ${SBATCH_ARGS[*]} --parsable $RUNNER"
else
    echo ">> submitting..."
    JOB=$(sbatch "${SBATCH_ARGS[@]}" --parsable "$RUNNER")
    echo ">> job id: $JOB"
    echo ">> stdout: $OUTDIR/slurm_${JOB}.log"
    echo ">> stderr: $OUTDIR/slurm_${JOB}.err.log"
fi

# Block driver until the job finishes (poll squeue every 60s, final
# state via sacct).  --no-wait skips this.
if [[ "$DRY_RUN" -ne 1 && "$WAIT" -eq 1 ]]; then
    echo ">> waiting for job $JOB ..."
    while squeue -j "$JOB" -h -o "%T" 2>/dev/null | grep -qE "PENDING|RUNNING|CONFIGURING|COMPLETING|RESIZING|SUSPENDED|REQUEUED"; do
        sleep 60
    done
    FINAL_STATE=$(sacct -j "$JOB" -X -n -o State 2>/dev/null | head -1 | tr -d ' ')
    echo ">> job $JOB final state: ${FINAL_STATE:-UNKNOWN}"
    case "$FINAL_STATE" in
        COMPLETED) ;;
        *)         echo ">> NON-SUCCESS state -- check $OUTDIR/slurm_${JOB}.err.log"; exit 1;;
    esac
fi

echo ">> done."
echo ">> after the job completes:"
echo "     cat $OUTDIR/bench.tsv"
echo "     python3 $OUTDIR/plot.py   # writes bench.png"
