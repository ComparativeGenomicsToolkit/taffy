#!/bin/bash
#
# taffy view RANDOM-SAMPLE benchmark driver -- SLURM
#
# Companion to taffy_view_bench_slurm.sh.  Where that one runs a single
# fixed-location size ladder, this one draws N random genome-wide hg38 regions
# at each of a few fixed sizes and reports the MEAN + min/max range over the
# samples (blockViz-style), so the numbers aren't hostage to one lucky/unlucky
# locus.
#
# Four tools, each extracting the SAME sampled regions to hg38-anchored leaf MAF
# (apples-to-apples -- see feedback_benchmark_apples_to_apples):
#   taf.tui   taffy view -U query --noAncestors -m   (universal .tui -> leaf MAF)
#   tai       taffy view -m                          (hg38-anchored .tai -> MAF)
#   hal       hal2maf                                (HAL -> hg38 MAF)
#   bb        bigMafToMaf                            (bigMaf -> hg38 MAF)
#
# Sampling happens AT RUNTIME on the cluster (the real .taf.gz/.tai live there;
# locally they are stubs).  A fixed seed makes the draw reproducible, and ALL
# four tools see the SAME regions (sampled once into regions.tsv).  Regions are
# drawn length-weighted over the canonical hg38 chroms (chr1..22,X,Y) so the
# sample is uniform over the genome; only chroms >= the size are eligible.
#
# This dataset names hg38 as GCA_000001405.15: the .tui/.tai coordinate space is
# `GCA_000001405.15.chrN`, hal2maf takes --refGenome GCA_000001405.15
# --refSequence chrN, and the bigMaf uses bare `chrN`.  --refGenome carries the
# label; the bare chrN is derived by stripping the prefix.
#
# Per cell we record wall seconds + max RSS (KB) + exit + timed-out flag +
# output byte count, tagged with (size, chrom, start, sample), to bench.tsv.
#
# Usage:
#   taffy_view_sample_bench_slurm.sh \
#       -u UNI.taf.gz  -m HG38.maf.gz  -H HAL.hal  -b BB.bb \
#       -o OUTDIR  -T 32  [options]

set -euo pipefail

UNI_TAF=""                    # universal TAF (.taf.gz with .tui sibling) -> taf.tui cells
HG38=""                       # hg38-anchored MAF (.maf.gz with .tai sibling) -> tai cells
HAL=""                        # HAL (hal2maf)
BB=""                         # genome-wide hg38 bigMaf (.bb, bigMafToMaf)
OUTDIR=""
T_TOTAL=32
TAFFY_T=""                    # bgzf threads per taffy cell; default 4 (floored at 1)
REF_GENOME="GCA_000001405.15" # the hg38 label in this dataset (.tui/.tai/hal refGenome)
SIZES_CSV="1000,100000,1000000,5000000,10000000"   # query sizes (bp)
N_SAMPLES=10                  # random regions drawn per size
SEED=20260620                 # fixed RNG seed -> reproducible, same regions for all tools
HAL_MAX_SIZE=""               # skip the hal cell for any size > this (bp).  Empty = no cap.
TIME_BUDGET=3600              # per-cell wall seconds (timeout sends SIGKILL)
HAL_TIME_BUDGET=""            # per-cell cap just for hal; defaults to TIME_BUDGET
SBATCH_TIME=24
SBATCH_MEM=64
TMP_GB=""
STAGE_LOCAL=1
PARTITION=""
ACCOUNT=""
DRY_RUN=0
WAIT=1
TAFFY="${TAFFY:-$(command -v taffy || true)}"
HAL2MAF="${HAL2MAF:-$(command -v hal2maf || true)}"
BIGMAF2MAF="${BIGMAF2MAF:-$(command -v bigMafToMaf || true)}"

usage() {
    cat >&2 <<EOF
taffy_view_sample_bench_slurm.sh -- bench 4 hg38-MAF extract tools on N random
                                    genome-wide regions per size, mean+range

Required:
  -u FILE       Universal TAF (.taf.gz with .tui sibling)  -> taf.tui cells
  -m FILE       hg38-anchored MAF (.maf.gz with .tai sibling) -> tai cells
  -H FILE       HAL file (for hal2maf)
  -b FILE       genome-wide hg38 bigMaf .bb (for bigMafToMaf)
  -o DIR        Output directory

Optional:
  -T INT        Total CPUs = concurrency thread-slot budget (default $T_TOTAL)
  --taffyT INT  bgzf threads per taffy cell (default 4, floored at 1)
  --refGenome NAME  hg38 label in the .tui/.tai/hal (default $REF_GENOME).
                    tui/tai region = \`NAME.chrN:start-end\`; hal --refGenome
                    NAME --refSequence chrN; bigMaf uses bare chrN.
  --sizes CSV   Query sizes in bp (default $SIZES_CSV)
  --nSamples INT  Random regions drawn per size (default $N_SAMPLES)
  --seed INT    RNG seed for the draw (default $SEED).  Same seed -> same
                regions; all four tools always see the same regions.
  --halMaxSize INT  Skip the hal2maf cell for any size > this (bp).  Default
                unset = hal runs every size.  (Other tools always run all.)
  --timeBudget SEC  Per-cell wall cap (default $TIME_BUDGET)
  --halTimeBudget SEC  Tighter/looser cap just for hal (default = --timeBudget)
  --time HRS    sbatch wall (default $SBATCH_TIME)
  --mem GB      sbatch mem (default $SBATCH_MEM)
  --tmp GB      Per-task local scratch (sbatch --tmp=N).  Default unset.
  --no-stage-local  Read inputs from network paths (no copy to \$TMPDIR)
  --partition X --account X
  --no-wait     Submit and detach (default: block until SLURM finishes)
  --dry-run     Print sbatch; do not submit
  -h            Help

Override binary paths via env: TAFFY, HAL2MAF, BIGMAF2MAF
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)                 UNI_TAF="$2"; shift 2;;
        -m)                 HG38="$2"; shift 2;;
        -H)                 HAL="$2"; shift 2;;
        -b)                 BB="$2"; shift 2;;
        -o)                 OUTDIR="$2"; shift 2;;
        -T)                 T_TOTAL="$2"; shift 2;;
        --taffyT)           TAFFY_T="$2"; shift 2;;
        --refGenome)        REF_GENOME="$2"; shift 2;;
        --sizes)            SIZES_CSV="$2"; shift 2;;
        --nSamples)         N_SAMPLES="$2"; shift 2;;
        --seed)             SEED="$2"; shift 2;;
        --halMaxSize)       HAL_MAX_SIZE="$2"; shift 2;;
        --timeBudget)       TIME_BUDGET="$2"; shift 2;;
        --halTimeBudget)    HAL_TIME_BUDGET="$2"; shift 2;;
        --time)             SBATCH_TIME="$2"; shift 2;;
        --mem)              SBATCH_MEM="$2"; shift 2;;
        --tmp)              TMP_GB="$2"; shift 2;;
        --no-stage-local)   STAGE_LOCAL=0; shift;;
        --partition)        PARTITION="$2"; shift 2;;
        --account)          ACCOUNT="$2"; shift 2;;
        --no-wait)          WAIT=0; shift;;
        --dry-run)          DRY_RUN=1; shift;;
        -h|--help)          usage 0;;
        *)                  echo "unknown arg: $1" >&2; usage 1;;
    esac
done

for v in UNI_TAF BB OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: missing required input \$$v" >&2; usage 1; }
done
# -m (hg38 MAF -> tai cell) and -H (HAL -> hal2maf cell) are OPTIONAL: omit
# either and that tool is simply not staged or run.  cmp5 of the slide bench
# uses only the universal .tui (taf.tui) and the bigMaf (bb).
[[ -n "$TAFFY"      ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$HAL2MAF"    ]] || { echo "WARN: hal2maf not found; that tool will be skipped" >&2; }
[[ -n "$BIGMAF2MAF" ]] || { echo "WARN: bigMafToMaf not found; that tool will be skipped" >&2; }
[[ -f "$UNI_TAF"        ]] || { echo "ERROR: $UNI_TAF not found" >&2; exit 1; }
[[ -f "${UNI_TAF}.tui"  ]] || { echo "ERROR: $UNI_TAF has no .tui sibling" >&2; exit 1; }
if [[ -n "$HG38" ]]; then
    [[ -f "$HG38"       ]] || { echo "ERROR: $HG38 not found" >&2; exit 1; }
    [[ -f "${HG38}.tai" ]] || { echo "ERROR: $HG38 has no .tai sibling" >&2; exit 1; }
fi
[[ -z "$HAL" || -f "$HAL" ]] || { echo "ERROR: $HAL not found" >&2; exit 1; }
[[ -f "$BB"         ]] || { echo "ERROR: $BB not found" >&2; exit 1; }

[[ "$SIZES_CSV" =~ ^[0-9]+(,[0-9]+)*$ ]] || { echo "ERROR: --sizes must be a CSV of integers (got '$SIZES_CSV')" >&2; exit 1; }
[[ "$N_SAMPLES" =~ ^[0-9]+$ && "$N_SAMPLES" -ge 1 ]] || { echo "ERROR: --nSamples must be a positive integer" >&2; exit 1; }
[[ "$SEED" =~ ^[0-9]+$ ]] || { echo "ERROR: --seed must be a non-negative integer" >&2; exit 1; }
if [[ -n "$HAL_MAX_SIZE" ]]; then
    [[ "$HAL_MAX_SIZE" =~ ^[0-9]+$ ]] || { echo "ERROR: --halMaxSize must be a non-negative integer (got '$HAL_MAX_SIZE')" >&2; exit 1; }
fi

if [[ -z "$TAFFY_T" ]]; then TAFFY_T=4; fi
[[ "$TAFFY_T" =~ ^[0-9]+$ ]] || { echo "ERROR: --taffyT must be a non-negative integer (got '$TAFFY_T')" >&2; exit 1; }
(( TAFFY_T >= 1 )) || TAFFY_T=1

# Bash size array from the CSV (drives the per-size waves; the sampler reads
# the CSV form directly).
IFS=',' read -r -a SIZES <<< "$SIZES_CSV"

mkdir -p "$OUTDIR" "$OUTDIR/logs"
echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> hal2maf:       ${HAL2MAF:-(skip)}"
echo ">> bigMafToMaf:   ${BIGMAF2MAF:-(skip)}"
echo ">> uni TAF:       $UNI_TAF  (+.tui)"
echo ">> hg38 MAF:      ${HG38:-(none -- tai cell skipped)}"
echo ">> HAL:           ${HAL:-(none -- hal2maf cell skipped)}"
echo ">> bigMaf:        $BB"
echo ">> ref genome:    $REF_GENOME (tui/tai NAME.chrN, hal refGenome, bb bare chrN)"
echo ">> sizes:         $SIZES_CSV bp"
echo ">> samples/size:  $N_SAMPLES   seed: $SEED"
echo ">> hal max size:  ${HAL_MAX_SIZE:-(none)} bp  (hal2maf cells skipped above this)"
echo ">> cpus/task:     $T_TOTAL (each taffy cell gets $TAFFY_T bgzf threads; throttle <= $T_TOTAL slots)"
echo ">> time budget:   $TIME_BUDGET s per cell${HAL_TIME_BUDGET:+ (hal: ${HAL_TIME_BUDGET}s)}"
echo ">> local-stage:   $([[ $STAGE_LOCAL -eq 1 ]] && echo "ON (copies to \$TMPDIR)" || echo "OFF (reads from network)")"

if [[ "$STAGE_LOCAL" -eq 1 ]]; then
    STAGE_BYTES=0
    for f in "$UNI_TAF" "${UNI_TAF}.tui" "$HG38" "${HG38}.tai" "$HAL" "$BB"; do
        [[ -f "$f" ]] && STAGE_BYTES=$(( STAGE_BYTES + $(stat -Lc %s "$f" 2>/dev/null || echo 0) ))
    done
    STAGE_GB=$(( STAGE_BYTES / (1024**3) ))
    echo ">> stage-in size: ~${STAGE_GB} GB total"
    if [[ -n "$TMP_GB" ]]; then
        echo ">> --tmp:         ${TMP_GB} GB per task (set explicitly)"
    else
        echo ">> --tmp:         not requested (pass --tmp $(( STAGE_GB + 50 )) only if your cluster enforces it)"
    fi
fi

# --- The region sampler (separate file, quoted heredoc -> no shell escaping).
SAMPLER="$OUTDIR/sample_regions.py"
cat > "$SAMPLER" <<'PYEOF'
#!/usr/bin/env python3
"""Sample N random genome-wide hg38 regions per size, length-weighted over the
canonical chroms, deterministically from a seed.  Reads a `taffy stats -s` dump
(name<ws>length per line), keeps `<refGenome>.chrN` canonical chroms (chr1..22,
X,Y -- no alt/unplaced contigs, which carry underscores and may be absent from
the bigMaf), and prints `size  chrom  start  end` (0-based half-open, bare chrN).
"""
import argparse, random, re

ap = argparse.ArgumentParser()
ap.add_argument("--stats", required=True)
ap.add_argument("--refGenome", required=True)
ap.add_argument("--sizes", required=True)        # CSV
ap.add_argument("--nSamples", type=int, required=True)
ap.add_argument("--seed", type=int, required=True)
a = ap.parse_args()

prefix = a.refGenome + "."
canon = re.compile(r"^chr([0-9]+|[XY])$")
chroms = []
with open(a.stats) as f:
    for line in f:
        p = line.split()
        if len(p) < 2 or not p[1].isdigit():
            continue
        name = p[0]
        if not name.startswith(prefix):
            continue
        bare = name[len(prefix):]
        if canon.match(bare):
            chroms.append((bare, int(p[1])))

if not chroms:
    raise SystemExit("sample_regions: no canonical %s.chrN sequences in %s" % (a.refGenome, a.stats))

sizes = [int(x) for x in a.sizes.split(",") if x]
rng = random.Random(a.seed)
for size in sizes:
    elig = [(c, L) for c, L in chroms if L >= size]
    if not elig:
        continue
    # weight by number of valid start positions -> uniform over the genome
    weights = [L - size + 1 for _, L in elig]
    for _ in range(a.nSamples):
        c, L = rng.choices(elig, weights=weights, k=1)[0]
        start = rng.randint(0, L - size)
        print("%d\t%s\t%d\t%d" % (size, c, start, start + size))
PYEOF
chmod +x "$SAMPLER"

# --- Generate the runner script (the thing sbatch executes). -----------
RUNNER="$OUTDIR/bench.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -uo pipefail
# Don't 'set -e': we want per-cell exits captured, not the job aborted.

UNI_TAF="$UNI_TAF"
HG38="$HG38"
HAL="$HAL"
BB="$BB"
OUTDIR="$OUTDIR"
T_TOTAL=$T_TOTAL
TAFFY_T=$TAFFY_T
REF_GENOME="$REF_GENOME"
SIZES_CSV="$SIZES_CSV"
N_SAMPLES=$N_SAMPLES
SEED=$SEED
TIME_BUDGET=$TIME_BUDGET
HAL_TIME_BUDGET=${HAL_TIME_BUDGET:-$TIME_BUDGET}
HAL_MAX_SIZE="$HAL_MAX_SIZE"
STAGE_LOCAL=$STAGE_LOCAL
TAFFY="$TAFFY"
HAL2MAF="$HAL2MAF"
BIGMAF2MAF="$BIGMAF2MAF"
SIZES=( ${SIZES[*]} )

BENCH_TSV="\$OUTDIR/bench.tsv"
LOGDIR="\$OUTDIR/logs"
mkdir -p "\$LOGDIR"

# --- Stage inputs to local scratch (\$TMPDIR or /tmp fallback). -----
# PARALLEL cp: launch every copy at once and wait for the lot, so the wall
# time is the single largest file (the HAL), not the sum.  Trap-cleanup so an
# aborted job doesn't leave TB of leftover scratch.
if [[ "\$STAGE_LOCAL" -eq 1 ]]; then
    SCRATCH="\${TMPDIR:-/tmp/taffy_view_sample_\${SLURM_JOB_ID:-\$\$}}"
    STAGE_DIR="\$SCRATCH/taffy_view_sample_stage_\${SLURM_JOB_ID:-\$\$}"
    mkdir -p "\$STAGE_DIR"
    trap 'rm -rf "\$STAGE_DIR" 2>/dev/null || true' EXIT
    stage_pids=()
    stage_bg() {
        local src="\$1"
        [[ -n "\$src" && -f "\$src" ]] || return 0
        local dst="\$STAGE_DIR/\$(basename "\$src")"
        echo "stage: \$src -> \$dst (\$(stat -Lc %s "\$src" 2>/dev/null || echo ?) bytes)" >&2
        ( t0=\$SECONDS; cp "\$src" "\$dst"; \\
          echo "       done: \$(basename "\$src") in \$((SECONDS - t0)) s" >&2 ) &
        stage_pids+=( \$! )
    }
    stage_bg "\$UNI_TAF"; stage_bg "\$UNI_TAF.tui"
    stage_bg "\$HG38";    stage_bg "\$HG38.tai"
    stage_bg "\$HAL"
    stage_bg "\$BB"
    stage_rc=0
    for p in "\${stage_pids[@]}"; do wait "\$p" || stage_rc=1; done
    [[ "\$stage_rc" -eq 0 ]] || { echo "ERROR: a stage-in cp returned non-zero" >&2; exit 1; }
    # VERIFY each staged copy matches its source byte count.  A cp to a full FS
    # can return 0 yet truncate (write succeeds, flush doesn't) -- this 2.4TB
    # stage (incl the 965GB HAL + 623GB bb) is exactly where that bit cmp5: the
    # staged hg38 MAF read failed cryptically.  Catch truncation here.
    for f in "\$UNI_TAF" "\$UNI_TAF.tui" "\$HG38" "\$HG38.tai" "\$HAL" "\$BB"; do
        [[ -f "\$f" ]] || continue
        ss=\$(stat -Lc %s "\$f" 2>/dev/null); ds=\$(stat -c %s "\$STAGE_DIR/\$(basename "\$f")" 2>/dev/null || echo -1)
        [[ "\$ss" == "\$ds" ]] || { echo "ERROR: staged \$(basename "\$f") is \$ds bytes, source is \$ss -- truncated (scratch full?)" >&2; exit 1; }
    done
    # Index files (.tai/.tui) must be NEWER than the data they index, or taffy
    # treats them as stale.  cp doesn't preserve mtimes and parallel staging
    # finishes the tiny .tai long before the huge .gz -- so the staged .tai ends
    # up OLDER than its .gz.  Touch the staged indexes last so they're newer.
    for ix in "\$STAGE_DIR"/*.tai "\$STAGE_DIR"/*.tui; do [[ -f "\$ix" ]] && touch "\$ix"; done
    UNI_TAF="\$STAGE_DIR/\$(basename "\$UNI_TAF")"
    [[ -n "\$HG38" ]] && HG38="\$STAGE_DIR/\$(basename "\$HG38")"
    [[ -n "\$HAL"  ]] && HAL="\$STAGE_DIR/\$(basename "\$HAL")"
    BB="\$STAGE_DIR/\$(basename "\$BB")"
    echo "stage: all inputs staged + size-verified to \$STAGE_DIR" >&2
fi

# --- Sample the regions once (same regions for every tool). ------------
# taffy stats -s gives every ref-seq name + length; the python sampler filters to
# canonical \$REF_GENOME.chrN and draws the windows.  Source is the hg38 MAF when
# given, else the universal TAF (it carries the same \$REF_GENOME.chrN leaf
# sequences) -- so the tui-only cmp5 (no -m) still samples.
CHROM_STATS="\$OUTDIR/hg38.stats.txt"
REGIONS="\$OUTDIR/regions.tsv"
STATS_SRC="\${HG38:-\$UNI_TAF}"
"\$TAFFY" stats -s -i "\$STATS_SRC" > "\$CHROM_STATS" 2> "\$OUTDIR/stats.err" || { echo "ERROR: taffy stats -s failed on \$STATS_SRC:" >&2; cat "\$OUTDIR/stats.err" >&2; exit 1; }
python3 "\$OUTDIR/sample_regions.py" --stats "\$CHROM_STATS" --refGenome "\$REF_GENOME" \\
    --sizes "\$SIZES_CSV" --nSamples "\$N_SAMPLES" --seed "\$SEED" > "\$REGIONS" \\
    || { echo "ERROR: region sampling failed" >&2; exit 1; }
echo "sampled \$(wc -l < "\$REGIONS") regions (\$N_SAMPLES per size, seed \$SEED) over sizes \$SIZES_CSV" >&2

# Truncate + write the header fresh (one job -> one bench.tsv).
printf "tool\tsize_bp\tchrom\tstart\tsample\twall_s\tpeak_rss_kb\texit\ttimed_out\tout_bytes\n" > "\$BENCH_TSV"

# Run one cell.  Args: tool size chrom start sample budget cmd...
# Pipes the tool's stdout through wc -c so we measure produced bytes, NOT the
# cost of writing a big MAF to (network) OUTDIR.  /usr/bin/time measures only
# the left side; PIPESTATUS[0] is the tool's real exit.
run_cell() {
    local tool="\$1" N="\$2" chrom="\$3" rstart="\$4" sample="\$5" budget="\$6"
    shift 6
    local stem="\${tool}_\${N}_\${sample}"
    local time_file="\$LOGDIR/time_\${stem}.txt"
    local out_file="\$LOGDIR/out_\${stem}"
    local err_file="\$LOGDIR/err_\${stem}.log"

    /usr/bin/time -q -f '%e %M' -o "\$time_file" \\
        timeout --signal=KILL "\$budget" "\$@" \\
        2> "\$err_file" | wc -c > "\$out_file"
    local rc=\${PIPESTATUS[0]}

    local wall rss out_bytes timed_out=0
    if [[ -s "\$time_file" ]]; then
        read -r wall rss < "\$time_file"
        if ! [[ "\$wall" =~ ^[0-9.]+\$ && "\$rss" =~ ^[0-9]+\$ ]]; then
            read -r wall rss < <(awk '/^[0-9.]+[ \t][0-9]+\$/ {l=\$0} END{print l}' "\$time_file")
            [[ -z "\$wall" ]] && { wall="NA"; rss="NA"; }
        fi
    else
        wall="NA"; rss="NA"
    fi
    if (( rc == 137 || rc == 124 )); then timed_out=1; fi
    out_bytes=\$(tr -dc '0-9' < "\$out_file" 2>/dev/null); [[ -n "\$out_bytes" ]] || out_bytes=0
    rm -f "\$out_file"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\n" \\
        "\$tool" "\$N" "\$chrom" "\$rstart" "\$sample" "\$wall" "\$rss" "\$rc" "\$timed_out" "\$out_bytes"
}

# --- Concurrency throttle: bound the SUM of running thread-slots by T_TOTAL
# so a size-wave (up to N_SAMPLES*4 cells) doesn't oversubscribe the alloc.
# taffy cell = TAFFY_T slots, hal/bb = 1.  Reset each wave (the wave's final
# 'wait' reaps every child).
THREAD_BUDGET=\$T_TOTAL
launched=0
declare -A pid_threads
acquire_slot() {
    local threads=\$1
    (( threads > THREAD_BUDGET )) && threads=\$THREAD_BUDGET
    while (( launched > 0 && launched + threads > THREAD_BUDGET )); do
        wait -n 2>/dev/null || true
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
register_pid() { pid_threads[\$1]=\$2; }

# --- Per-size wave: every sampled region of this size, all 4 tools, throttled.
for N in "\${SIZES[@]}"; do
    echo "=== size wave N=\$N ===" >&2
    t_wave=\$SECONDS
    unset pids rowfiles; declare -A pids rowfiles   # CLEAR stale keys: a declare-A on an existing
                                                    # assoc array does NOT reset it, so a wave where
                                                    # hal is skipped (size > HAL_MAX_SIZE) would else
                                                    # re-append the prior wave's hal rowfiles.
    launched=0
    unset pid_threads; declare -A pid_threads

    sample=0
    while IFS=\$'\t' read -r rsize chrom rstart rend; do
        [[ "\$rsize" == "\$N" ]] || continue
        sample=\$(( sample + 1 ))
        region="\${REF_GENOME}.\${chrom}:\${rstart}-\${rend}"

        # taf.tui: universal -> leaf hg38 MAF
        key="taf.tui_\${sample}"
        rowfiles[\$key]="\$LOGDIR/row_taf.tui_\${N}_\${sample}.tsv"
        acquire_slot \$TAFFY_T
        ( run_cell taf.tui "\$N" "\$chrom" "\$rstart" "\$sample" "\$TIME_BUDGET" \\
            "\$TAFFY" view -i "\$UNI_TAF" -r "\$region" -U query --noAncestors -m -T "\$TAFFY_T" \\
          ) > "\${rowfiles[\$key]}" &
        pids[\$key]=\$!; register_pid \$! \$TAFFY_T

        # tai: hg38-anchored MAF (only when a --maf/.tai was provided)
        if [[ -n "\$HG38" ]]; then
        key="tai_\${sample}"
        rowfiles[\$key]="\$LOGDIR/row_tai_\${N}_\${sample}.tsv"
        acquire_slot \$TAFFY_T
        ( run_cell tai "\$N" "\$chrom" "\$rstart" "\$sample" "\$TIME_BUDGET" \\
            "\$TAFFY" view -i "\$HG38" -r "\$region" -m -T "\$TAFFY_T" \\
          ) > "\${rowfiles[\$key]}" &
        pids[\$key]=\$!; register_pid \$! \$TAFFY_T
        fi

        # hal: hal2maf (respect HAL_MAX_SIZE)
        if [[ -n "\$HAL" ]] && [[ -n "\$HAL2MAF" ]] && { [[ -z "\$HAL_MAX_SIZE" ]] || (( N <= HAL_MAX_SIZE )); }; then
            key="hal_\${sample}"
            rowfiles[\$key]="\$LOGDIR/row_hal_\${N}_\${sample}.tsv"
            acquire_slot 1
            ( run_cell hal "\$N" "\$chrom" "\$rstart" "\$sample" "\$HAL_TIME_BUDGET" \\
                "\$HAL2MAF" --refGenome "\$REF_GENOME" --refSequence "\$chrom" \\
                           --start "\$rstart" --length "\$N" "\$HAL" /dev/stdout \\
              ) > "\${rowfiles[\$key]}" &
            pids[\$key]=\$!; register_pid \$! 1
        fi

        # bb: bigMafToMaf -> hg38 MAF (bare chrN; 0-based -start/-end)
        if [[ -n "\$BIGMAF2MAF" ]]; then
            key="bb_\${sample}"
            rowfiles[\$key]="\$LOGDIR/row_bb_\${N}_\${sample}.tsv"
            acquire_slot 1
            ( run_cell bb "\$N" "\$chrom" "\$rstart" "\$sample" "\$TIME_BUDGET" \\
                "\$BIGMAF2MAF" -chrom="\$chrom" -start="\$rstart" -end="\$rend" "\$BB" stdout \\
              ) > "\${rowfiles[\$key]}" &
            pids[\$key]=\$!; register_pid \$! 1
        fi
    done < "\$REGIONS"

    for tool in "\${!pids[@]}"; do
        wait "\${pids[\$tool]}" || true
        [[ -s "\${rowfiles[\$tool]}" ]] && cat "\${rowfiles[\$tool]}" >> "\$BENCH_TSV"
    done
    echo "=== size wave N=\$N took \$((SECONDS - t_wave)) s ===" >&2
done

echo "bench done.  TSV: \$BENCH_TSV" >&2
EOF
chmod +x "$RUNNER"

# --- Companion plot: mean line + min/max band over samples, per tool. -----
PLOT="$OUTDIR/plot.py"
cat > "$PLOT" <<'PY'
#!/usr/bin/env python3
"""Random-sample view bench: for each tool, the MEAN over the sampled regions
with a shaded min/max band, vs query size (blockViz-style).  Three panels:
wall time, peak RSS, output size.  Every cell emits hg38-anchored leaf MAF, so
all three axes are directly comparable.

Timed-out / errored / zero-output samples are excluded from the aggregate; a
stderr note reports how many were dropped per (tool,size).  If EVERY sample of a
tool at a size dropped (e.g. hal2maf timing out at the top size), that point is
simply absent -- the line stops, like a cap.
"""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

bench_dir = os.path.dirname(os.path.abspath(__file__))
rows = []
with open(os.path.join(bench_dir, "bench.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            r["size_bp"]     = int(r["size_bp"])
            r["wall_s"]      = float(r["wall_s"]) if r["wall_s"] != "NA" else None
            r["peak_rss_kb"] = float(r["peak_rss_kb"]) if r["peak_rss_kb"] != "NA" else None
            r["timed_out"]   = int(r["timed_out"])
            r["exit"]        = int(r["exit"])
            r["out_bytes"]   = int(r["out_bytes"])
            rows.append(r)
        except (ValueError, KeyError):
            continue

label = {
    "taf.tui": "taffy view .tui (universal)",
    "tai":     "taffy view .tai",
    "hal":     "hal2maf",
    "bb":      "bigMaf (bigMafToMaf)",
}
colors = {"taf.tui": "#9467bd", "tai": "#2ca02c", "hal": "#d62728", "bb": "#ff7f0e"}
order = ["taf.tui", "tai", "hal", "bb"]
present = {r["tool"] for r in rows}
tools = [t for t in order if t in present] + sorted(present - set(order))

def valid(r):
    return (not r["timed_out"] and r["exit"] == 0 and r["out_bytes"] > 0
            and r["wall_s"] is not None and r["peak_rss_kb"] is not None)

def agg(tool, field, scale):
    by = {}
    for r in rows:
        if r["tool"] != tool or not valid(r):
            continue
        by.setdefault(r["size_bp"], []).append(r[field] * scale)
    out = []
    for s in sorted(by):
        vs = by[s]
        out.append((s, sum(vs) / len(vs), min(vs), max(vs), len(vs)))
    return out

# Report dropped (timed-out / errored) samples per (tool,size).
drop = {}
for r in rows:
    if not valid(r):
        drop.setdefault((r["tool"], r["size_bp"]), 0)
        drop[(r["tool"], r["size_bp"])] += 1
for (t, s), n in sorted(drop.items()):
    print(f"note: {t} size={s}: {n} sample(s) dropped (timeout/error/empty)", file=sys.stderr)

fig, axes = plt.subplots(1, 3, figsize=(21, 6))
fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.12, wspace=0.24)
panels = [("wall_s", 1.0, "wall time", "seconds"),
          ("peak_rss_kb", 1 / 1024.0, "peak RSS", "MB"),
          ("out_bytes", 1 / 1e6, "output size", "MB")]

for ax, (field, scale, title, ylab) in zip(axes, panels):
    nonhal_hi = 0.0   # cap y to the non-hal2maf tools (let a slow hal2maf run off-chart)
    for tool in tools:
        pts = agg(tool, field, scale)
        if not pts:
            continue
        xs   = [p[0] for p in pts]
        mean = [p[1] for p in pts]
        lo   = [p[2] for p in pts]
        hi   = [p[3] for p in pts]
        if tool != "hal":
            nonhal_hi = max(nonhal_hi, max(hi))
        c = colors.get(tool)
        ax.plot(xs, mean, "o-", color=c, label=label.get(tool, tool))
        ax.fill_between(xs, lo, hi, color=c, alpha=0.15)
    # Linear axes (browser-query view); cap y to the non-hal2maf tools so a slow
    # hal2maf curve runs off the top rather than squashing everything else.
    ax.set_xlabel("query size (bp)")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    if nonhal_hi > 0:
        ax.set_ylim(0, nonhal_hi * 1.15)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1e6:g}M" if x >= 1e6 else (f"{x/1e3:g}k" if x >= 1e3 else f"{x:g}")))
    ax.legend(fontsize=8)

nsamples = max((p[4] for t in tools for p in agg(t, "wall_s", 1.0)), default=0)
fig.suptitle(f"577-way hg38 MAF extraction: random genome-wide regions, "
             f"mean ± range over up to {nsamples} samples/size", fontsize=12, y=0.98)
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
    -J taffy_view_sample_bench
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
