#!/bin/bash
#
# taffy index -u (universal .tui builder) -- SLURM wrapper with local-stage
#
# Single-task SLURM job: copies the input MAF/TAF to $TMPDIR, runs
# `taffy index -u` writing the .tui (+ spill files) to scratch, then copies
# the .tui back to the requested output path.  Trap-cleanup removes scratch
# on any exit (success or failure).
#
# Why local-stage?  Phase 1 of the indexer is a streaming scan -- on a
# vertebrate-scale input that's TB of sequential reads from the network FS.
# Phase 2 spills per-genome temp files (also TB-class total).  Doing both
# from scratch instead of the source FS keeps the network I/O bounded to a
# single sequential cp at start and one short cp at end, instead of an
# unbounded random pattern across the whole run.
#
# Usage:
#   taffy_index_slurm.sh -i INPUT.uni.maf.gz [-o OUTPUT.tui] [options]
#
# The default output is sibling of the input: INPUT.tui.

set -euo pipefail

INPUT=""
OUTPUT=""
T_THREADS=1               # taffy index -T  (phase-2 OMP workers; bgzf
                          # auto-capped at TAFFY_MAX_BGZF_THREADS=8 inside)
SBATCH_TIME=48            # hours
SBATCH_MEM=256            # GB
TMP_GB=""                 # sbatch --tmp=N when set; default unset
GENOME_NAMES=""           # optional --genomeNames passthrough
PARTITION=""
ACCOUNT=""
DRY_RUN=0
WAIT=1                    # block driver until SLURM job finishes; --no-wait to detach
TAFFY="${TAFFY:-$(command -v taffy || true)}"

usage() {
    cat >&2 <<EOF
taffy_index_slurm.sh -- single-task SLURM job that builds a .tui in scratch
                       and copies the result to its final location.

Required:
  -i FILE             Input universal MAF or TAF
                      (must be a cactus-hal2maf --universal output)

Optional:
  -o FILE             Output .tui path (default: <INPUT>.tui)
  -T N                Worker threads (default $T_THREADS).  Drives BOTH
                      phase-2 OMP workers AND bgzf decompression.  bgzf
                      side is auto-capped at TAFFY_MAX_BGZF_THREADS=8
                      inside taffy itself (past 8 the htslib mt-pool
                      overhead dominates the decompress win).  Phase-2
                      side uses the full N: each in-flight worker holds
                      ONE genome's runs[] in RAM (~tens of GB for
                      vertebrate-scale giants), so memory peak roughly
                      N * 40 GB.  On a 1.5 TB host with the 577-way:
                      N=8 safe, N=16 borderline, N>=32 may OOM.
  --tmpDir DIR        Override the in-job scratch base.  Default \$TMPDIR.
                      Forwarded to \`taffy index --tmpDir\` so phase-1 spill
                      files also land in scratch (NOT next to the output).
  --genomeNames FILE  Passed through to \`taffy index --genomeNames\`.
                      Use when seq names contain >1 '.' and you need to
                      disambiguate genome vs sequence.
  --time HRS          sbatch --time (default $SBATCH_TIME)
  --mem GB            sbatch --mem  (default $SBATCH_MEM)
  --tmp GB            sbatch --tmp= local scratch requirement (optional).
                      Required on clusters that enforce \`--tmp\`; size to
                      ~(2 * input size) so input + spill files fit.
  --partition X --account X
  --no-wait           Submit and detach (default: driver blocks until
                      SLURM completes).  When blocking, the driver
                      polls squeue every 60 s and exits with non-zero
                      if the job's final state isn't COMPLETED.
                      stdout/stderr split: <OUTPUT>.slurm_<jobid>.log
                      and .err.log.
  --dry-run           Print the sbatch command, don't submit
  -h                  Help

Override taffy binary via env: TAFFY=/path/to/taffy
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)              INPUT="$2"; shift 2;;
        -o)              OUTPUT="$2"; shift 2;;
        -T)              T_THREADS="$2"; shift 2;;
        --tmpDir)        TMPDIR_OVERRIDE="$2"; shift 2;;
        --genomeNames)   GENOME_NAMES="$2"; shift 2;;
        --time)          SBATCH_TIME="$2"; shift 2;;
        --mem)           SBATCH_MEM="$2"; shift 2;;
        --tmp)           TMP_GB="$2"; shift 2;;
        --partition)     PARTITION="$2"; shift 2;;
        --account)       ACCOUNT="$2"; shift 2;;
        --no-wait)       WAIT=0; shift;;
        --dry-run)       DRY_RUN=1; shift;;
        -h|--help)       usage 0;;
        *)               echo "unknown arg: $1" >&2; usage 1;;
    esac
done

[[ -n "$INPUT" ]] || { echo "ERROR: -i required" >&2; usage 1; }
[[ -n "$TAFFY" ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -f "$INPUT" ]] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
[[ -z "$GENOME_NAMES" || -f "$GENOME_NAMES" ]] || {
    echo "ERROR: --genomeNames file not found: $GENOME_NAMES" >&2; exit 1; }
# Default output is sibling of input: INPUT.tui.
if [[ -z "$OUTPUT" ]]; then OUTPUT="${INPUT}.tui"; fi
# Resolve to absolute paths so the runner script doesn't depend on cwd.
INPUT=$(readlink -f "$INPUT")
OUTPUT=$(readlink -f "$OUTPUT" || echo "$OUTPUT")
mkdir -p "$(dirname "$OUTPUT")"

INPUT_BYTES=$(stat -c %s "$INPUT" 2>/dev/null || stat -f %z "$INPUT")
INPUT_GB=$(( INPUT_BYTES / (1024**3) ))

echo ">> input:           $INPUT (${INPUT_GB} GB)"
echo ">> output:          $OUTPUT"
echo ">> taffy:           $TAFFY"
echo ">> threads:         $T_THREADS (phase-2 OMP; bgzf auto-capped at 8 inside taffy)"
[[ -n "${TMPDIR_OVERRIDE:-}" ]] && echo ">> tmpdir override: $TMPDIR_OVERRIDE"
[[ -n "$GENOME_NAMES" ]] && echo ">> genome names:    $GENOME_NAMES"
[[ -n "$TMP_GB" ]] && echo ">> --tmp request:   ${TMP_GB} GB"
if [[ -z "$TMP_GB" ]]; then
    echo ">> hint:            if your cluster advertises --tmp, consider --tmp $(( INPUT_GB * 2 + 50 )) (input + spill budget); otherwise omit"
fi

mkdir -p "$(dirname "$OUTPUT")"
LOG_DIR=$(dirname "$OUTPUT")
RUNNER=$(mktemp -p "$LOG_DIR" .tui_runner_XXXXXX.sh)

cat > "$RUNNER" <<EOF
#!/bin/bash
set -euo pipefail
INPUT="$INPUT"
OUTPUT="$OUTPUT"
TAFFY="$TAFFY"
T_THREADS="$T_THREADS"
GENOME_NAMES="$GENOME_NAMES"

# Resolve scratch.  SLURM sets \$TMPDIR per task; if not set (running
# outside SLURM), use the override or fall back to /tmp.
# Per-job subdir under SCRATCH so concurrent index jobs on the same
# node (or another taffy script's stage) don't race on the stage path.
SCRATCH="${TMPDIR_OVERRIDE:-\${TMPDIR:-/tmp/taffy_index_\${SLURM_JOB_ID:-\$\$}}}"
STAGE_DIR="\$SCRATCH/taffy_index_stage_\${SLURM_JOB_ID:-\$\$}"
mkdir -p "\$STAGE_DIR"
# Trap-cleanup on exit (success OR failure).  Phase-1 spill files alone
# can be TB-class; leaving them on shared scratch is rude.
trap 'rm -rf "\$STAGE_DIR" 2>/dev/null || true' EXIT

BASENAME=\$(basename "\$INPUT")
LOCAL_INPUT="\$STAGE_DIR/\$BASENAME"
LOCAL_TUI="\$STAGE_DIR/\$BASENAME.tui"

echo "[\$(date +%H:%M:%S)] stage-in: \$INPUT -> \$LOCAL_INPUT"
t0=\$SECONDS
cp "\$INPUT" "\$LOCAL_INPUT"
echo "[\$(date +%H:%M:%S)] stage-in done in \$((SECONDS - t0)) s"

GENOME_NAMES_FLAG=""
if [[ -n "\$GENOME_NAMES" ]]; then
    # Stage the genome-names file too (cheap, keeps the run self-contained).
    LOCAL_GN="\$STAGE_DIR/\$(basename "\$GENOME_NAMES")"
    cp "\$GENOME_NAMES" "\$LOCAL_GN"
    GENOME_NAMES_FLAG="-n \$LOCAL_GN"
fi

echo "[\$(date +%H:%M:%S)] indexing"
t0=\$SECONDS
# --tmpDir => phase-1 spills also land in scratch, not next to the output.
# -l INFO so the slurm log captures phase-1 progress ticks (every 600 s)
# + phase-2 per-genome timings -- without this the log is silent through
# multi-hour runs.
"\$TAFFY" index -i "\$LOCAL_INPUT" -u \\
    -T "\$T_THREADS" \\
    --tmpDir "\$STAGE_DIR" \\
    -l INFO \\
    \$GENOME_NAMES_FLAG
echo "[\$(date +%H:%M:%S)] index done in \$((SECONDS - t0)) s"

echo "[\$(date +%H:%M:%S)] stage-out: \$LOCAL_TUI -> \$OUTPUT"
t0=\$SECONDS
# .tmp + rename so a half-copied file never lands at the final path (some
# downstream tools test for existence rather than completeness).
cp "\$LOCAL_TUI" "\$OUTPUT.tmp"
mv "\$OUTPUT.tmp" "\$OUTPUT"
echo "[\$(date +%H:%M:%S)] stage-out done in \$((SECONDS - t0)) s"
echo "[\$(date +%H:%M:%S)] DONE: \$OUTPUT (\$(stat -c %s "\$OUTPUT") bytes)"
EOF
chmod +x "$RUNNER"

SBATCH_ARGS=(
    --cpus-per-task=$T_THREADS
    --mem="${SBATCH_MEM}G"
    --time="${SBATCH_TIME}:00:00"
    --output="${OUTPUT}.slurm_%j.log"
    --error="${OUTPUT}.slurm_%j.err.log"
    -J taffy_index
)
[[ -n "$TMP_GB"    ]] && SBATCH_ARGS+=( --tmp="${TMP_GB}G" )
[[ -n "$PARTITION" ]] && SBATCH_ARGS+=( --partition="$PARTITION" )
[[ -n "$ACCOUNT"   ]] && SBATCH_ARGS+=( --account="$ACCOUNT" )

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo ">> DRY RUN -- would submit:"
    echo "sbatch ${SBATCH_ARGS[*]} --parsable $RUNNER"
    echo
    echo ">> generated runner: $RUNNER"
    echo "   (preserved for inspection; rm when done)"
else
    echo ">> submitting..."
    JOB=$(sbatch "${SBATCH_ARGS[@]}" --parsable "$RUNNER")
    echo ">> job id: $JOB"
    echo ">> runner: $RUNNER"
    echo ">> stdout: ${OUTPUT}.slurm_${JOB}.log"
    echo ">> stderr: ${OUTPUT}.slurm_${JOB}.err.log"
fi

# Block driver until the job finishes if --no-wait wasn't passed.
# Uses squeue polling + sacct final-state check (squeue drops completed
# jobs after a short retention window so squeue alone can race the
# exit).
if [[ "$DRY_RUN" -ne 1 && "$WAIT" -eq 1 ]]; then
    echo ">> waiting for job $JOB ..."
    while squeue -j "$JOB" -h -o "%T" 2>/dev/null | grep -qE "PENDING|RUNNING|CONFIGURING|COMPLETING|RESIZING|SUSPENDED|REQUEUED"; do
        sleep 60
    done
    FINAL_STATE=$(sacct -j "$JOB" -X -n -o State 2>/dev/null | head -1 | tr -d ' ')
    echo ">> job $JOB final state: ${FINAL_STATE:-UNKNOWN}"
    case "$FINAL_STATE" in
        COMPLETED) ;;
        *)         echo ">> NON-SUCCESS state -- check ${OUTPUT}.slurm_${JOB}.err.log"; exit 1;;
    esac
fi

echo ">> done."
