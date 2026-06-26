#!/bin/bash
#
# Universal per-species depth bigWig for a BIG universal MAF (e.g. 577-way) -- SLURM
# =================================================================================
# Shards universal-column space into ranges, runs one
#   taffy depth --bin BIN --minLeaves 1 --columnRange LO-HI --perSpecies --bigwig
# per shard (sbatch array, each emitting a 64-bit per-species VECTOR bigWig slice
# on the single uni0 axis), then a dependent job MERGES the per-shard bigWigs into
# one with `taffy merge-bigwig`.
#
# 64-bit axis: the raw column [0,T) lives on a single chrom uni0 (no more 2e9
# uni0..uniK chunking, no wigToBigWig).  Each shard is a disjoint, bin-aligned
# column slice; bigWig is binary so the assembly MERGES (re-streams + re-indexes),
# not `cat`.  The output carries an <out>.names sidecar (one leaf per component).
#
# For a small alignment (<= a fish subtree) use uni_depth_bigwig_local.sh.
#
# Usage: uni_depth_bigwig_slurm.sh -i UNI.maf.gz -o OUT.bw [options]
#   -i FILE        universal MAF/TAF (its <FILE>.tui must already exist)
#   -o FILE        output bigWig; shards + .names + .sizes land beside it, so point
#                  it at a ROOMY fs (the full 577 per-species .bw is a few-to-tens
#                  of GB), NOT a near-full /home
#   -b INT         bin width in bp (default 1000)
#   -s INT         shard size in universal columns (default 1e8; must be a
#                  multiple of BIN so shard edges fall on bin boundaries)
#   -T INT         depth threads per shard (default 8)
#   --mem GB       sbatch --mem per shard task (default 24; per-species N-wide
#                  buffers + gerp scratch -- raise for the 577 backbone)
#   --time HRS     sbatch --time per shard task (default 12; raise for dense
#                  577 backbone shards -- a TIMEOUT strands the whole build)
#   --partition P  --account A    passed through to sbatch
#   --dry-run      write the job scripts and print sbatch lines, don't submit
#   --wait         submit, then BLOCK until the build finishes (poll the merge job)
# Env: TAFFY overrides the binary.
#
# Requirements: taffy (`depth` + `merge-bigwig` subcommands), sbatch.

set -euo pipefail

INPUT=""; OUT=""; BIN=1000; SHARD=100000000; THREADS=8
MEM_GB=24; TIME_HRS=12; PARTITION=""; ACCOUNT=""; DRYRUN=0; WAIT=0
TAFFY="${TAFFY:-taffy}"

usage() { sed -n '20,38p' "$0" >&2; exit "${1:-1}"; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT=$2; shift 2;;
    -o) OUT=$2; shift 2;;
    -b) BIN=$2; shift 2;;
    -s) SHARD=$2; shift 2;;
    -T) THREADS=$2; shift 2;;
    --mem) MEM_GB=$2; shift 2;;
    --time) TIME_HRS=$2; shift 2;;
    --partition) PARTITION=$2; shift 2;;
    --account) ACCOUNT=$2; shift 2;;
    --dry-run) DRYRUN=1; shift;;
    --wait) WAIT=1; shift;;
    -h|--help) usage 0;;
    *) echo "ERROR: unknown arg: $1" >&2; usage;;
  esac
done

[ -n "$INPUT" ] && [ -n "$OUT" ] || usage
[ -e "$INPUT" ]       || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
[ -e "$INPUT.tui" ]   || { echo "ERROR: $INPUT.tui missing -- build it first (taffy index -u)" >&2; exit 1; }
command -v "$TAFFY" >/dev/null 2>&1 || { echo "ERROR: taffy not found (set \$TAFFY)" >&2; exit 1; }
(( BIN > 0 ))   || { echo "ERROR: -b ($BIN) must be > 0" >&2; exit 1; }
(( SHARD > 0 )) || { echo "ERROR: -s ($SHARD) must be > 0" >&2; exit 1; }
(( SHARD % BIN == 0 )) || { echo "ERROR: -s ($SHARD) must be a multiple of -b ($BIN) so shard edges fall on bin boundaries" >&2; exit 1; }

T=$("$TAFFY" stats -i "$INPUT" -u)
NSHARD=$(( (T + SHARD - 1) / SHARD ))
OUTDIR=$(dirname "$OUT"); mkdir -p "$OUTDIR"
BASE=$(basename "$OUT"); BASE=${BASE%.bw}; BASE=${BASE%.bigWig}
SHARDDIR="$OUTDIR/$BASE.shards"; mkdir -p "$SHARDDIR"
RUNNER="$OUTDIR/.$BASE.runner.sh"
ASSEMBLE="$OUTDIR/.$BASE.assemble.sh"
echo ">> T=$T  bin=$BIN bp  shard=$SHARD cols  shards=$NSHARD  (per-species, single uni0 axis)" >&2

# ---- per-shard per-species depth bigWig ----
# Shard 0 also writes the global chrom-sizes (--sizes; a single `uni0<TAB>T` line,
# a pure function of T) for reference; the merge reads T from the shard headers.
cat > "$RUNNER" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=${MEM_GB}G
#SBATCH --time=${TIME_HRS}:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
LO=\$(( SLURM_ARRAY_TASK_ID * $SHARD ))
HI=\$(( LO + $SHARD )); (( HI > $T )) && HI=$T
(( LO >= HI )) && exit 0
SIZES_ARG=""
(( SLURM_ARRAY_TASK_ID == 0 )) && SIZES_ARG="--sizes $OUTDIR/$BASE.sizes"
"$TAFFY" depth -i "$INPUT" --bin $BIN --minLeaves 1 \\
  -T \$SLURM_CPUS_PER_TASK --columnRange \$LO-\$HI \$SIZES_ARG \\
  --perSpecies --bigwig "$SHARDDIR/shard.\$(printf '%05d' \$SLURM_ARRAY_TASK_ID).bw"
EOF

# ---- assembly: merge the per-shard bigWigs (waits for the whole array) ----
cat > "$ASSEMBLE" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=8:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
# Each shard is a disjoint uni0 column slice; the zero-padded shard filenames sort
# in column order, so the glob is already in order.  merge-bigwig re-streams them
# into one bigWig (rebuilding the R-tree index + zoom) and writes $OUT.names.
"$TAFFY" merge-bigwig -o "$OUT" "$SHARDDIR"/shard.*.bw
echo "made $OUT (\$(du -h "$OUT"|cut -f1)); $OUT.names (\$(wc -l < "$OUT.names") components)"
EOF

if (( DRYRUN )); then
  echo "[dry-run] wrote $RUNNER and $ASSEMBLE" >&2
  echo "[dry-run] sbatch --array=0-$((NSHARD-1)) $RUNNER" >&2
  echo "[dry-run] sbatch --dependency=afterok:<arrayid> $ASSEMBLE" >&2
  exit 0
fi
ARRAY_ID=$(sbatch --parsable --array=0-$((NSHARD-1)) "$RUNNER")
echo ">> depth array submitted: job $ARRAY_ID ($NSHARD shards)" >&2
ASM_ID=$(sbatch --parsable --dependency=afterok:"$ARRAY_ID" "$ASSEMBLE")
echo ">> merge submitted (afterok:$ARRAY_ID): job $ASM_ID -> $OUT" >&2

if (( WAIT )); then
  # Block until the merge job leaves the queue.  These are real SLURM jobs (array
  # + afterok merge), so the build still finishes if THIS poll process dies -- run
  # it under nohup/tmux for a multi-hour build; --wait just adds the blocking.  If
  # a shard fails, afterok cancels the merge -> it leaves the queue and $OUT is
  # never written, so the -s check reports the failure.
  echo ">> [--wait] blocking until job $ASM_ID completes (polling every 30s)..." >&2
  sleep 5
  while [[ -n "$(squeue -h -j "$ASM_ID" 2>/dev/null)" ]]; do sleep 30; done
  if [[ -s "$OUT" ]]; then
    echo ">> done -> $OUT ($(du -h "$OUT" 2>/dev/null | cut -f1))" >&2
  else
    echo ">> BUILD FAILED: $OUT not written -- inspect the depth array ($ARRAY_ID) and merge ($ASM_ID)" >&2
    exit 1
  fi
fi
