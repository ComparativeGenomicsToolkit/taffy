#!/bin/bash
#
# Universal-depth bigWig for a BIG universal MAF (e.g. 577-way) -- SLURM
# =====================================================================
# Shards universal-column space into ranges, runs one
#   taffy depth --depth --bin BIN --columnRange LO-HI
# per shard (sbatch array), then a dependent job concatenates the per-shard
# bedGraphs into the universal-column depth bigWig.
#
# The output is on the UNIVERSAL-COLUMN axis: the raw column [0,T) 2e9-chunked
# into chroms uni0,uni1,... (see TUI_UNI_CHUNK).  Each shard emits an already-
# sorted, monotone, bigWig-ready bedGraph slice; a bin never straddles a shard
# boundary (shards are bin-aligned) nor a chunk boundary (2e9 is a multiple of
# BIN), so assembly is just `cat` in shard order -- NO sort, NO merge.
#
# The chrom-sizes file (uni0..uniK, each 2e9 except the last = T mod 2e9) is a
# pure function of T, so shard 0 emits it via `taffy depth --sizes` and the
# assembly reuses it.
#
# For a small alignment (<= a fish subtree) use uni_depth_bigwig_local.sh.
#
# Usage: uni_depth_bigwig_slurm.sh -i UNI.maf.gz -o OUT.bw [options]
#   -i FILE        universal MAF/TAF (its <FILE>.tui must already exist)
#   -o FILE        output bigWig; shards + .bg + .sizes land beside it, so point
#                  it at a ROOMY fs (~9GB peak for the 577), NOT a near-full /home
#   -b INT         bin width in bp (default 1000; must divide 2e9)
#   -s INT         shard size in universal columns (default 1e8; must be a
#                  multiple of BIN so shard edges fall on bin boundaries)
#   -T INT         depth threads per shard (default 8)
#   --mem GB       sbatch --mem per shard task (default 16)
#   --time HRS     sbatch --time per shard task (default 8; raise to 12+ for the
#                  dense 577 backbone shards -- a TIMEOUT strands the whole build)
#   --partition P  --account A    passed through to sbatch
#   --dry-run      write the job scripts and print sbatch lines, don't submit
#   --wait         submit, then BLOCK until the build finishes (poll the assembly
#                  job); else submit array + afterok assembly and return at once
# Env: TAFFY, WIGTOBIGWIG override the binaries.
#
# Requirements: taffy (`depth` subcommand), wigToBigWig (UCSC), sbatch.

set -euo pipefail

INPUT=""; OUT=""; BIN=1000; SHARD=100000000; THREADS=8
MEM_GB=16; TIME_HRS=8; PARTITION=""; ACCOUNT=""; DRYRUN=0; WAIT=0
TAFFY="${TAFFY:-taffy}"; WIGTOBIGWIG="${WIGTOBIGWIG:-wigToBigWig}"

usage() { sed -n '20,37p' "$0" >&2; exit "${1:-1}"; }
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
command -v "$WIGTOBIGWIG" >/dev/null 2>&1 || echo "WARN: wigToBigWig not on this PATH; the assembly job needs it on the compute node" >&2
(( BIN > 0 ))   || { echo "ERROR: -b ($BIN) must be > 0" >&2; exit 1; }
(( 2000000000 % BIN == 0 )) || { echo "ERROR: -b ($BIN) must divide 2000000000 (the uni-axis chunk size)" >&2; exit 1; }
(( SHARD > 0 )) || { echo "ERROR: -s ($SHARD) must be > 0" >&2; exit 1; }
(( SHARD % BIN == 0 )) || { echo "ERROR: -s ($SHARD) must be a multiple of -b ($BIN) so shard edges fall on bin boundaries" >&2; exit 1; }

T=$("$TAFFY" stats -i "$INPUT" -u)
NSHARD=$(( (T + SHARD - 1) / SHARD ))
OUTDIR=$(dirname "$OUT"); mkdir -p "$OUTDIR"
BASE=$(basename "$OUT"); BASE=${BASE%.bw}; BASE=${BASE%.bigWig}
SHARDDIR="$OUTDIR/$BASE.shards"; mkdir -p "$SHARDDIR"
RUNNER="$OUTDIR/.$BASE.runner.sh"
ASSEMBLE="$OUTDIR/.$BASE.assemble.sh"
echo ">> T=$T  bin=$BIN bp  shard=$SHARD cols  shards=$NSHARD  (universal-column uni axis)" >&2

# ---- per-shard depth ----
# Shard 0 also writes the global uni0..uniK chrom-sizes (--sizes; a pure function
# of T, independent of the shard's column range), which the assembly reuses.
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
  --depth "$SHARDDIR/shard.\$(printf '%05d' \$SLURM_ARRAY_TASK_ID).bg"
EOF

# ---- assembly (waits for the whole array) ----
cat > "$ASSEMBLE" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=8:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
# Universal-column axis: each shard is its own column slice, already monotone and
# bin/chunk-aligned, so just concatenate in shard (column) order -- NO data sort,
# NO merge.  Shard files are zero-padded by column index, so the glob expands in
# column order.  The uni0..uniK chrom sizes were written by shard 0 (--sizes).
cat "$SHARDDIR"/shard.*.bg > "$OUTDIR/$BASE.bg"
"$WIGTOBIGWIG" "$OUTDIR/$BASE.bg" "$OUTDIR/$BASE.sizes" "$OUT"
echo "made $OUT (\$(du -h "$OUT"|cut -f1))"
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
echo ">> assembly submitted (afterok:$ARRAY_ID): job $ASM_ID -> $OUT" >&2

if (( WAIT )); then
  # Block until the assembly job leaves the queue.  These are real SLURM jobs
  # (array + afterok assembly), so the build still finishes if THIS poll process
  # dies -- run it under nohup/tmux for a multi-hour build; --wait just adds the
  # blocking.  If a shard fails, afterok cancels the assembly -> it leaves the
  # queue and $OUT is never written, so the -s check reports the failure.
  echo ">> [--wait] blocking until job $ASM_ID completes (polling every 30s)..." >&2
  sleep 5
  while [[ -n "$(squeue -h -j "$ASM_ID" 2>/dev/null)" ]]; do sleep 30; done
  if [[ -s "$OUT" ]]; then
    echo ">> done -> $OUT ($(du -h "$OUT" 2>/dev/null | cut -f1))" >&2
  else
    echo ">> BUILD FAILED: $OUT not written -- inspect the depth array ($ARRAY_ID) and assembly ($ASM_ID)" >&2
    exit 1
  fi
fi
