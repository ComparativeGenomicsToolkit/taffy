#!/bin/bash
#
# Universal-depth bigWig for a BIG universal MAF (e.g. 577-way) -- SLURM
# =====================================================================
# Shards universal-column space into BIN-aligned ranges, runs one
#   taffy gerp --universal --depthOnly --bin BIN --columnRange LO-HI
# per shard (sbatch array), then a dependent job concatenates the per-shard
# bedGraphs and wigToBigWig's them into the uni<chunk> depth bigWig that
# `taffy lift --bigwig` queries.
#
# Why this and not gerp_shard_slurm.sh: that one emits the per-column RS+depth
# wigs (named ancestor coords).  This one emits the BINNED, integer-universal-
# column (uni<chunk>) depth -- the value source for the query shim.  Shard
# boundaries are multiples of BIN, so gerp's --bin alignment guard accepts every
# shard, no bin is split across shards, and the per-shard bedGraphs concatenate
# in column order = already sorted (wigToBigWig needs no external sort).
#
# For a small alignment (≲ a fish subtree) use uni_depth_bigwig_local.sh.
#
# Usage: uni_depth_bigwig_slurm.sh -i UNI.maf.gz -o OUT.bw [options]
#   -i FILE        universal MAF/TAF (its <FILE>.tui must already exist)
#   -o FILE        output bigWig (shards + .bg + .uni.sizes land beside it)
#   -b INT         bin width in universal columns; must divide 4e9 (default 1000)
#   -s INT         shard size in columns; must be a multiple of -b (default 1e8)
#   -T INT         gerp threads per shard (default 8)
#   --mem GB       sbatch --mem per shard task (default 16)
#   --time HRS     sbatch --time per shard task (default 4)
#   --partition P  --account A    passed through to sbatch
#   --dry-run      write the job scripts and print sbatch lines, don't submit
# Env: TAFFY, WIGTOBIGWIG override the binaries.
#
# Requirements: taffy (gerp --bin, commit 8f9c5d7+), wigToBigWig (UCSC), sbatch.

set -euo pipefail

INPUT=""; OUT=""; BIN=1000; SHARD=100000000; THREADS=8
MEM_GB=16; TIME_HRS=4; PARTITION=""; ACCOUNT=""; DRYRUN=0
TAFFY="${TAFFY:-taffy}"; WIGTOBIGWIG="${WIGTOBIGWIG:-wigToBigWig}"

usage() { sed -n '17,30p' "$0" >&2; exit "${1:-1}"; }
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
    -h|--help) usage 0;;
    *) echo "ERROR: unknown arg: $1" >&2; usage;;
  esac
done

[ -n "$INPUT" ] && [ -n "$OUT" ] || usage
[ -e "$INPUT" ]       || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
[ -e "$INPUT.tui" ]   || { echo "ERROR: $INPUT.tui missing -- build it first (taffy index -u, or taffy_index_slurm.sh)" >&2; exit 1; }
command -v "$TAFFY" >/dev/null 2>&1 || { echo "ERROR: taffy not found (set \$TAFFY)" >&2; exit 1; }
command -v "$WIGTOBIGWIG" >/dev/null 2>&1 || echo "WARN: wigToBigWig not on this PATH; the assembly job needs it on the compute node" >&2
(( 4000000000 % BIN == 0 )) || { echo "ERROR: -b ($BIN) must divide 4000000000" >&2; exit 1; }
(( SHARD % BIN == 0 ))      || { echo "ERROR: -s ($SHARD) must be a multiple of -b ($BIN)" >&2; exit 1; }

T=$("$TAFFY" stats -i "$INPUT" -u)
TA=$(( T / BIN * BIN ))                 # bin-aligned end (drops the final <BIN cols)
NSHARD=$(( (TA + SHARD - 1) / SHARD ))
OUTDIR=$(dirname "$OUT"); mkdir -p "$OUTDIR"
BASE=$(basename "$OUT"); BASE=${BASE%.bw}; BASE=${BASE%.bigWig}
SHARDDIR="$OUTDIR/$BASE.shards"; mkdir -p "$SHARDDIR"
RUNNER="$OUTDIR/.$BASE.runner.sh"
ASSEMBLE="$OUTDIR/.$BASE.assemble.sh"
echo ">> T=$T  bin-aligned end=$TA  bin=$BIN  shard=$SHARD cols  shards=$NSHARD" >&2

# ---- per-shard gerp ----
cat > "$RUNNER" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=${MEM_GB}G
#SBATCH --time=${TIME_HRS}:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
LO=\$(( SLURM_ARRAY_TASK_ID * $SHARD ))
HI=\$(( LO + $SHARD )); (( HI > $TA )) && HI=$TA
(( LO >= HI )) && exit 0
"$TAFFY" gerp -i "$INPUT" --universal --depthOnly --bin $BIN --minLeaves 1 \\
  -T \$SLURM_CPUS_PER_TASK --columnRange \$LO-\$HI \\
  -o /dev/null -D "$SHARDDIR/shard.\$(printf '%05d' \$SLURM_ARRAY_TASK_ID).bg"
EOF

# ---- assembly (waits for the whole array) ----
cat > "$ASSEMBLE" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
cat "$SHARDDIR"/shard.*.bg > "$OUTDIR/$BASE.bg"
awk -v T=$T -v C=4000000000 'BEGIN{for(c=0;c*C<T;c++){s=((c+1)*C<=T)?C:T-c*C;print "uni"c"\t"s}}' > "$OUTDIR/$BASE.uni.sizes"
"$WIGTOBIGWIG" "$OUTDIR/$BASE.bg" "$OUTDIR/$BASE.uni.sizes" "$OUT"
echo "made $OUT (\$(du -h "$OUT"|cut -f1))"
EOF

if (( DRYRUN )); then
  echo "[dry-run] wrote $RUNNER and $ASSEMBLE" >&2
  echo "[dry-run] sbatch --array=0-$((NSHARD-1)) $RUNNER" >&2
  echo "[dry-run] sbatch --dependency=afterok:<arrayid> $ASSEMBLE" >&2
  exit 0
fi
ARRAY_ID=$(sbatch --parsable --array=0-$((NSHARD-1)) "$RUNNER")
echo ">> gerp array submitted: job $ARRAY_ID ($NSHARD shards)" >&2
ASM_ID=$(sbatch --parsable --dependency=afterok:"$ARRAY_ID" "$ASSEMBLE")
echo ">> assembly submitted (afterok:$ARRAY_ID): job $ASM_ID -> $OUT" >&2
