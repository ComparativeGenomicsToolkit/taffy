#!/bin/bash
#
# Universal-depth bigWig for a BIG universal MAF (e.g. 577-way) -- SLURM
# =====================================================================
# Shards universal-column space into ranges, runs one
#   taffy depth --coords named --depth --bin BIN --columnRange LO-HI
# per shard (sbatch array), then a dependent job sorts + merges the per-shard
# records into the NAMED row-0 (ancestor) depth bigWig.
#
# Named coords (the only option past 2^31 columns) key each record on the row-0
# ancestor seq, which recurs through the column order -- so a (name,bin) can be
# split across shard boundaries.  Shards therefore need NO bin alignment; the
# assembly's sort + per-(name,bin) merge (uni_depth_merge.awk) reconciles them.
#
# For a small alignment (<= a fish subtree) use uni_depth_bigwig_local.sh.
#
# Usage: uni_depth_bigwig_slurm.sh -i UNI.maf.gz -o OUT.bw [options]
#   -i FILE        universal MAF/TAF (its <FILE>.tui must already exist)
#   -o FILE        output bigWig (shards + .bg + .sizes land beside it)
#   -b INT         bin width in bp (default 1000)
#   -s INT         shard size in universal columns (default 1e8)
#   -T INT         depth threads per shard (default 8)
#   --mem GB       sbatch --mem per shard task (default 16)
#   --time HRS     sbatch --time per shard task (default 4)
#   --partition P  --account A    passed through to sbatch
#   --dry-run      write the job scripts and print sbatch lines, don't submit
# Env: TAFFY, WIGTOBIGWIG override the binaries.
#
# Requirements: taffy (`depth` subcommand), wigToBigWig (UCSC), sbatch.

set -euo pipefail

INPUT=""; OUT=""; BIN=1000; SHARD=100000000; THREADS=8
MEM_GB=16; TIME_HRS=4; PARTITION=""; ACCOUNT=""; DRYRUN=0
TAFFY="${TAFFY:-taffy}"; WIGTOBIGWIG="${WIGTOBIGWIG:-wigToBigWig}"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() { sed -n '17,27p' "$0" >&2; exit "${1:-1}"; }
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
[ -e "$INPUT.tui" ]   || { echo "ERROR: $INPUT.tui missing -- build it first (taffy index -u)" >&2; exit 1; }
command -v "$TAFFY" >/dev/null 2>&1 || { echo "ERROR: taffy not found (set \$TAFFY)" >&2; exit 1; }
command -v "$WIGTOBIGWIG" >/dev/null 2>&1 || echo "WARN: wigToBigWig not on this PATH; the assembly job needs it on the compute node" >&2
(( BIN > 0 ))   || { echo "ERROR: -b ($BIN) must be > 0" >&2; exit 1; }
(( SHARD > 0 )) || { echo "ERROR: -s ($SHARD) must be > 0" >&2; exit 1; }

T=$("$TAFFY" stats -i "$INPUT" -u)
NSHARD=$(( (T + SHARD - 1) / SHARD ))
OUTDIR=$(dirname "$OUT"); mkdir -p "$OUTDIR"
BASE=$(basename "$OUT"); BASE=${BASE%.bw}; BASE=${BASE%.bigWig}
SHARDDIR="$OUTDIR/$BASE.shards"; mkdir -p "$SHARDDIR"
RUNNER="$OUTDIR/.$BASE.runner.sh"
ASSEMBLE="$OUTDIR/.$BASE.assemble.sh"
echo ">> T=$T  bin=$BIN bp  shard=$SHARD cols  shards=$NSHARD  (named row-0 coords)" >&2

# ---- per-shard depth ----
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
"$TAFFY" depth -i "$INPUT" --coords named --bin $BIN --minLeaves 1 \\
  -T \$SLURM_CPUS_PER_TASK --columnRange \$LO-\$HI \\
  --depth "$SHARDDIR/shard.\$(printf '%05d' \$SLURM_ARRAY_TASK_ID).bg"
EOF

# ---- assembly (waits for the whole array) ----
cat > "$ASSEMBLE" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=8:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
# Named coords: a (row-0 seq, bin) can be split across shard boundaries, so sort
# by (name,start) then sum (sum,cnt) per (name,bin) into a mean-depth bedGraph.
# Chrom sizes (ancestor refChrs + leaves) come from the .tui via stats -s.
cat "$SHARDDIR"/shard.*.bg | LC_ALL=C sort -k1,1 -k2,2n -T "$OUTDIR" | awk -v N=$BIN -f "$SCRIPTDIR/uni_depth_merge.awk" > "$OUTDIR/$BASE.bg"
"$TAFFY" stats -i "$INPUT" -s > "$OUTDIR/$BASE.sizes"
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
