#!/bin/bash
#
# taffy gerp shard driver -- SLURM
#
# Splits a universal MAF (with .tui) into N roughly-equal-bp shards by
# anchor chrom, submits an sbatch array that runs `taffy gerp -r` on each
# shard, and optionally concatenates the per-shard wig outputs.
#
# Outputs (under $OUTDIR):
#
#   chrom_sizes.tsv               # raw `taffy stats -s` output
#   manifest.tsv                  # shard_id <TAB> chrom <TAB> start <TAB> end
#   regions/shard_<K>.regions     # one "chrom:start-end" per line
#   shard_<K>.rs.wig.gz           # per-shard outputs (after the array runs)
#   shard_<K>.depth.wig.gz
#   all.rs.wig.gz                 # final concatenation (if --concat)
#   all.depth.wig.gz
#   slurm_<jobid>_<k>.log         # one log per array task
#
# Usage:
#   gerp_shard_slurm.sh -i UNI.maf.gz -o OUTDIR [options]
#
# Requirements: taffy (with gerp + stats + the .tui index for $INPUT),
# python3, sbatch.

set -euo pipefail

INPUT=""
OUTDIR=""
N_SHARDS=10
THREADS=10
TIME_HOURS=24
MEM_GB=32
MIN_CHROM_BP=0    # drop chroms smaller than this (recommended >0 on big inputs)
TREE=""
PARTITION=""
ACCOUNT=""
DO_CONCAT=1
DRY_RUN=0
TAFFY="${TAFFY:-$(command -v taffy || true)}"

usage() {
    cat >&2 <<EOF
gerp_shard_slurm.sh -- run \`taffy gerp\` on a universal MAF in N SLURM shards

Required:
  -i FILE       Input universal MAF (.maf.gz with a sibling .tui index)
  -o DIR        Output directory (created if missing)

Optional:
  -n INT        Number of shards (default 10)
  -T INT        Threads per shard (default 10; --cpus-per-task=INT)
  --time HRS    Per-task wall limit in hours (default 24)
  --mem GB      Per-task memory in GB (default 32)
  --min-chrom-bp N
                Drop chroms smaller than N bp from the manifest entirely.
                Default 0 (include all).  Useful for trimming long-tail
                tiny anchor chroms whose alignment contribution is
                negligible.  With \`taffy gerp -R\` (the runner uses one
                invocation per shard, so tui_load is amortised) the main
                remaining per-region cost is the tui_query + iterator
                setup.  10000 is a reasonable starting point on
                vertebrate-scale alignments; 0 is fine if you want
                whole-genome coverage.
  --tree FILE   Newick tree override (default: \`# hal\` from MAF header)
  --partition X SLURM partition (--partition=X)
  --account X   SLURM account (--account=X)
  --no-concat   Skip the post-array concatenation job
  --dry-run     Print the sbatch commands but do not submit
  -h            Help
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)            INPUT="$2"; shift 2;;
        -o)            OUTDIR="$2"; shift 2;;
        -n)            N_SHARDS="$2"; shift 2;;
        -T)            THREADS="$2"; shift 2;;
        --time)        TIME_HOURS="$2"; shift 2;;
        --mem)         MEM_GB="$2"; shift 2;;
        --min-chrom-bp) MIN_CHROM_BP="$2"; shift 2;;
        --tree)        TREE="$2"; shift 2;;
        --partition)   PARTITION="$2"; shift 2;;
        --account)     ACCOUNT="$2"; shift 2;;
        --no-concat)   DO_CONCAT=0; shift;;
        --dry-run)     DRY_RUN=1; shift;;
        -h|--help)     usage 0;;
        *)             echo "unknown arg: $1" >&2; usage 1;;
    esac
done

[[ -n "$INPUT"  ]] || { echo "ERROR: -i required" >&2; usage 1; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: -o required" >&2; usage 1; }
[[ -n "$TAFFY"  ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -f "$INPUT"      ]] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
[[ -f "${INPUT}.tui" ]] || { echo "ERROR: $INPUT has no .tui sibling -- run \`taffy index -u\` first" >&2; exit 1; }

mkdir -p "$OUTDIR" "$OUTDIR/regions"
echo ">> output dir: $OUTDIR"
echo ">> taffy:      $TAFFY"
echo ">> input:      $INPUT (.tui present)"
echo ">> shards:     $N_SHARDS,  threads/shard: $THREADS"

# --- Step 1: list anchor chroms + sizes. -----------------------------
CHROM_TSV="$OUTDIR/chrom_sizes.tsv"
if [[ ! -s "$CHROM_TSV" ]]; then
    echo ">> listing chrom sizes via \`taffy stats -s\` ..."
    "$TAFFY" stats -i "$INPUT" -s > "$CHROM_TSV"
fi
N_CHROMS=$(wc -l < "$CHROM_TSV")
TOTAL_BP=$(awk '{s+=$2} END{print s+0}' "$CHROM_TSV")
echo ">> $N_CHROMS chroms, $TOTAL_BP total bp"

# --- Step 2: partition via LPT (Longest Processing Time first). -------
# Sort chroms by size desc, greedily assign each to the currently-smallest
# shard.  Near-optimal balance; O(N log K) with a heap (we use list-min
# since N_SHARDS is small).
MANIFEST="$OUTDIR/manifest.tsv"
python3 - "$CHROM_TSV" "$N_SHARDS" "$MANIFEST" "$MIN_CHROM_BP" <<'PY'
import sys
import heapq
src, n_shards, dst, min_bp = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
chroms = []
n_dropped = 0
dropped_bp = 0
with open(src) as f:
    for line in f:
        parts = line.split()
        if len(parts) < 2: continue
        sz = int(parts[1])
        if sz < min_bp:
            n_dropped += 1
            dropped_bp += sz
            continue
        chroms.append((parts[0], sz))
chroms.sort(key=lambda kv: -kv[1])  # biggest first
heap = [(0, k) for k in range(n_shards)]
heapq.heapify(heap)
assignments = [[] for _ in range(n_shards)]
for c, sz in chroms:
    cur_size, k = heapq.heappop(heap)
    assignments[k].append((c, sz))
    heapq.heappush(heap, (cur_size + sz, k))
with open(dst, "w") as out:
    out.write("shard\tchrom\tstart\tend\n")
    for k, items in enumerate(assignments):
        for c, sz in items:
            out.write(f"{k}\t{c}\t0\t{sz}\n")
shard_sizes  = [sum(sz for _, sz in items) for items in assignments]
shard_counts = [len(items) for items in assignments]
print(f"shard bp range: min={min(shard_sizes):,}  max={max(shard_sizes):,}  "
      f"max/min={max(shard_sizes)/max(min(shard_sizes),1):.2f}", file=sys.stderr)
print(f"regions per shard: min={min(shard_counts)}  max={max(shard_counts)}",
      file=sys.stderr)
if n_dropped > 0:
    print(f"dropped {n_dropped} chroms below --min-chrom-bp ({dropped_bp:,} bp, "
          f"{100*dropped_bp/(dropped_bp+sum(shard_sizes)):.3f}% of total)",
          file=sys.stderr)
PY

# --- Step 3: emit one regions-file per shard. -------------------------
rm -f "$OUTDIR"/regions/shard_*.regions
awk -v dir="$OUTDIR/regions" 'NR>1 {print $2":"$3"-"$4 >> (dir "/shard_" $1 ".regions")}' \
    "$MANIFEST"

# --- Step 4: write the per-shard runner. ------------------------------
# One `taffy gerp -R` invocation per shard -- the .tui loads ONCE for the
# shard instead of once per region, which is the whole point of -R.
RUNNER="$OUTDIR/run_shard.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -euo pipefail
K=\${SLURM_ARRAY_TASK_ID:-0}
INPUT="$INPUT"
OUTDIR="$OUTDIR"
THREADS="$THREADS"
TAFFY="$TAFFY"
TREE_FLAG=$([[ -n "$TREE" ]] && echo "-t $TREE" || echo "")

REGIONS="\$OUTDIR/regions/shard_\${K}.regions"
RS_OUT="\$OUTDIR/shard_\${K}.rs.wig.gz"
DEPTH_OUT="\$OUTDIR/shard_\${K}.depth.wig.gz"

# Idempotency: if both outputs exist + non-empty, skip.
if [[ -s "\$RS_OUT" && -s "\$DEPTH_OUT" ]]; then
    echo "shard \${K}: outputs present, skipping"
    exit 0
fi

# Write to .tmp then rename so a failed run doesn't leave a half-written
# .wig.gz that future idempotency checks would mistake for completed work.
"\$TAFFY" gerp -i "\$INPUT" -R "\$REGIONS" -c -T "\$THREADS" \\
    \$TREE_FLAG \\
    -o "\$RS_OUT".tmp -D "\$DEPTH_OUT".tmp
mv "\$RS_OUT".tmp    "\$RS_OUT"
mv "\$DEPTH_OUT".tmp "\$DEPTH_OUT"
echo "shard \${K}: done (\$(wc -l < \$REGIONS) regions)"
EOF
chmod +x "$RUNNER"

# --- Step 5: assemble + submit sbatch array. --------------------------
SBATCH_ARGS=(
    --array=0-$((N_SHARDS - 1))
    --cpus-per-task="$THREADS"
    --mem="${MEM_GB}G"
    --time="${TIME_HOURS}:00:00"
    --output="$OUTDIR/slurm_%A_%a.log"
    -J taffy_gerp_shard
)
[[ -n "$PARTITION" ]] && SBATCH_ARGS+=( --partition="$PARTITION" )
[[ -n "$ACCOUNT"   ]] && SBATCH_ARGS+=( --account="$ACCOUNT" )

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo ">> DRY RUN -- would submit:"
    echo "sbatch ${SBATCH_ARGS[*]} --parsable $RUNNER"
    ARRAY_JOB=DRY
else
    echo ">> submitting array..."
    ARRAY_JOB=$(sbatch "${SBATCH_ARGS[@]}" --parsable "$RUNNER")
    echo ">> array job id: $ARRAY_JOB"
fi

# --- Step 6: optional concat job (depends on array success). -----------
if [[ "$DO_CONCAT" -eq 1 ]]; then
    CONCAT_SCRIPT="$OUTDIR/concat.sh"
    cat > "$CONCAT_SCRIPT" <<EOF
#!/bin/bash
set -euo pipefail
cd "$OUTDIR"
cat shard_*.rs.wig.gz    > all.rs.wig.gz.tmp    && mv all.rs.wig.gz.tmp    all.rs.wig.gz
cat shard_*.depth.wig.gz > all.depth.wig.gz.tmp && mv all.depth.wig.gz.tmp all.depth.wig.gz
echo "concat done: \$(ls -la all.*.wig.gz)"
EOF
    chmod +x "$CONCAT_SCRIPT"

    CONCAT_ARGS=(
        --dependency=afterok:"$ARRAY_JOB"
        --cpus-per-task=2
        --mem=4G
        --time=1:00:00
        --output="$OUTDIR/slurm_concat_%j.log"
        -J taffy_gerp_concat
    )
    [[ -n "$PARTITION" ]] && CONCAT_ARGS+=( --partition="$PARTITION" )
    [[ -n "$ACCOUNT"   ]] && CONCAT_ARGS+=( --account="$ACCOUNT" )

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo ">> DRY RUN -- would submit:"
        echo "sbatch ${CONCAT_ARGS[*]} $CONCAT_SCRIPT"
    else
        CONCAT_JOB=$(sbatch "${CONCAT_ARGS[@]}" --parsable "$CONCAT_SCRIPT")
        echo ">> concat job id: $CONCAT_JOB"
    fi
fi

echo ">> done."
