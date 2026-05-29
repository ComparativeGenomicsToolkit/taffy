#!/bin/bash
#
# taffy gerp shard driver -- SLURM
#
# Splits a universal MAF's universal-column space into N equal ranges,
# submits an sbatch array running `taffy gerp --columnRange LO-HI` per
# shard, and optionally concatenates the per-shard wig outputs.
#
# Sharding model
# --------------
# T = `taffy stats -i FILE -u`  (total universal columns from the .tui)
# Shard k of N owns columns [k*T/N, (k+1)*T/N).  Each universal column
# belongs to exactly one shard (tui_extract_iterator clips physical blocks
# at the shard boundary, so a block straddling two shards becomes two
# sub-blocks, with disjoint columns).  No chrom enumeration, no bin
# packing, no double counting from leaf vs. anchor chrom names.
#
# Outputs (under $OUTDIR):
#
#   manifest.tsv                  # shard_id <TAB> col_lo <TAB> col_hi
#   shard_<K>.rs.wig.gz           # per-shard outputs (after the array runs)
#   shard_<K>.depth.wig.gz
#   all.rs.wig.gz                 # final concatenation (if --concat)
#   all.depth.wig.gz
#   slurm_<jobid>_<k>.log         # one log per array task
#
# Usage:
#   gerp_shard_slurm.sh -i UNI.maf.gz -o OUTDIR [options]
#
# Requirements: taffy (with gerp + the .tui index for $INPUT), sbatch.

set -euo pipefail

INPUT=""
OUTDIR=""
N_SHARDS=10
THREADS=10
TIME_HOURS=24
MEM_GB=32
TREE=""
PARTITION=""
ACCOUNT=""
DO_CONCAT=1
DRY_RUN=0
TAFFY="${TAFFY:-$(command -v taffy || true)}"

usage() {
    cat >&2 <<EOF
gerp_shard_slurm.sh -- run \`taffy gerp\` on a universal MAF in N SLURM shards
                       (column-range sharding via the .tui's T)

Required:
  -i FILE       Input universal MAF (.maf.gz with a sibling .tui index)
  -o DIR        Output directory (created if missing)

Optional:
  -n INT        Number of shards (default 10)
  -T INT        Threads per shard (default 10; --cpus-per-task=INT)
  --time HRS    Per-task wall limit in hours (default 24)
  --mem GB      Per-task memory in GB (default 32)
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

mkdir -p "$OUTDIR"
echo ">> output dir: $OUTDIR"
echo ">> taffy:      $TAFFY"
echo ">> input:      $INPUT (.tui present)"
echo ">> shards:     $N_SHARDS,  threads/shard: $THREADS"

# --- Step 1: total universal column count T from the .tui. ------------
echo ">> reading T via \`taffy stats -u\` ..."
T_TOTAL=$("$TAFFY" stats -i "$INPUT" -u)
[[ "$T_TOTAL" =~ ^[0-9]+$ ]] || {
    echo "ERROR: taffy stats -u did not return a positive integer: $T_TOTAL" >&2
    exit 1
}
echo ">> T = $T_TOTAL universal columns"

# --- Step 2: partition T into N equal column ranges. ------------------
# Use python for the integer math + manifest write -- shell arithmetic
# is too fragile at T > 2^32.
MANIFEST="$OUTDIR/manifest.tsv"
python3 - "$T_TOTAL" "$N_SHARDS" "$MANIFEST" <<'PY'
import sys
T, N, dst = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
with open(dst, "w") as out:
    out.write("shard\tcol_lo\tcol_hi\n")
    for k in range(N):
        lo = (k     * T) // N
        hi = ((k + 1) * T) // N
        out.write(f"{k}\t{lo}\t{hi}\n")
print(f"shard size: ~{T // N:,} columns each (last shard +{T - (T // N) * N} for remainder)",
      file=sys.stderr)
PY

# --- Step 3: write the per-shard runner. ------------------------------
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
MANIFEST="$OUTDIR/manifest.tsv"

# Lookup this shard's column range from the manifest (1-indexed, skip header).
read SHARD_ID COL_LO COL_HI < <(awk -v k="\$K" 'NR>1 && \$1==k {print \$1, \$2, \$3; exit}' "\$MANIFEST")
if [[ -z "\${COL_LO:-}" ]]; then
    echo "shard \${K}: no manifest row found" >&2
    exit 1
fi

RS_OUT="\$OUTDIR/shard_\${K}.rs.wig.gz"
DEPTH_OUT="\$OUTDIR/shard_\${K}.depth.wig.gz"

# Idempotency: if both outputs exist + non-empty, skip.
if [[ -s "\$RS_OUT" && -s "\$DEPTH_OUT" ]]; then
    echo "shard \${K}: outputs present, skipping"
    exit 0
fi

# .tmp -> rename so a failed run doesn't leave a half-written .wig.gz
# that future idempotency checks would mistake for completed work.
"\$TAFFY" gerp -i "\$INPUT" --columnRange "\${COL_LO}-\${COL_HI}" -c -T "\$THREADS" \\
    \$TREE_FLAG \\
    -o "\$RS_OUT".tmp -D "\$DEPTH_OUT".tmp
mv "\$RS_OUT".tmp    "\$RS_OUT"
mv "\$DEPTH_OUT".tmp "\$DEPTH_OUT"
echo "shard \${K}: done (columns [\${COL_LO}, \${COL_HI}))"
EOF
chmod +x "$RUNNER"

# --- Step 4: assemble + submit sbatch array. --------------------------
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

# --- Step 5: optional concat job (depends on array success). -----------
# BGZF is append-safe: each shard's .wig.gz ends in an EOF block which is
# also valid mid-stream, and zcat / htslib read concatenated bgzip
# transparently.  The concat preserves universal-column order if shards
# emit in order, which they do (each shard's columns are disjoint from
# the others' and we cat in shard-id order).
if [[ "$DO_CONCAT" -eq 1 ]]; then
    CONCAT_SCRIPT="$OUTDIR/concat.sh"
    cat > "$CONCAT_SCRIPT" <<EOF
#!/bin/bash
set -euo pipefail
cd "$OUTDIR"
# Sort by numeric shard id so concat order matches column order.
ls -1 shard_*.rs.wig.gz    | sort -t_ -k2,2n | xargs cat > all.rs.wig.gz.tmp    && mv all.rs.wig.gz.tmp    all.rs.wig.gz
ls -1 shard_*.depth.wig.gz | sort -t_ -k2,2n | xargs cat > all.depth.wig.gz.tmp && mv all.depth.wig.gz.tmp all.depth.wig.gz
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
