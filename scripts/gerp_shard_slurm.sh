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
# Local-stage mode (default ON)
# -----------------------------
# At vertebrate scale the .tui is a 92 GB random-access index against a
# 1.4 TB .maf.gz; the per-column query pattern hammers the network FS.
# By default each task copies INPUT + INPUT.tui to $TMPDIR up-front (one
# sequential bulk read per task) then runs gerp against the local copies.
# Use fewer shards with more cores each to amortise the stage-in cost
# (e.g. -n 8 -T 32 instead of -n 32 -T 8).  --no-stage-local reads from
# the original network path -- only sensible for small inputs or local
# filesystems.
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
TMP_GB=""            # passed to sbatch --tmp= when set (MB on most schedulers)
STAGE_LOCAL=1
TREE=""
PARTITION=""
ACCOUNT=""
DO_CONCAT=1
DRY_RUN=0
WAIT=1                # block driver until SLURM finishes; --no-wait to detach
TAFFY="${TAFFY:-$(command -v taffy || true)}"

usage() {
    cat >&2 <<EOF
gerp_shard_slurm.sh -- run \`taffy gerp\` on a universal MAF in N SLURM shards
                       (column-range sharding; default stages input to \$TMPDIR)

Required:
  -i FILE       Input universal MAF (.maf.gz with a sibling .tui index)
  -o DIR        Output directory (created if missing)

Optional:
  -n INT        Number of shards (default 10; with local-stage prefer
                fewer-larger shards, e.g. -n 8 -T 32, to amortise the
                stage-in cost per task)
  -T INT        Threads per shard (default 10; --cpus-per-task=INT)
  --time HRS    Per-task wall limit in hours (default 24)
  --mem GB      Per-task memory in GB (default 32)
  --tmp GB      Per-task local-scratch requirement (--tmp=N to sbatch).
                Default unset.  Required when local-stage is on AND
                your cluster enforces \`--tmp\`; size to (1.5 * input
                bytes + output budget) so the stage + outputs fit.
  --no-stage-local
                Disable copying INPUT + .tui to \$TMPDIR; read from the
                original path.  Saves the stage time but every per-
                column .tui query hits the network FS.  Only sensible
                on local-filesystem inputs.
  --tree FILE   Newick tree override (default: \`# hal\` from MAF header).
                Also staged to \$TMPDIR when local-stage is on.
  --partition X SLURM partition (--partition=X)
  --account X   SLURM account (--account=X)
  --no-concat   Skip the post-array concatenation job
  --no-wait     Submit and detach (default: driver blocks until SLURM
                returns).  Without --no-wait the driver passes -W to
                sbatch so it exits only after the array (and concat,
                if any) completes -- handy when wrapping this script
                in a Makefile or a pipeline.  stderr from each task
                lands in a separate .err.log next to the .log file.
  --dry-run     Print the sbatch commands but do not submit
  -h            Help
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)              INPUT="$2"; shift 2;;
        -o)              OUTDIR="$2"; shift 2;;
        -n)              N_SHARDS="$2"; shift 2;;
        -T)              THREADS="$2"; shift 2;;
        --time)          TIME_HOURS="$2"; shift 2;;
        --mem)           MEM_GB="$2"; shift 2;;
        --tmp)           TMP_GB="$2"; shift 2;;
        --no-stage-local) STAGE_LOCAL=0; shift;;
        --tree)          TREE="$2"; shift 2;;
        --partition)     PARTITION="$2"; shift 2;;
        --account)       ACCOUNT="$2"; shift 2;;
        --no-concat)     DO_CONCAT=0; shift;;
        --no-wait)       WAIT=0; shift;;
        --dry-run)       DRY_RUN=1; shift;;
        -h|--help)       usage 0;;
        *)               echo "unknown arg: $1" >&2; usage 1;;
    esac
done

[[ -n "$INPUT"  ]] || { echo "ERROR: -i required" >&2; usage 1; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: -o required" >&2; usage 1; }
[[ -n "$TAFFY"  ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -f "$INPUT"      ]] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
[[ -f "${INPUT}.tui" ]] || { echo "ERROR: $INPUT has no .tui sibling -- run \`taffy index -u\` first" >&2; exit 1; }
[[ -z "$TREE" || -f "$TREE" ]] || { echo "ERROR: --tree file not found: $TREE" >&2; exit 1; }

mkdir -p "$OUTDIR"
echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> input:         $INPUT (.tui present)"
echo ">> shards:        $N_SHARDS,  threads/shard: $THREADS"
echo ">> local-stage:   $([[ $STAGE_LOCAL -eq 1 ]] && echo "ON (copies to \$TMPDIR)" || echo "OFF (reads from network)")"
[[ -n "$TMP_GB" ]] && echo ">> --tmp request: ${TMP_GB} GB per task"

# Belt-and-suspenders sizing hint when local-stage is on.
if [[ "$STAGE_LOCAL" -eq 1 ]]; then
    INPUT_BYTES=$(stat -c %s "$INPUT" 2>/dev/null || stat -f %z "$INPUT")
    TUI_BYTES=$(stat -c %s "$INPUT.tui" 2>/dev/null || stat -f %z "$INPUT.tui")
    STAGE_GB=$(( (INPUT_BYTES + TUI_BYTES) / (1024**3) ))
    echo ">> stage-in size: ~${STAGE_GB} GB per task (input + .tui)"
    if [[ -z "$TMP_GB" ]]; then
        echo ">> hint:          if your cluster advertises --tmp, consider --tmp $(( STAGE_GB + 100 )) (input + small output budget); otherwise omit"
    fi
fi

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
n_emit, n_skip = 0, 0
with open(dst, "w") as out:
    out.write("shard\tcol_lo\tcol_hi\n")
    for k in range(N):
        lo = (k     * T) // N
        hi = ((k + 1) * T) // N
        # T < N rounds several shards to lo == hi.  Skip them in the manifest
        # so the runner doesn't fire an empty `taffy gerp --columnRange`.
        # SLURM will still run the array slot but the runner's manifest
        # lookup returns no row and the runner exits 0 cleanly.
        if hi <= lo:
            n_skip += 1
            continue
        out.write(f"{k}\t{lo}\t{hi}\n")
        n_emit += 1
print(f"shard size: ~{T // N:,} columns each (last shard +{T - (T // N) * N} for remainder)",
      file=sys.stderr)
if n_skip > 0:
    print(f"skipped {n_skip} empty shards (T < N_SHARDS); manifest has {n_emit} rows",
          file=sys.stderr)
PY

# --- Step 3: write the per-shard runner. ------------------------------
# When STAGE_LOCAL is on, each task copies INPUT + .tui (+ optional tree)
# to $TMPDIR up front then runs gerp locally.  $TMPDIR is set by SLURM's
# job prolog on every reasonable cluster; we fall back to /tmp if not.
RUNNER="$OUTDIR/run_shard.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -euo pipefail
K=\${SLURM_ARRAY_TASK_ID:-0}
INPUT="$INPUT"
OUTDIR="$OUTDIR"
THREADS="$THREADS"
TAFFY="$TAFFY"
TREE="$TREE"
STAGE_LOCAL=$STAGE_LOCAL
MANIFEST="$OUTDIR/manifest.tsv"

# Lookup this shard's column range from the manifest (skip header).
read SHARD_ID COL_LO COL_HI < <(awk -v k="\$K" 'NR>1 && \$1==k {print \$1, \$2, \$3; exit}' "\$MANIFEST")
if [[ -z "\${COL_LO:-}" ]]; then
    # T < N_SHARDS pruned this shard out of the manifest; nothing to do.
    echo "shard \${K}: empty (T < N_SHARDS), skipping"
    exit 0
fi

RS_OUT="\$OUTDIR/shard_\${K}.rs.wig.gz"
DEPTH_OUT="\$OUTDIR/shard_\${K}.depth.wig.gz"

# Idempotency: if both outputs exist + non-empty, skip without staging.
if [[ -s "\$RS_OUT" && -s "\$DEPTH_OUT" ]]; then
    echo "shard \${K}: outputs present, skipping"
    exit 0
fi

# Resolve the scratch dir.  SLURM sets \$TMPDIR per task; if not set
# (running outside SLURM, or unusual cluster), make our own under /tmp.
SCRATCH="\${TMPDIR:-/tmp/taffy_gerp_\${SLURM_JOB_ID:-\$\$}_\${K}}"
mkdir -p "\$SCRATCH"
# Per-task subdir so colocated array shards on the same node (shared
# node-local \$TMPDIR, common on clusters where TMPDIR is /data/tmp
# without a per-task suffix) don't race on the stage path.  Job id +
# array task id is unique across the whole cluster lifetime.
STAGE_DIR="\$SCRATCH/taffy_stage_\${SLURM_JOB_ID:-\$\$}_\${K}"
# Clean up our scratch on exit -- SLURM's prolog usually does this too,
# but be explicit so failure modes don't leave 1.5 TB lying around.
# IMPORTANT: only rm OUR per-task subdir, not the entire SCRATCH (which
# is shared with sibling tasks on the same node in the colocated case).
trap 'rm -rf "\$STAGE_DIR" 2>/dev/null || true' EXIT

if [[ "\$STAGE_LOCAL" -eq 1 ]]; then
    mkdir -p "\$STAGE_DIR"
    BASENAME=\$(basename "\$INPUT")
    LOCAL_INPUT="\$STAGE_DIR/\$BASENAME"
    echo "shard \${K}: staging INPUT + .tui to \$STAGE_DIR ..."
    t0=\$SECONDS
    cp "\$INPUT"     "\$LOCAL_INPUT"
    cp "\$INPUT.tui" "\$LOCAL_INPUT.tui"
    if [[ -n "\$TREE" ]]; then
        LOCAL_TREE="\$STAGE_DIR/\$(basename "\$TREE")"
        cp "\$TREE" "\$LOCAL_TREE"
    fi
    echo "shard \${K}: stage-in took \$((SECONDS - t0)) s"
    GERP_INPUT="\$LOCAL_INPUT"
    GERP_TREE="\${LOCAL_TREE:-}"
else
    GERP_INPUT="\$INPUT"
    GERP_TREE="\$TREE"
fi
TREE_FLAG=\$([[ -n "\$GERP_TREE" ]] && echo "-t \$GERP_TREE" || echo "")

# Outputs land in scratch first; cp -> mv (rename within OUTDIR) at the
# end so partial writes never make it to a path that the idempotency
# check would mistake for completed work.
RS_TMP="\$SCRATCH/shard_\${K}.rs.wig.gz"
DEPTH_TMP="\$SCRATCH/shard_\${K}.depth.wig.gz"

echo "shard \${K}: running gerp on columns [\${COL_LO}, \${COL_HI}) ..."
t0=\$SECONDS
# -l INFO so the slurm log shows phase-1 progress + paralog counts +
# blocks/cols summary as it runs (without it the shard is silent for hours).
"\$TAFFY" gerp -i "\$GERP_INPUT" --columnRange "\${COL_LO}-\${COL_HI}" -c -T "\$THREADS" \\
    -l INFO \\
    \$TREE_FLAG \\
    -o "\$RS_TMP" -D "\$DEPTH_TMP"
echo "shard \${K}: gerp took \$((SECONDS - t0)) s"

echo "shard \${K}: staging out to \$OUTDIR ..."
t0=\$SECONDS
cp "\$RS_TMP"    "\$RS_OUT".tmp    && mv "\$RS_OUT".tmp    "\$RS_OUT"
cp "\$DEPTH_TMP" "\$DEPTH_OUT".tmp && mv "\$DEPTH_OUT".tmp "\$DEPTH_OUT"
echo "shard \${K}: stage-out took \$((SECONDS - t0)) s"

# Drop local copies eagerly even though trap will sweep on exit -- frees
# scratch for the next task to land on this node.
rm -f "\$RS_TMP" "\$DEPTH_TMP"
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
    --error="$OUTDIR/slurm_%A_%a.err.log"
    -J taffy_gerp_shard
)
[[ -n "$TMP_GB"    ]] && SBATCH_ARGS+=( --tmp="${TMP_GB}G" )
[[ -n "$PARTITION" ]] && SBATCH_ARGS+=( --partition="$PARTITION" )
[[ -n "$ACCOUNT"   ]] && SBATCH_ARGS+=( --account="$ACCOUNT" )

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo ">> DRY RUN -- would submit:"
    echo "sbatch ${SBATCH_ARGS[*]} --parsable $RUNNER"
    ARRAY_JOB=DRY
else
    echo ">> submitting array..."
    ARRAY_JOB=$(sbatch "${SBATCH_ARGS[@]}" --parsable "$RUNNER")
    echo ">> array job id:  $ARRAY_JOB"
    echo ">> per-task logs: $OUTDIR/slurm_${ARRAY_JOB}_<task>.log (stdout)"
    echo ">>                $OUTDIR/slurm_${ARRAY_JOB}_<task>.err.log (stderr)"
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
        --error="$OUTDIR/slurm_concat_%j.err.log"
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

# --- Step 7: optionally block until everything finishes. --------------
# `sbatch --wait <job>` blocks the caller until <job> finishes, with the
# job's exit code as the caller's exit code.  Waiting on the concat
# (which has afterok-dep on the array) covers the array too -- it can't
# run until the array succeeds, and won't be reported done until itself
# finishes.  If --no-concat, wait on the array instead.
if [[ "$DRY_RUN" -ne 1 && "$WAIT" -eq 1 ]]; then
    if [[ "$DO_CONCAT" -eq 1 ]]; then
        WAIT_JOB="$CONCAT_JOB"
        WAIT_KIND="concat (which depends on array $ARRAY_JOB)"
    else
        WAIT_JOB="$ARRAY_JOB"
        WAIT_KIND="array"
    fi
    echo ">> waiting for $WAIT_KIND job $WAIT_JOB ..."
    # --wait without resubmitting: re-attach with `sattach`-style, but
    # the cleanest way is a fresh sbatch that depends on the target.
    # Cheaper: just poll squeue.  Polling avoids spawning an extra job
    # and respects the partition / account constraints already in place.
    while squeue -j "$WAIT_JOB" -h -o "%T" 2>/dev/null | grep -qE "PENDING|RUNNING|CONFIGURING|COMPLETING|RESIZING|SUSPENDED|REQUEUED"; do
        sleep 60
    done
    # Read the final state via sacct (squeue drops completed jobs after a
    # short retention window, so squeue alone can race the exit).
    FINAL_STATE=$(sacct -j "$WAIT_JOB" -X -n -o State 2>/dev/null | head -1 | tr -d ' ')
    echo ">> $WAIT_KIND final state: ${FINAL_STATE:-UNKNOWN}"
    case "$FINAL_STATE" in
        COMPLETED) ;;
        *)         echo ">> NON-SUCCESS state -- check $OUTDIR/slurm_*.err.log"; exit 1;;
    esac
fi

echo ">> done."
