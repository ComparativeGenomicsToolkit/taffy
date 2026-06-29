#!/bin/bash
#
# Per-reference bigMafSummary for a BIG universal alignment (e.g. 577-way) -- SLURM
# =================================================================================
# Produces ONE bed3+4 bigMafSummary per leaf genome (`taffy summary --allRefs`)
# from a universal TAF/MAF, via column-range sharding.  Two dependent SLURM arrays:
#
#   A) SCAN  (array of N shards): each
#        taffy summary --allRefs --shard i/N -i UNI -o SCRATCH
#      .tui-seeks its disjoint column slice [i*T/N,(i+1)*T/N), scores every leaf
#      present, and emits that leaf's RAW single-covered records as text to
#      SCRATCH/<ref>.shard<i>.recs  (no sort/merge in the shard).
#
#   B) MERGE (array of M jobs, afterok A): each
#        taffy summary --shardMerge --refShard j/M -o SCRATCH
#      handles 1/M of the references (sorted + interleaved), gathering each ref's
#      <ref>.shard*.recs across ALL N shards, bucket-by-chrom + sort+merge ->
#      SCRATCH/<ref>.bed.  The per-reference sort is global, so each <ref>.bed
#      equals a non-sharded `taffy summary -r <ref>` run.
#
#   C) FINALIZE (afterok B): mv SCRATCH/*.bed -> OUTDIR, rm the *.recs temps.
#
# No input prep: reuses the existing UNI.tui (the seek index).  bigBed conversion
# (bedToBigBed per ref, needs each ref's chrom.sizes) is a SEPARATE step.
#
# SIZING (577-way ballpark -- read this before launching):
#   * Total compute ~4000-8000 core-hours (D^2 scoring; single-threaded scan, so
#     wall ~= total / concurrent slots; ~6-11 h at ~700 concurrent cores).
#   * Per shard task is SINGLE-THREADED (the scan); -T only threads the MERGE.
#   * SCRATCH must be HUGE: the raw .recs temps are ~tens of TB (the full D^2
#     record set hits disk before MERGE crunches it to ~1-2 TB of beds).  They are
#     removed only after MERGE succeeds.  Put SCRATCH on a roomy parallel FS.
#   * Temp file count ~ N * (refs-per-shard); with N=725 that's ~100-400 k files.
#     Smaller -s => more shards => shorter tasks + finer load-balance but MORE
#     files; larger -s => fewer files but longer, lumpier tasks.  Block density
#     varies along the column axis, so equal-width shards are uneven -- keep N well
#     above your concurrent-slot count so the scheduler absorbs the stragglers.
#   * MERGE --tmp is node-local (one ref's packed buckets at a time, ~tens of GB
#     for a big ref like the human/backbone); keep it off the shared FS.
#
# Usage: uni_summary_shard_slurm.sh -i UNI.taf.gz -o OUTDIR [options]
#   -i FILE        universal TAF/MAF (its <FILE>.tui MUST already exist)
#   -o DIR         output directory for the final per-reference .bed files
#   --scratch DIR  working dir for the .recs temps + in-progress beds (default:
#                  OUTDIR; needs ~tens of TB -- point at a roomy parallel FS)
#   -s INT         columns per shard (default 1e8 -> ~725 shards on a 577-way);
#                  sets N = ceil(T / s)
#   -M INT         number of parallel MERGE jobs (default 64)
#   -T INT         merge threads per job (default 8; per-chrom sort+merge workers)
#   --localtmp DIR node-local scratch for MERGE per-chrom buckets (default
#                  $TMPDIR or /tmp); must hold ~the biggest single ref's buckets
#   --scan-mem GB  --scan-time HRS   sbatch --mem/--time per SCAN shard (default 16G / 24h)
#   --merge-mem GB --merge-time HRS  sbatch --mem/--time per MERGE job (default 32G / 12h)
#   --partition P  --account A       passed through to every sbatch
#   --dry-run      write the job scripts + print the sbatch lines, don't submit
#   --wait         submit, then BLOCK until FINALIZE leaves the queue
# Env: TAFFY overrides the binary.
#
# Requirements: taffy (summary --allRefs/--shard/--shardMerge/--refShard), sbatch.

set -euo pipefail

INPUT=""; OUTDIR=""; SCRATCH=""; SHARD=100000000; NMERGE=64; THREADS=8
LOCALTMP=""; SCAN_MEM=16; SCAN_TIME=24; MERGE_MEM=32; MERGE_TIME=12
PARTITION=""; ACCOUNT=""; DRYRUN=0; WAIT=0
TAFFY="${TAFFY:-taffy}"

usage() { sed -n '46,71p' "$0" >&2; exit "${1:-1}"; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT=$2; shift 2;;
    -o) OUTDIR=$2; shift 2;;
    --scratch) SCRATCH=$2; shift 2;;
    -s) SHARD=$2; shift 2;;
    -M) NMERGE=$2; shift 2;;
    -T) THREADS=$2; shift 2;;
    --localtmp) LOCALTMP=$2; shift 2;;
    --scan-mem) SCAN_MEM=$2; shift 2;;
    --scan-time) SCAN_TIME=$2; shift 2;;
    --merge-mem) MERGE_MEM=$2; shift 2;;
    --merge-time) MERGE_TIME=$2; shift 2;;
    --partition) PARTITION=$2; shift 2;;
    --account) ACCOUNT=$2; shift 2;;
    --dry-run) DRYRUN=1; shift;;
    --wait) WAIT=1; shift;;
    -h|--help) usage 0;;
    *) echo "ERROR: unknown arg: $1" >&2; usage;;
  esac
done

[ -n "$INPUT" ] && [ -n "$OUTDIR" ] || usage
[ -e "$INPUT" ]     || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
[ -e "$INPUT.tui" ] || { echo "ERROR: $INPUT.tui missing -- the shard seek needs it (build with taffy index -u)" >&2; exit 1; }
command -v "$TAFFY" >/dev/null 2>&1 || { echo "ERROR: taffy not found (set \$TAFFY)" >&2; exit 1; }
(( SHARD > 0 ))  || { echo "ERROR: -s ($SHARD) must be > 0" >&2; exit 1; }
(( NMERGE > 0 )) || { echo "ERROR: -M ($NMERGE) must be > 0" >&2; exit 1; }
[ -n "$SCRATCH" ] || SCRATCH="$OUTDIR"
[ -n "$LOCALTMP" ] || LOCALTMP='${TMPDIR:-/tmp}'   # expanded inside the job, on the node
mkdir -p "$OUTDIR" "$SCRATCH"

T=$("$TAFFY" stats -i "$INPUT" -u)
[[ "$T" =~ ^[0-9]+$ ]] || { echo "ERROR: 'taffy stats -i $INPUT -u' did not return an integer T (got: $T)" >&2; exit 1; }
NSHARD=$(( (T + SHARD - 1) / SHARD ))
(( NSHARD >= 1 )) || NSHARD=1
(( NMERGE > NSHARD*0 )) || true
RUNNER="$SCRATCH/.uni_summary.scan.sh"
MERGER="$SCRATCH/.uni_summary.merge.sh"
FINALIZE="$SCRATCH/.uni_summary.finalize.sh"
echo ">> T=$T columns  shard=$SHARD cols  scan_shards=$NSHARD  merge_jobs=$NMERGE  -> $OUTDIR" >&2

# ---- A) SCAN: one column-shard per array task -> per-ref raw .recs ----
cat > "$RUNNER" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=${SCAN_MEM}G
#SBATCH --time=${SCAN_TIME}:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
"$TAFFY" summary --allRefs --shard \$SLURM_ARRAY_TASK_ID/$NSHARD \\
  -i "$INPUT" -o "$SCRATCH"
EOF

# ---- B) MERGE: 1/M of the references per array task (afterok A) ----
cat > "$MERGER" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=${MERGE_MEM}G
#SBATCH --time=${MERGE_TIME}:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
"$TAFFY" summary --shardMerge --refShard \$SLURM_ARRAY_TASK_ID/$NMERGE \\
  -o "$SCRATCH" --threads \$SLURM_CPUS_PER_TASK --tmp "$LOCALTMP"
EOF

# ---- C) FINALIZE: collect beds to OUTDIR, drop the .recs temps (afterok B) ----
cat > "$FINALIZE" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=4:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
nbed=\$(ls "$SCRATCH"/*.bed 2>/dev/null | wc -l)
(( nbed > 0 )) || { echo "FINALIZE: no .bed files in $SCRATCH -- merge failed" >&2; exit 1; }
if [ "$(readlink -f "$SCRATCH")" != "$(readlink -f "$OUTDIR")" ]; then
  mv "$SCRATCH"/*.bed "$OUTDIR"/
fi
rm -f "$SCRATCH"/*.shard*.recs
echo "FINALIZE: \$nbed per-reference beds in $OUTDIR; .recs temps removed"
EOF

if (( DRYRUN )); then
  echo "[dry-run] wrote $RUNNER, $MERGER, $FINALIZE" >&2
  echo "[dry-run] sbatch --array=0-$((NSHARD-1)) $RUNNER" >&2
  echo "[dry-run] sbatch --dependency=afterok:<A> --array=0-$((NMERGE-1)) $MERGER" >&2
  echo "[dry-run] sbatch --dependency=afterok:<B> $FINALIZE" >&2
  exit 0
fi
A_ID=$(sbatch --parsable --array=0-$((NSHARD-1)) "$RUNNER")
echo ">> SCAN array submitted: job $A_ID ($NSHARD shards)" >&2
B_ID=$(sbatch --parsable --dependency=afterok:"$A_ID" --array=0-$((NMERGE-1)) "$MERGER")
echo ">> MERGE array submitted (afterok:$A_ID): job $B_ID ($NMERGE jobs)" >&2
F_ID=$(sbatch --parsable --dependency=afterok:"$B_ID" "$FINALIZE")
echo ">> FINALIZE submitted (afterok:$B_ID): job $F_ID -> $OUTDIR" >&2

if (( WAIT )); then
  # Real SLURM jobs (two arrays + a finalize); the build still completes if THIS
  # poll dies -- run under nohup/tmux for a multi-hour build.  afterok means any
  # shard/merge failure cancels downstream, so FINALIZE leaves the queue without
  # writing -> the bed-count check reports the failure.
  echo ">> [--wait] blocking until job $F_ID completes (polling every 60s)..." >&2
  sleep 5
  while [[ -n "$(squeue -h -j "$F_ID" 2>/dev/null)" ]]; do sleep 60; done
  nbed=$(ls "$OUTDIR"/*.bed 2>/dev/null | wc -l)
  if (( nbed > 0 )); then
    echo ">> done -> $nbed per-reference beds in $OUTDIR" >&2
  else
    echo ">> BUILD FAILED: no beds in $OUTDIR -- inspect scan ($A_ID), merge ($B_ID), finalize ($F_ID)" >&2
    exit 1
  fi
fi
