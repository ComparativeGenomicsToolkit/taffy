#!/bin/bash
#
# Per-reference bigMafSummary for a BIG universal alignment (e.g. 577-way) -- SLURM
# =================================================================================
# Produces ONE bed3+4 bigMafSummary bigBed per leaf genome from a universal TAF/MAF
# using PER-REFERENCE sharding: NO shared scratch for the intermediate records.
#
# Why per-reference (not column-range): a universal alignment is ordered on
# universal columns, but a leaf reference is SHUFFLED through that order, so a
# reference's records are scattered across the whole file.  Column-sharding (the
# older uni_summary_shard_slurm.sh) therefore had to spill the entire D^2 record
# set (~tens of TB) to a SHARED FS before it could merge per reference.  Instead,
# each job here scans the whole file once, masters only its 1/M slice of the
# references (`taffy summary --allRefs --refSubset j/M`), and merges them to
# NODE-LOCAL disk -- only the finished ~2 GB bed (or its bigBed) ever touches the
# shared FS.  Cost: the file is read M times (I/O traded for scratch).  Bonus: a
# single-reference master scores ref x D, not D^2, so there is no dense-root
# straggler -- the per-block work is even.
#
# Two dependent SLURM stages:
#   A) SUMMARIZE (array of M jobs): job j runs
#        taffy summary --allRefs --refSubset j/M -i UNI -o LOCAL/beds --tmp LOCAL/tmp
#      mastering the ~refs/M leaves with name-hash %M==j, writing each <ref>.bed +
#      <ref>.chrom.sizes (from UNI.tui) to NODE-LOCAL disk, then (unless
#      --no-bigbed) bedToBigBed each straight to OUTDIR/<ref>.summary.bb.  The raw
#      per-reference records live and die under LOCAL/tmp (node-local), bounded by
#      --mergeMem; nothing intermediate hits the shared FS.
#   C) FINALIZE (afterok A): count the OUTDIR/*.summary.bb (or *.bed) and report.
#
# No input prep beyond the existing UNI.tui (chrom sizes for bedToBigBed).
#
# SIZING (577-way ballpark -- read before launching):
#   * Total compute ~= the same Sum(D^2) scoring as one fan-out scan, but spread
#     over M passes; wall ~= total / concurrent slots.  Single-ref masters mean NO
#     dense-root straggler.
#   * I/O: the file is decompressed+scanned M times (~M * file-size of reads).
#     Smaller M = fewer reads but each job masters MORE refs => MORE node-local
#     disk; larger M = less node-local but more reads.  Tune M to your node-local.
#   * NODE-LOCAL --tmp must hold the job's masters' raw records AT ONCE during the
#     scan: ~ (refs/M) references' packed records, dominated by the deepest ref
#     (~70 GB for a 577-way backbone genome).  With refs=577 and M=64 that is ~9
#     refs => budget a few hundred GB of node-local; raise M if that is too much.
#   * Peak RAM per job ~= --mergeMem (merge budget) + the biggest single bigBed
#     build (sequential, so max not sum); --mem covers both.
#   * Shared FS holds only the final per-reference bigBeds (~1-2 TB total).
#
# Usage: uni_summary_refsubset_slurm.sh -i UNI.taf.gz -o OUTDIR [options]
#   -i FILE        universal TAF/MAF (its <FILE>.tui MUST already exist)
#   -o DIR         output directory for the per-reference .summary.bb (or .bed)
#   -M INT         number of passes/jobs; each masters refs/M references (default 64)
#   -T INT         threads per job (default 8; --cpus-per-task)
#   --mergeMem GB  per-job in-RAM chrom-merge budget passed to taffy (default 8)
#   --localtmp DIR node-local scratch root (default $TMPDIR or /tmp); per-task subdir
#                  holds one pass's raw records + beds -- keep OFF the shared FS
#   --no-bigbed    stop at .bed (write beds to OUTDIR; skip bedToBigBed)
#   --bedToBigBed PATH  the bedToBigBed binary (default: bedToBigBed on PATH)
#   --as FILE      bed AS schema for bedToBigBed (default: a built-in mafSummary.as)
#   --mem GB       --time HRS   sbatch --mem/--time per SUMMARIZE job (default 24G / 24h)
#   --partition P  --account A  passed through to every sbatch
#   --dry-run      write the job scripts + print the sbatch lines, don't submit
#   --wait         submit, then BLOCK until FINALIZE leaves the queue
# Env: TAFFY overrides the binary.
#
# Requirements: taffy (summary --allRefs --refSubset), sbatch, and (unless
# --no-bigbed) bedToBigBed.

set -euo pipefail

INPUT=""; OUTDIR=""; NPASS=64; THREADS=8; MERGEMEM=8; LOCALTMP=""
MEM=24; TIME=24; BIGBED=1; BTB="bedToBigBed"; ASFILE=""
PARTITION=""; ACCOUNT=""; DRYRUN=0; WAIT=0
TAFFY="${TAFFY:-taffy}"

usage() { sed -n '52,72p' "$0" >&2; exit "${1:-1}"; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT=$2; shift 2;;
    -o) OUTDIR=$2; shift 2;;
    -M) NPASS=$2; shift 2;;
    -T) THREADS=$2; shift 2;;
    --mergeMem) MERGEMEM=$2; shift 2;;
    --localtmp) LOCALTMP=$2; shift 2;;
    --no-bigbed) BIGBED=0; shift;;
    --bedToBigBed) BTB=$2; shift 2;;
    --as) ASFILE=$2; shift 2;;
    --mem) MEM=$2; shift 2;;
    --time) TIME=$2; shift 2;;
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
[ -e "$INPUT.tui" ] || { echo "ERROR: $INPUT.tui missing -- needed for chrom.sizes (build with taffy index -u)" >&2; exit 1; }
command -v "$TAFFY" >/dev/null 2>&1 || { echo "ERROR: taffy not found (set \$TAFFY)" >&2; exit 1; }
(( NPASS > 0 )) || { echo "ERROR: -M ($NPASS) must be > 0" >&2; exit 1; }
if (( BIGBED )); then
  command -v "$BTB" >/dev/null 2>&1 || [ -x "$BTB" ] || { echo "ERROR: bedToBigBed not found ('$BTB'); set --bedToBigBed PATH or use --no-bigbed" >&2; exit 1; }
  if [ -n "$ASFILE" ]; then [ -e "$ASFILE" ] || { echo "ERROR: --as file not found: $ASFILE" >&2; exit 1; }; fi
fi
[ -n "$LOCALTMP" ] || LOCALTMP='${TMPDIR:-/tmp}'   # expanded inside the job, on the node
mkdir -p "$OUTDIR"

# A built-in mafSummary.as so no kent checkout is needed on the cluster (overridable via --as).
# Lives in OUTDIR (shared) so every node sees the same path.
if (( BIGBED )) && [ -z "$ASFILE" ]; then
  ASFILE="$OUTDIR/.mafSummary.as"
  cat > "$ASFILE" <<'AS'
table mafSummary
"Positions and scores for alignment blocks"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string src;        "Sequence name or database of alignment"
    float  score;      "Floating point score."
    char[1] leftStatus;  "Gap/break annotation for preceding block"
    char[1] rightStatus; "Gap/break annotation for following block"
    )
AS
fi

SUMRUN="$OUTDIR/.uni_summary.summarize.sh"
FINALIZE="$OUTDIR/.uni_summary.finalize.sh"
echo ">> per-reference summary: M=$NPASS passes, -T$THREADS, --mergeMem ${MERGEMEM}G, bigbed=$BIGBED -> $OUTDIR" >&2

# ---- A) SUMMARIZE: one ref-subset per array task, node-local raw -> bigBed ----
cat > "$SUMRUN" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=${MEM}G
#SBATCH --time=${TIME}:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
export OMP_WAIT_POLICY=passive   # OMP scan threads sleep (not spin) between dense blocks
LOCAL="$LOCALTMP/uni_sum_\$SLURM_ARRAY_TASK_ID.\$\$"
mkdir -p "\$LOCAL/beds" "\$LOCAL/tmp"
trap 'rm -rf "\$LOCAL"' EXIT
# scan the whole file once, master only this pass's 1/M of the references, merge
# to node-local disk (raw records bounded by --mergeMem, never touch shared FS):
"$TAFFY" summary --allRefs --refSubset \$SLURM_ARRAY_TASK_ID/$NPASS \\
  -i "$INPUT" -o "\$LOCAL/beds" --tmp "\$LOCAL/tmp" \\
  --threads \$SLURM_CPUS_PER_TASK --mergeMem $MERGEMEM
shopt -s nullglob
made=0
for bed in "\$LOCAL"/beds/*.bed; do
  ref=\$(basename "\$bed" .bed)
  if [ "$BIGBED" = 1 ]; then
    sizes="\$LOCAL/beds/\$ref.chrom.sizes"
    [ -s "\$sizes" ] || { echo "SUMMARIZE: missing \$sizes (is $INPUT.tui complete?)" >&2; exit 1; }
    # the bed is chrom-contiguous + within-chrom start-sorted already (no resort)
    "$BTB" -type=bed3+4 -as="$ASFILE" -tab "\$bed" "\$sizes" "$OUTDIR/\$ref.summary.bb"
  else
    mv "\$bed" "$OUTDIR/\$ref.bed"
  fi
  made=\$((made+1))
done
echo "SUMMARIZE pass \$SLURM_ARRAY_TASK_ID/$NPASS: \$made references -> $OUTDIR"
(( made > 0 )) || { echo "SUMMARIZE: pass produced 0 references -- check the input/tree" >&2; exit 1; }
EOF

# ---- C) FINALIZE: count outputs, report ----
cat > "$FINALIZE" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
if [ "$BIGBED" = 1 ]; then
  n=\$(ls "$OUTDIR"/*.summary.bb 2>/dev/null | wc -l)
  (( n > 0 )) || { echo "FINALIZE: no .summary.bb in $OUTDIR -- the summarize array failed" >&2; exit 1; }
  rm -f "$OUTDIR/.mafSummary.as"
  echo "FINALIZE: \$n per-reference .summary.bb in $OUTDIR"
else
  n=\$(ls "$OUTDIR"/*.bed 2>/dev/null | wc -l)
  (( n > 0 )) || { echo "FINALIZE: no .bed in $OUTDIR -- the summarize array failed" >&2; exit 1; }
  echo "FINALIZE: \$n per-reference .bed in $OUTDIR"
fi
EOF

if (( DRYRUN )); then
  echo "[dry-run] wrote $SUMRUN, $FINALIZE" >&2
  echo "[dry-run] sbatch --array=0-$((NPASS-1)) $SUMRUN" >&2
  echo "[dry-run] sbatch --dependency=afterok:<A> $FINALIZE" >&2
  exit 0
fi
A_ID=$(sbatch --parsable --array=0-$((NPASS-1)) "$SUMRUN")
echo ">> SUMMARIZE array submitted: job $A_ID ($NPASS passes)" >&2
F_ID=$(sbatch --parsable --dependency=afterok:"$A_ID" "$FINALIZE")
echo ">> FINALIZE submitted (afterok:$A_ID): job $F_ID -> $OUTDIR" >&2

if (( WAIT )); then
  # Real SLURM jobs; the build still completes if THIS poll dies -- run under
  # nohup/tmux for a multi-hour build.  afterok means any pass failure cancels
  # FINALIZE, so it leaves the queue without writing -> the count check reports it.
  echo ">> [--wait] blocking until job $F_ID completes (polling every 60s)..." >&2
  sleep 5
  while [[ -n "$(squeue -h -j "$F_ID" 2>/dev/null)" ]]; do sleep 60; done
  if (( BIGBED )); then n=$(ls "$OUTDIR"/*.summary.bb 2>/dev/null | wc -l); else n=$(ls "$OUTDIR"/*.bed 2>/dev/null | wc -l); fi
  if (( n > 0 )); then
    echo ">> done -> $n per-reference $( ((BIGBED)) && echo bigBeds || echo beds) in $OUTDIR" >&2
  else
    echo ">> BUILD FAILED: no outputs in $OUTDIR -- inspect summarize ($A_ID), finalize ($F_ID)" >&2
    exit 1
  fi
fi
