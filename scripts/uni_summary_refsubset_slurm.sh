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
#   * I/O: the file is decompressed+scanned M times (~M * file-size of reads); stage
#     it node-local (default) and OFF a slow NFS.
#   * RAM IS THE LIMIT (the raw is NEVER spilled -- it's merged in RAM as it scans):
#     peak RAM per job ~= (refs/M) x the merged BED size (~2-3 GB/ref on a 577-way).
#     So smaller M = MORE RAM (more refs resident) + fewer passes; larger M = less
#     RAM + more passes/reads.  M IS THE RAM KNOB: --mergeMem only sets the compaction
#     CADENCE (how often each ref window-bins), NOT a hard cap -- a deep ref still
#     grows to ~2x its merged size regardless.  So MEASURE peak RSS with a real
#     `--allRefs --refSubset 0/M` pass (NOT -r) and size SLURM --mem to it; if it is
#     too big, RAISE M (fewer refs/pass) -- do not expect --mergeMem to bound it.
#   * NODE-LOCAL --tmp now holds only the STAGED input copy (~the .taf.gz size) + the
#     per-ref beds -- NOT the raw records (the ~90 TB the old disk design spilled is
#     gone).  A few hundred GB, not TBs.
#   * M also CAPS CORES: the scan threads across a block's masters, so a pass uses at
#     most ~refs/M cores -- pick M so refs/M >= -T.  Decode is single-threaded (htslib
#     bgzf_mt leaks), so on a huge file the decode wall is a floor -- factor it in.
#   * Shared FS holds only the final per-reference bigBeds (~1-2 TB total).
#
# Usage: uni_summary_refsubset_slurm.sh -i UNI.taf.gz -o OUTDIR [options]
#   -i FILE        universal TAF/MAF (its <FILE>.tui MUST already exist)
#   -o DIR         output directory for the per-reference .summary.bb (or .bed)
#   -M INT         number of passes/jobs; each masters refs/M references (default 64)
#   -T INT         threads per job (default 8; --cpus-per-task)
#   --mergeMem GB  per-job in-RAM chrom-merge budget passed to taffy (default 8)
#   --localtmp DIR node-local scratch root (default $TMPDIR or /tmp); holds the STAGED
#                  input copy + this pass's beds (raw is in RAM) -- keep OFF shared FS
#   --no-bigbed    stop at .bed (write beds to OUTDIR; skip bedToBigBed)
#   --no-stage     don't copy the input to node-local first (default: stage it, so
#                  the M concurrent passes don't all stream the file off shared FS)
#   --only RANGE   run only these array task(s), e.g. 0 or 0-3 -- a fragment probe
#                  (one pass; skips FINALIZE; read MaxRSS/Elapsed from sacct)
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
PARTITION=""; ACCOUNT=""; DRYRUN=0; WAIT=0; STAGE=1; ONLY=""
TAFFY="${TAFFY:-taffy}"

usage() { sed -n '54,74p' "$0" >&2; exit "${1:-1}"; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT=$2; shift 2;;
    -o) OUTDIR=$2; shift 2;;
    -M) NPASS=$2; shift 2;;
    -T) THREADS=$2; shift 2;;
    --mergeMem) MERGEMEM=$2; shift 2;;
    --localtmp) LOCALTMP=$2; shift 2;;
    --no-bigbed) BIGBED=0; shift;;
    --no-stage) STAGE=0; shift;;
    --only) ONLY=$2; shift 2;;
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
#SBATCH --output=$OUTDIR/summarize-%A_%a.log
${PARTITION:+#SBATCH --partition=$PARTITION}
${ACCOUNT:+#SBATCH --account=$ACCOUNT}
set -euo pipefail
export OMP_WAIT_POLICY=passive   # OMP scan threads sleep (not spin) between dense blocks
LOCAL="$LOCALTMP/uni_sum_\$SLURM_ARRAY_TASK_ID.\$\$"
mkdir -p "\$LOCAL/beds"
trap 'rm -rf "\$LOCAL"' EXIT
# STAGE the input to node-local disk (default; --no-stage to skip): the M passes
# otherwise all stream the whole file off the shared FS at once.  Copy the .taf.gz;
# SYMLINK the .tui (its read is just the genome roster for chrom.sizes, not 100 GB).
IN="$INPUT"
if [ "$STAGE" = 1 ]; then
  echo "staging $INPUT -> \$LOCAL/in.taf.gz ..."
  cp "$INPUT" "\$LOCAL/in.taf.gz"
  ln -sf "$INPUT.tui" "\$LOCAL/in.taf.gz.tui"
  IN="\$LOCAL/in.taf.gz"
fi
# scan the whole file once, master only this pass's 1/M of the references, and
# merge IN RAM (the raw records are NEVER spilled to disk -- only the beds written;
# RAM peaks at ~2x this pass's refs' merged size, capped by SLURM --mem, not disk):
"$TAFFY" summary --allRefs --refSubset \$SLURM_ARRAY_TASK_ID/$NPASS \\
  -i "\$IN" -o "\$LOCAL/beds" \\
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
# an empty subset (no leaf's name-hash landed in this pass) is a legitimate no-op,
# NOT a failure -- still drop a per-pass sentinel so FINALIZE can confirm all $NPASS
# passes actually ran (raising M makes empty passes EXPECTED, not an error).
mkdir -p "$OUTDIR/.passes"
touch "$OUTDIR/.passes/\$SLURM_ARRAY_TASK_ID.done"
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
# COMPLETENESS: every one of the $NPASS passes must have left its sentinel (an
# empty subset still does).  FINALIZE is afterany, so this runs even when a pass
# failed -- reporting the shortfall instead of silently shipping a partial browser.
# A bare n>0 check would pass a 568/577 result; this catches it.
done_n=\$(ls "$OUTDIR"/.passes/*.done 2>/dev/null | wc -l)
if (( done_n != $NPASS )); then
  miss=""
  for j in \$(seq 0 $((NPASS-1))); do [ -f "$OUTDIR/.passes/\$j.done" ] || miss="\${miss:+\$miss,}\$j"; done
  echo "FINALIZE: only \$done_n/$NPASS passes completed -- FAILED pass(es): \$miss" >&2
  echo "  output is INCOMPLETE; rerun just those: sbatch --array=\$miss $SUMRUN ; then re-run FINALIZE" >&2
  exit 1
fi
if [ "$BIGBED" = 1 ]; then n=\$(ls "$OUTDIR"/*.summary.bb 2>/dev/null | wc -l); rm -f "$OUTDIR/.mafSummary.as"; else n=\$(ls "$OUTDIR"/*.bed 2>/dev/null | wc -l); fi
rm -rf "$OUTDIR/.passes"
echo "FINALIZE: all $NPASS passes complete -> \$n per-reference outputs in $OUTDIR"
EOF

ARR="${ONLY:-0-$((NPASS-1))}"
if (( DRYRUN )); then
  echo "[dry-run] wrote $SUMRUN$([ -z "$ONLY" ] && echo ", $FINALIZE")" >&2
  echo "[dry-run] sbatch --array=$ARR $SUMRUN" >&2
  [ -z "$ONLY" ] && echo "[dry-run] sbatch --dependency=afterany:<A> $FINALIZE  # afterany so a pass failure is REPORTED" >&2
  exit 0
fi
A_ID=$(sbatch --parsable --array="$ARR" "$SUMRUN")
echo ">> SUMMARIZE array submitted: job $A_ID (passes $ARR)" >&2

# ONE-FRAGMENT probe (--only): completeness is meaningless for a subset, so no FINALIZE.
if [ -n "$ONLY" ]; then
  echo ">> --only $ONLY: fragment probe, FINALIZE skipped.  When it finishes:" >&2
  echo "     sacct -j $A_ID --format=JobID,State,MaxRSS,Elapsed   # MaxRSS -> --mem, Elapsed -> per-pass wall" >&2
  echo "     less $OUTDIR/summarize-${A_ID}_*.log ; ls $OUTDIR/*.summary.bb   # per-ref rows + outputs" >&2
  exit 0
fi

F_ID=$(sbatch --parsable --dependency=afterany:"$A_ID" "$FINALIZE")   # afterany: run even if a pass fails, so the shortfall is reported
echo ">> FINALIZE submitted (afterany:$A_ID): job $F_ID -> $OUTDIR" >&2

if (( WAIT )); then
  # Real SLURM jobs; the build still completes if THIS poll dies -- run under
  # nohup/tmux for a multi-hour build.  FINALIZE is afterany, so it runs even when a
  # pass fails: it removes OUTDIR/.passes/ only when all $NPASS passes completed,
  # else leaves it and exits non-zero with the failed-pass list.
  echo ">> [--wait] blocking until job $F_ID completes (polling every 60s)..." >&2
  sleep 5
  while [[ -n "$(squeue -h -j "$F_ID" 2>/dev/null)" ]]; do sleep 60; done
  if (( BIGBED )); then n=$(ls "$OUTDIR"/*.summary.bb 2>/dev/null | wc -l); else n=$(ls "$OUTDIR"/*.bed 2>/dev/null | wc -l); fi
  if [ -d "$OUTDIR/.passes" ]; then
    echo ">> BUILD INCOMPLETE: not all $NPASS passes finished ($n outputs so far) -- see FINALIZE ($F_ID) for the failed-pass list; inspect summarize ($A_ID)" >&2
    exit 1
  fi
  echo ">> done -> $n per-reference $( ((BIGBED)) && echo bigBeds || echo beds) in $OUTDIR" >&2
fi
