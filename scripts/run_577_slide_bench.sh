#!/bin/bash
#
# run_577_slide_bench.sh -- one wrapper to drive all SIX comparisons for
# the 577-way / hg38-source conference slide.  Shared config block; invokes
# one bench driver per comparison.  NEVER submits sbatch itself -- it calls
# each driver with --dry-run (default) or lets each driver submit (--submit).
#
# THE SIX COMPARISONS (all full-577-way, hg38 source/reference):
#   1. MAF extraction   : tui (view -U query --noAncestors) vs tai (view) vs hal (hal2maf)
#                         vs bigmaf (bigMafToMaf)   -> taffy_view_bench_slurm.sh  [all -> hg38 leaf MAF]
#   2. base liftover    : taffy lift on the BASE .tui vs halLiftover
#                                                    -> taffy_lift_bench_slurm.sh
#   3. blockViz query   : taffyBlockViz vs halBlockViz
#                                                    -> taffy_blockviz_bench_slurm.sh
#   4. chained liftover : taffy lift on the CHAINED .tui vs UCSC liftOver
#                                                    -> taffy_lift_bench_slurm.sh
#   5. view sample      : random genome-wide MAF extraction (the #1 tools)
#                                                    -> taffy_view_sample_bench_slurm.sh
#   6. depth summary    : taffy lift --bigwig (depth bigWig) vs bigMafSummary
#                                                    -> uni_depth_summary_bench_slurm.sh
#
# THE HAL CAP, applied to the existing (#1/#2) scripts
# ----------------------------------------------------
# #1 (view) builds its OWN internal log-decade ladder and has no --no-hal,
# so a two-pass split is awkward (it would duplicate the taffy/tai/bb
# cells).  Instead a small, targeted --halMaxSize knob was ADDED to
# taffy_view_bench_slurm.sh: the hal2maf cell is simply not launched for
# ladder sizes above the cap (no run-then-timeout), while taffy/tai/bb run
# the full ladder.  One clean pass.  (The blockViz driver #3 has the same
# native --halMaxSize knob.)
#
# #2 (base lift) is a size x species matrix taking an explicit -S ladder.
# The lift driver has its OWN native --halMaxSize knob (same pattern as the
# view / blockViz drivers): the halLiftover cell is simply NOT launched for
# any (species,size) where size > the cap, while the taffy/liftOver cells
# run the full ladder.  So #2 is ONE clean pass over the full LIFT_SIZES
# ladder with --halMaxSize -- no two-pass split (the old split re-copied
# the 98 GB .tui for the second pass for nothing).
#
# #4 (chained lift, vs UCSC liftOver) has NO hal tool (liftOver is
# chain-based), so it needs no cap -- it runs the full ladder in one pass
# with --no-hal.
#
# THE CHAINED-.tui RESOLUTION, for #4   (symlink + a tiny driver option)
# ----------------------------------------------------------------------
# taffy_lift_bench_slurm.sh resolves its .tui as the sibling of the
# -u/--uniTaf path: it validates BOTH `-f "$UNI"` and `-f "$UNI.tui"`, and
# at runtime it stubs `$UNI` to a 0-byte file + `taffy lift` reads only
# `$UNI.tui`.  The chained .tui (..._g10000.tui) has no .uni.taf.gz sibling.
# We satisfy the .tui resolution with a SYMLINK:
#     : > $WORK/chained.uni.taf.gz                    # 0-byte stub
#     ln -sf <CHAINED_TUI>  $WORK/chained.uni.taf.gz.tui
#     --uniTaf $WORK/chained.uni.taf.gz
#
# BUT the symlink alone is NOT sufficient: before launching cells the lift
# driver enumerates REF chroms with `taffy stats -s -i <source>`, and
# `taffy stats -s` reads the SOURCE TAF/MAF HEADER (verified), not the .tui
# -- so a 0-byte stub crashes it (taf.c check_input_format assertion).  We
# therefore ALSO added a tiny `--statsSrc FILE` option to the lift driver:
# a readable TAF/MAF used ONLY for that chrom enumeration, decoupled from
# the -u/--uniTaf .tui that drives the actual lift.  The wrapper passes
# `--statsSrc <base uni source>` for #4 (the base uni shares hg38's chrom
# set with the chained .tui).  So #4 = symlink (for the .tui) + --statsSrc
# (for chrom enumeration).  Chosen over a header-bearing fake stub because
# that would require regenerating the exact chrom header from the source
# anyway; --statsSrc is the minimal honest decoupling.  REPORTED in summary.
#
# Usage:
#   run_577_slide_bench.sh [--test] [--submit] [config overrides...]
# Default is DRY-RUN: prints the sbatch command each driver would submit.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ======================================================================
# SHARED CONFIG BLOCK -- edit these (or override via CLI flags) for the
# real cluster run.  Everything big is a cluster-side path supplied here.
# ======================================================================
# -- Inputs (all cluster-side for the real run; locals shown as hints) --
BASE_TUI_MAF=""        # base universal MAF whose .tui sibling drives #1 view + #2 base-lift
                       #   (e.g. /path/vgp-577way-v1.uni.maf.gz, with .tui sibling)
BASE_TUI_TAF=""        # OR the TAF-anchored universal (.uni.taf.gz + .tui); set whichever you have
CHAINED_TUI=""         # chained .tui for #3 blockViz + #4 chained-lift
                       #   (e.g. /path/vgp-577way-v1.uni.taf.gz.chained_g10000.tui)
HAL=""                 # full 577-way .hal (hal2maf / halLiftover / blockVizBed)
HG38_MAF=""            # hg38-anchored MAF (.maf.gz + .tai) for #1 view's tai cells
BIGBED=""              # bigmaf bigBed for #1 view's bb cells
DEPTH_BW=""            # universal-depth bigWig for #6 (uni_depth_bigwig_slurm.sh output)
SUMMARY_BB=""          # hg38 bigMafSummary bigBed for #6 (bigBedToBed baseline)
CHAINS_DIR=""          # dir of UCSC chains named <REF>_vs_<genome_id>.chain.gz (#2 + #4)
TREE=""                # species tree .nwk (plots' divergence x-axis)
PANEL="$HERE/vgp577_hg38_panel.tsv"   # target-genome panel (-L)

# -- Reference / query geometry --
REF="GCA_000001405.15"          # hg38 genome id (matches the universal-MAF / chain prefix)
# Naming: hg38 is the accession GCA_000001405.15 in the universal MAF, the
# hg38-anchored MAF, the HAL and the chained tui -- EXCEPT the bigBed, whose
# chrom is the bare UCSC name (chr20).
VIEW_REF_CHROM="GCA_000001405.15.chr20"   # taffy view -r (genome.seq) for #1 tui+tai cells
HAL_GENOME="GCA_000001405.15"   # hal2maf --refGenome / blockViz tSpecies
HAL_SEQ="chr20"                 # hal2maf --refSequence / bare tChrom
BB_CHROM="chr20"                # bigMafToMaf chrom for #1 (the bigMaf = bare UCSC name)
BLOCKVIZ_CHROM="chr20"          # bare tChrom for #3 (both blockViz tools)
VIEW_START=1000000              # #1 view window start (skip the chr20 telomere prefix)
VIEW_NORM=0                     # #1 _norm cells (view|taffy norm): 0=OFF this run, set 1 to re-enable

# -- Ladders / panel knobs --
VIEW_MAX_SIZE=""                # #1 view --maxSize cap (empty = chrom end from .tai)
BLOCKVIZ_SIZES="1000,100000,500000,1000000,10000000,64000000"   # #3 ladder: 1k 100k 500k 1M 10M whole-chr20
LIFT_SIZES="1000,100000,1000000,10000000"   # #4 (chained tui --fast vs liftOver) ladder (1k 100k 1M 10M)
LIFT_PERCOL_SIZES="1000,100000,1000000"     # #2 (full tui per-column: hal vs tui) base-level ladder (1k 100k 1M)
LIFT_N_INTERVALS=100            # #2 + #4 random intervals per (species,size) cell
HAL_MAX_SIZE=10000000           # the cap: hal tools skipped above this (bp) in #1/#3/#5.  Raised
                                # 500k -> 10M so hal2maf reaches 1M/10M (a few more data points);
                                # HAL_TIME_BUDGET governs whether the big ones complete vs time out.
# -- #5 random-sample view bench --
VIEW_SAMPLE_SIZES="100,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000"  # #5 sizes (bp): linear 100->1M, 11 pts
VIEW_SAMPLE_N=10                # #5 random genome-wide regions per size
VIEW_SAMPLE_SEED=20260620       # #5 RNG seed (reproducible; same regions for every tool)
# -- #6 depth-summary bench (zoom-out: depth bigWig lift vs bigMafSummary) --
DEPTH_SUMMARY_SIZES="1000000,20000000,40000000,60000000,80000000,100000000,120000000,140000000,160000000,180000000,200000000"  # #6 sizes (bp): linear 1M->200M, 11 pts
DEPTH_SUMMARY_N=10              # #6 random hg38 regions per size
DEPTH_SUMMARY_SEED=20260620     # #6 RNG seed

# -- SLURM / runtime --
T_TOTAL=32
TIME_BUDGET=3600
HAL_TIME_BUDGET=3000          # per-cell cap for hal2maf in #1/#5 (was 600).  Raised so hal2maf
                              # at 1M (~4min) and 10M (~40min) can complete; bigger ones time out.
SBATCH_TIME=240               # 10 days (240h -> --time=240:00:00); long unattended run
SBATCH_MEM=512                # generous headroom: 32 concurrent cells, liftOver ~2-4 GB ea
TMP_GB=""
MAX_OUTPUT_BLOCKS=500           # #3 taffyBlockViz --maxOutputBlocks
PARTITION="long"
ACCOUNT=""
OUTROOT="$PWD/577_slide_bench"  # parent of the four per-comparison OUTDIRs

# -- Which comparisons to run (all on by default) --
DO_1=1; DO_2=1; DO_3=1; DO_4=1; DO_5=1; DO_6=1

# -- Mode --
SUBMIT=0                        # default DRY-RUN; --submit lets each driver actually sbatch
TEST=0                          # --test: tiny smoke config (see below)
NO_WAIT=0                       # pass --no-wait through to drivers when submitting
NO_STAGE=0                      # --no-stage-local: pass through to every driver (read
                                # inputs from network instead of copying to local scratch)

# Tool/bin env passthroughs (each driver also honours these from env).
# Resolved once here on the submit host and threaded into the jobs as
# absolute paths (the compute node shares the filesystem).  Set $TAFFY (or
# put the freshly-built taffy first on PATH) to control which binary runs.
TAFFY="${TAFFY:-$(command -v taffy || true)}"
# taffyBlockVizTest is built alongside taffy (same bin/), so derive it from
# the resolved taffy; blockVizBed comes from PATH (the hal bin dir).  No
# hardcoded absolute path -- that only ever matched one machine.
TAFFYBLOCKVIZ="${TAFFYBLOCKVIZ:-$([[ -n "$TAFFY" ]] && echo "${TAFFY%/*}/taffyBlockVizTest" || command -v taffyBlockVizTest 2>/dev/null || true)}"
BLOCKVIZBED="${BLOCKVIZBED:-$(command -v blockVizBed 2>/dev/null || true)}"

usage() {
    cat >&2 <<EOF
run_577_slide_bench.sh -- drive all four 577-way / hg38 slide comparisons

  --test            Tiny smoke config: 1 chrom, 2 small sizes (1000,100000),
                    2 species, short --timeBudget, small --maxOutputBlocks.
                    Implies --no-stage-local (reads inputs in place; a smoke
                    must not spend ~2h staging 4TB).  Still DRY-RUN unless you
                    also pass --submit.
  --submit          Let each driver actually submit (default: --dry-run only,
                    nothing is sent to SLURM).
  --no-wait         Pass --no-wait to each driver (submit + detach).
  --no-stage-local  Pass --no-stage-local to every driver: read inputs from
                    the network instead of copying them to local scratch.
                    Skips the big stage-in (fast functional smoke), but the
                    timings go disk-bound -- not for the real measurement run.
  --only LIST       Comma list of comparisons to run: 1,2,3,4,5,6 (default all).
                    5 = cmp5_view_sample (random genome-wide view bench);
                    6 = cmp6_depth_summary (depth bigWig vs bigMafSummary).
  -o DIR            Output root (default $OUTROOT); each comparison gets a
                    subdir cmp1_view / cmp2_baselift / cmp3_blockviz /
                    cmp4_chainlift / cmp5_view_sample / cmp6_depth_summary
                    under it.

  Input paths (REQUIRED for the comparisons you run):
    --baseTui FILE     base universal MAF (.uni.maf.gz + .tui)  [#1 #2 #4*]
    --baseTuiTaf FILE  OR base universal TAF (.uni.taf.gz + .tui) [#1 #2 #4*]
                       (*#4 uses it ONLY as the readable --statsSrc for
                        chrom enumeration; the chained .tui drives the lift)
    --chainedTui FILE  chained .tui (..._g10000.tui)            [#3 #4]
    --hal FILE         full 577-way .hal                        [#1 #2 #3]
    --hg38Maf FILE     hg38-anchored MAF (.maf.gz + .tai)       [#1]
    --bigbed FILE      bigmaf bigBed                            [#1]
    --depthBw FILE     universal-depth bigWig                   [#6]
    --summaryBb FILE   hg38 bigMafSummary bigBed                [#6]
    --chains DIR       UCSC chains dir                          [#2 #4]
    --tree FILE        species tree .nwk                        [#1 #2 #4 plots]
    --panel FILE       target panel -L (default $PANEL)

  Geometry / ladder overrides (sensible defaults baked in):
    --ref ID  --viewRefChrom NAME  --halGenome NAME  --halSeq NAME
    --bbChrom NAME  --blockvizChrom NAME  --viewStart INT  --viewMaxSize INT
    --blockvizSizes CSV  --liftSizes CSV  --nIntervals INT  --halMaxSize N
    --T INT  --timeBudget SEC  --halTimeBudget SEC  --time HRS  --mem GB
    --tmp GB  --maxOutputBlocks N  --partition X  --account X
  -h / --help

Every input path is cluster-side for the real run.  Known LOCAL hints
(for --dry-run smoke only):
  base .tui    /home/hickey/dev/work/unitaf/vgp-577way-v1.uni.taf.gz.tui
               (sibling of vgp-577way-v1.uni.taf.gz)
  chained .tui /home/hickey/dev/work/unitaf/vgp-577way-v1.uni.taf.gz.chained_g10000.tui
  tree         /home/hickey/dev/work/paffy-comp/vgp-577way.nwk
  subtree hals /home/hickey/dev/work/unitaf/vgp-577way-v1-MuridaeAnc3.hal
               /home/hickey/dev/work/unitaf/8-t2t-apes-2023v2.hal
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --test)           TEST=1; shift;;
        --submit)         SUBMIT=1; shift;;
        --no-wait)        NO_WAIT=1; shift;;
        --no-stage-local) NO_STAGE=1; shift;;
        --only)           ONLY="$2"; shift 2;;
        -o)               OUTROOT="$2"; shift 2;;
        --baseTui)        BASE_TUI_MAF="$2"; shift 2;;
        --baseTuiTaf)     BASE_TUI_TAF="$2"; shift 2;;
        --chainedTui)     CHAINED_TUI="$2"; shift 2;;
        --hal)            HAL="$2"; shift 2;;
        --hg38Maf)        HG38_MAF="$2"; shift 2;;
        --bigbed)         BIGBED="$2"; shift 2;;
        --depthBw)        DEPTH_BW="$2"; shift 2;;
        --summaryBb)      SUMMARY_BB="$2"; shift 2;;
        --chains)         CHAINS_DIR="$2"; shift 2;;
        --tree)           TREE="$2"; shift 2;;
        --panel)          PANEL="$2"; shift 2;;
        --ref)            REF="$2"; shift 2;;
        --viewRefChrom)   VIEW_REF_CHROM="$2"; shift 2;;
        --halGenome)      HAL_GENOME="$2"; shift 2;;
        --halSeq)         HAL_SEQ="$2"; shift 2;;
        --bbChrom)        BB_CHROM="$2"; shift 2;;
        --blockvizChrom)  BLOCKVIZ_CHROM="$2"; shift 2;;
        --viewStart)      VIEW_START="$2"; shift 2;;
        --viewNorm)       VIEW_NORM="$2"; shift 2;;
        --viewMaxSize)    VIEW_MAX_SIZE="$2"; shift 2;;
        --blockvizSizes)  BLOCKVIZ_SIZES="$2"; shift 2;;
        --liftSizes)      LIFT_SIZES="$2"; shift 2;;
        --nIntervals)     LIFT_N_INTERVALS="$2"; shift 2;;
        --halMaxSize)     HAL_MAX_SIZE="$2"; shift 2;;
        --T)              T_TOTAL="$2"; shift 2;;
        --timeBudget)     TIME_BUDGET="$2"; shift 2;;
        --halTimeBudget)  HAL_TIME_BUDGET="$2"; shift 2;;
        --time)           SBATCH_TIME="$2"; shift 2;;
        --mem)            SBATCH_MEM="$2"; shift 2;;
        --tmp)            TMP_GB="$2"; shift 2;;
        --maxOutputBlocks) MAX_OUTPUT_BLOCKS="$2"; shift 2;;
        --partition)      PARTITION="$2"; shift 2;;
        --account)        ACCOUNT="$2"; shift 2;;
        -h|--help)        usage 0;;
        *)                echo "unknown arg: $1" >&2; usage 1;;
    esac
done

# --only N,M,... selects which comparisons run.
if [[ -n "${ONLY:-}" ]]; then
    DO_1=0; DO_2=0; DO_3=0; DO_4=0; DO_5=0; DO_6=0
    IFS=',' read -r -a _only <<< "$ONLY"
    for c in "${_only[@]}"; do
        case "$c" in
            1) DO_1=1;; 2) DO_2=1;; 3) DO_3=1;; 4) DO_4=1;; 5) DO_5=1;; 6) DO_6=1;;
            *) echo "ERROR: --only takes 1..6 (got '$c')" >&2; exit 1;;
        esac
    done
fi

# --test: shrink everything to a 1-chrom, 2-size, 2-species smoke run.
if [[ "$TEST" -eq 1 ]]; then
    echo ">> --test mode: tiny smoke config (reads inputs in place -- no 4TB stage-in)" >&2
    NO_STAGE=1                  # a smoke must not spend ~2h staging 4TB; its tiny
                                # 1-100kb queries read fine straight off the shared FS
    BLOCKVIZ_SIZES="1000,100000"
    LIFT_SIZES="1000,100000"
    VIEW_SAMPLE_SIZES="1000,100000"   # #5 smoke: 2 small sizes
    VIEW_SAMPLE_N=2                    # #5 smoke: 2 random regions per size
    DEPTH_SUMMARY_SIZES="1000000,10000000"   # #6 smoke: 2 small zoom-out sizes
    DEPTH_SUMMARY_N=2                  # #6 smoke: 2 random regions per size
    VIEW_MAX_SIZE=100000        # cap #1's view ladder for the smoke (else it runs the
                                # whole-chrom log-decade ladder + needs the chrom length)
    LIFT_N_INTERVALS=5
    HAL_MAX_SIZE=100000
    TIME_BUDGET=300
    HAL_TIME_BUDGET=120
    MAX_OUTPUT_BLOCKS=50
    T_TOTAL=8
    SBATCH_MEM=32
    SBATCH_TIME=1               # a smoke run must not request 24 h
    # 2-species panel written to a temp file used as -L.
    TEST_PANEL="$OUTROOT/test_panel.tsv"
fi

mkdir -p "$OUTROOT"

# Two small target species for --test (must have chains + be in the tree).
if [[ "$TEST" -eq 1 ]]; then
    cat > "$TEST_PANEL" <<'EOF'
GCF_011100685.1	Canis_lupus_familiaris	Dog
GCF_016700215.2	Gallus_gallus	Chicken
EOF
    PANEL="$TEST_PANEL"
fi

# Driver paths.
VIEW_DRV="$HERE/taffy_view_bench_slurm.sh"
LIFT_DRV="$HERE/taffy_lift_bench_slurm.sh"
BLOCKVIZ_DRV="$HERE/taffy_blockviz_bench_slurm.sh"
SAMPLE_DRV="$HERE/taffy_view_sample_bench_slurm.sh"
DEPTH_SUMMARY_DRV="$HERE/uni_depth_summary_bench_slurm.sh"
for d in "$VIEW_DRV" "$LIFT_DRV" "$BLOCKVIZ_DRV" "$SAMPLE_DRV" "$DEPTH_SUMMARY_DRV"; do
    [[ -f "$d" ]] || { echo "ERROR: driver not found: $d" >&2; exit 1; }
done

# Common flags appended to every driver invocation.
COMMON_FLAGS=()
[[ "$SUBMIT" -eq 0 ]] && COMMON_FLAGS+=( --dry-run )
[[ "$NO_WAIT" -eq 1 ]] && COMMON_FLAGS+=( --no-wait )
[[ "$NO_STAGE" -eq 1 ]] && COMMON_FLAGS+=( --no-stage-local )
[[ -n "$PARTITION" ]]  && COMMON_FLAGS+=( --partition "$PARTITION" )
[[ -n "$ACCOUNT"   ]]  && COMMON_FLAGS+=( --account "$ACCOUNT" )
[[ -n "$TMP_GB"    ]]  && COMMON_FLAGS+=( --tmp "$TMP_GB" )

# Resolve the base universal source flag (-u MAF or --uniTaf TAF).
BASE_SRC_FLAG=()
if   [[ -n "$BASE_TUI_MAF" ]]; then BASE_SRC_FLAG=( -u "$BASE_TUI_MAF" )
elif [[ -n "$BASE_TUI_TAF" ]]; then BASE_SRC_FLAG=( --uniTaf "$BASE_TUI_TAF" )
fi

echo "================================================================"
echo " 577-way / hg38 slide benchmark wrapper"
echo "   mode:        $([[ $SUBMIT -eq 1 ]] && echo SUBMIT || echo 'DRY-RUN (no sbatch)')$([[ $TEST -eq 1 ]] && echo ' [--test]')"
echo "   out root:    $OUTROOT"
echo "   panel:       $PANEL"
echo "   halMaxSize:  $HAL_MAX_SIZE bp"
echo "   run:         #1=$DO_1 #2=$DO_2 #3=$DO_3 #4=$DO_4 #5=$DO_5 #6=$DO_6"
echo "================================================================"

# Helper: echo + run a driver invocation (array-safe).
run_driver() {
    echo
    echo ">>> $*"
    "$@"
}

# ======================================================================
# #1 -- MAF extraction (view): tui vs tai vs hal vs bigmaf.  Size ladder.
#       The view driver builds its OWN log-decade ladder internally and
#       runs every tool per size; taffy/tai/bb run the FULL ladder.  The
#       hal2maf cell is capped natively by --halMaxSize (added to the view
#       driver): sizes above the cap simply DON'T launch hal2maf -- no
#       run-then-timeout.  Single clean pass.
# ======================================================================
if [[ "$DO_1" -eq 1 ]]; then
    O1="$OUTROOT/cmp1_view"
    # MAF extraction: feed BOTH tui sources (maf- and taf-anchored) when
    # supplied, and force MAF output from every tool (--mafOnly) so the
    # comparison is apples-to-apples on MAF (drops the _taf / tai_taf cells).
    V1_SRC=()
    [[ -n "$BASE_TUI_MAF" ]] && V1_SRC+=( -u "$BASE_TUI_MAF" )
    [[ -n "$BASE_TUI_TAF" ]] && V1_SRC+=( --uniTaf "$BASE_TUI_TAF" )
    [[ ${#V1_SRC[@]} -gt 0 ]] || { echo "ERROR(#1): need --baseTui and/or --baseTuiTaf" >&2; exit 1; }
    for v in HG38_MAF HAL BIGBED; do
        [[ -n "${!v}" ]] || { echo "ERROR(#1): missing --$(echo "$v" | tr 'A-Z_' 'a-z')" >&2; exit 1; }
    done
    V1_FLAGS=( "${V1_SRC[@]}"
        -m "$HG38_MAF" -H "$HAL" -b "$BIGBED" -o "$O1"
        --refChrom "$VIEW_REF_CHROM" --halGenome "$HAL_GENOME"
        --halSeq "$HAL_SEQ" --bbChrom "$BB_CHROM" --start "$VIEW_START"
        --halMaxSize "$HAL_MAX_SIZE" --mafOnly
        -T "$T_TOTAL" --timeBudget "$TIME_BUDGET"
        --time "$SBATCH_TIME" --mem "$SBATCH_MEM" )
    [[ -n "$HAL_TIME_BUDGET" ]] && V1_FLAGS+=( --halTimeBudget "$HAL_TIME_BUDGET" )
    [[ -n "$VIEW_MAX_SIZE"   ]] && V1_FLAGS+=( --maxSize "$VIEW_MAX_SIZE" )
    [[ "$VIEW_NORM" -eq 0    ]] && V1_FLAGS+=( --no-norm )   # _norm cells deactivated this run
    run_driver env TAFFY="$TAFFY" bash "$VIEW_DRV" "${V1_FLAGS[@]}" "${COMMON_FLAGS[@]}"
    echo ">> #1: ${#V1_SRC[@]}-source MAF extraction (maf.tui${BASE_TUI_TAF:+ + taf.tui}) all -> MAF; hal2maf capped at --halMaxSize=$HAL_MAX_SIZE; taffy/tai/bb full ladder" >&2
fi

# ======================================================================
# #2 -- base-level lift: hal vs tui (FULL .tui, PER-COLUMN path), no liftOver.
#       The apples-to-apples base-resolution head-to-head -- halLiftover and
#       taffy lift both emit every aligned segment.  --column-walk runs the
#       per-column O(columns) path (the analog of hal's per-base walk);
#       --no-liftover drops the chain-block tool (that's the broad #4 baseline,
#       not base-level).  Small ladder (LIFT_PERCOL_SIZES, <=500kb) where both
#       per-column taffy and halLiftover are feasible; hal capped via --halMaxSize.
# ======================================================================
if [[ "$DO_2" -eq 1 ]]; then
    O2="$OUTROOT/cmp2_baselift"
    [[ ${#BASE_SRC_FLAG[@]} -gt 0 ]] || { echo "ERROR(#2): need --baseTui or --baseTuiTaf" >&2; exit 1; }
    for v in HAL CHAINS_DIR TREE; do
        [[ -n "${!v}" ]] || { echo "ERROR(#2): missing input for \$$v" >&2; exit 1; }
    done
    run_driver env TAFFY="$TAFFY" bash "$LIFT_DRV" "${BASE_SRC_FLAG[@]}" \
        -H "$HAL" -c "$CHAINS_DIR" -t "$TREE" -o "$O2" \
        -r "$REF" -L "$PANEL" -S "$LIFT_PERCOL_SIZES" -N "$LIFT_N_INTERVALS" \
        --column-walk --no-liftover --halMaxSize "${LIFT_PERCOL_SIZES##*,}" \
        -T "$T_TOTAL" --timeBudget "$TIME_BUDGET" \
        --time "$SBATCH_TIME" --mem "$SBATCH_MEM" "${COMMON_FLAGS[@]}"
    echo ">> #2 base-level: hal vs tui (full tui, per-column) into $O2 (-S '$LIFT_PERCOL_SIZES'); no liftOver; halLiftover <= --halMaxSize=$HAL_MAX_SIZE" >&2
fi

# ======================================================================
# #3 -- blockViz: taffyBlockViz vs halBlockViz.  NEW driver does the
#       per-size hal cap natively (--halMaxSize) -- single invocation.
# ======================================================================
if [[ "$DO_3" -eq 1 ]]; then
    O3="$OUTROOT/cmp3_blockviz"
    for v in CHAINED_TUI HAL; do
        [[ -n "${!v}" ]] || { echo "ERROR(#3): missing input for \$$v" >&2; exit 1; }
    done
    run_driver env TAFFYBLOCKVIZ="$TAFFYBLOCKVIZ" BLOCKVIZBED="$BLOCKVIZBED" bash "$BLOCKVIZ_DRV" \
        -u "$CHAINED_TUI" -H "$HAL" -L "$PANEL" -S "$BLOCKVIZ_SIZES" \
        --halMaxSize "$HAL_MAX_SIZE" -r "$BLOCKVIZ_CHROM" --tSpecies "$HAL_GENOME" \
        --start 0 --maxOutputBlocks "$MAX_OUTPUT_BLOCKS" \
        -o "$O3" -T "$T_TOTAL" --timeBudget "$TIME_BUDGET" \
        --time "$SBATCH_TIME" --mem "$SBATCH_MEM" "${COMMON_FLAGS[@]}"
    # blockViz has no --tree of its own; drop the tree into the cmp3 dir so
    # plot.py can label the panel's nearest/most-distant target (else it just
    # prints the genome count).
    [[ -n "$TREE" && -r "$TREE" ]] && cp -f "$TREE" "$O3/tree.nwk" 2>/dev/null \
        && echo ">> #3 tree -> $O3/tree.nwk (for plot.py panel annotation)" >&2
fi

# ======================================================================
# #4 -- chained liftover: taffy lift on the CHAINED .tui vs UCSC liftOver.
#       NO hal tool (liftOver is chain-based) -> --no-hal, full ladder,
#       one pass.  Chained .tui resolved via the symlink trick (see header).
# ======================================================================
if [[ "$DO_4" -eq 1 ]]; then
    O4="$OUTROOT/cmp4_chainlift"
    for v in CHAINED_TUI CHAINS_DIR TREE; do
        [[ -n "${!v}" ]] || { echo "ERROR(#4): missing input for \$$v" >&2; exit 1; }
    done
    [[ -f "$CHAINED_TUI" ]] || { echo "ERROR(#4): chained .tui not found: $CHAINED_TUI" >&2; exit 1; }
    # #4 needs a readable uni SOURCE for the lift driver's `taffy stats -s`
    # chrom enumeration (the chained-.tui stub is 0-byte and would crash it).
    # The base uni shares hg38's chrom set, so reuse it via --statsSrc.
    STATS_SRC=""
    if   [[ -n "$BASE_TUI_MAF" ]]; then STATS_SRC="$BASE_TUI_MAF"
    elif [[ -n "$BASE_TUI_TAF" ]]; then STATS_SRC="$BASE_TUI_TAF"
    fi
    [[ -n "$STATS_SRC" ]] || { echo "ERROR(#4): need --baseTui or --baseTuiTaf as the readable --statsSrc for chrom enumeration (the chained .tui has no readable source)" >&2; exit 1; }
    # Symlink so the lift driver's .tui-sibling resolution finds the chained
    # .tui without a MAF.  WORK dir holds the stub + symlink.
    WORK="$O4/_tui_link"
    mkdir -p "$WORK"
    CHAINED_STUB="$WORK/chained.uni.taf.gz"
    : > "$CHAINED_STUB"                       # 0-byte stub; taffy lift reads only the .tui
    ln -sf "$(readlink -f "$CHAINED_TUI")" "${CHAINED_STUB}.tui"
    echo ">> #4 chained-.tui resolution: stub $CHAINED_STUB + symlink ${CHAINED_STUB}.tui -> $CHAINED_TUI" >&2
    echo ">> #4 chrom enumeration via --statsSrc $STATS_SRC (readable; shares hg38 chroms)" >&2
    # --uniTaf so cells are labelled 'taf.tui' (the chained .tui is
    # TAF-anchored).  --statsSrc points chrom enumeration at the readable
    # base source.  --no-hal because #4's baseline is UCSC liftOver, not
    # halLiftover; -H is then optional.
    run_driver env TAFFY="$TAFFY" bash "$LIFT_DRV" \
        --uniTaf "$CHAINED_STUB" --statsSrc "$STATS_SRC" \
        -c "$CHAINS_DIR" -t "$TREE" -o "$O4" \
        -r "$REF" -L "$PANEL" -S "$LIFT_SIZES" -N "$LIFT_N_INTERVALS" \
        --no-hal -T "$T_TOTAL" --timeBudget "$TIME_BUDGET" \
        --time "$SBATCH_TIME" --mem "$SBATCH_MEM" "${COMMON_FLAGS[@]}"
fi

# ======================================================================
# #5 -- random-sample view bench: N random genome-wide hg38 regions per size,
#       mean+range over the samples (blockViz-style).  Four tools, SAME regions,
#       all -> hg38-anchored leaf MAF: taf.tui / tai / hal2maf / bigMafToMaf.
#       Complements #1 (fixed-location log-decade ladder) with a representative,
#       error-barred view at a few interesting sizes.  Reuses the same inputs;
#       hal2maf pushed to --halMaxSize (10M) so it appears at the bigger sizes.
# ======================================================================
if [[ "$DO_5" -eq 1 ]]; then
    O5="$OUTROOT/cmp5_view_sample"
    # taf.tui preferred (the user's "taf tui view"); fall back to the maf-backed
    # universal if only that was supplied.
    V5_SRC=""
    if   [[ -n "$BASE_TUI_TAF" ]]; then V5_SRC="$BASE_TUI_TAF"
    elif [[ -n "$BASE_TUI_MAF" ]]; then V5_SRC="$BASE_TUI_MAF"
    fi
    [[ -n "$V5_SRC" ]] || { echo "ERROR(#5): need --baseTuiTaf (or --baseTui) for the universal .tui source" >&2; exit 1; }
    for v in HG38_MAF HAL BIGBED; do
        [[ -n "${!v}" ]] || { echo "ERROR(#5): missing --$(echo "$v" | tr 'A-Z_' 'a-z')" >&2; exit 1; }
    done
    V5_FLAGS=( -u "$V5_SRC" -m "$HG38_MAF" -H "$HAL" -b "$BIGBED" -o "$O5"
        --refGenome "$HAL_GENOME"
        --sizes "$VIEW_SAMPLE_SIZES" --nSamples "$VIEW_SAMPLE_N" --seed "$VIEW_SAMPLE_SEED"
        --halMaxSize "$HAL_MAX_SIZE"
        -T "$T_TOTAL" --timeBudget "$TIME_BUDGET"
        --time "$SBATCH_TIME" --mem "$SBATCH_MEM" )
    [[ -n "$HAL_TIME_BUDGET" ]] && V5_FLAGS+=( --halTimeBudget "$HAL_TIME_BUDGET" )
    run_driver env TAFFY="$TAFFY" bash "$SAMPLE_DRV" "${V5_FLAGS[@]}" "${COMMON_FLAGS[@]}"
    echo ">> #5: random-sample view ($VIEW_SAMPLE_N regions/size, sizes $VIEW_SAMPLE_SIZES; taf.tui/tai/hal2maf/bigMafToMaf -> hg38 MAF); hal2maf to --halMaxSize=$HAL_MAX_SIZE" >&2
fi

# ======================================================================
# #6 -- zoom-out summary latency: depth bigWig lift vs bigMafSummary.
#       The high-zoom companion to #5 (1Mb-200Mb).  Times `taffy lift
#       --bigwig` on the CHAINED .tui over the universal-depth bigWig
#       ("ours") against bigMafSummary via bigBedToBed ("theirs"); wall
#       time only.  Reuses the chained .tui (#3/#4) for the lift and the
#       hg38 MAF (#1) as the region-sampling --chromSrc.
# ======================================================================
if [[ "$DO_6" -eq 1 ]]; then
    O6="$OUTROOT/cmp6_depth_summary"
    [[ -n "$CHAINED_TUI" ]] || { echo "ERROR(#6): need --chainedTui (the lift opens <it>.tui)" >&2; exit 1; }
    [[ -n "$DEPTH_BW"    ]] || { echo "ERROR(#6): need --depthBw (the universal-depth bigWig)" >&2; exit 1; }
    [[ -n "$SUMMARY_BB"  ]] || { echo "ERROR(#6): need --summaryBb (the hg38 bigMafSummary bigBed)" >&2; exit 1; }
    [[ -n "$HG38_MAF"    ]] || { echo "ERROR(#6): need --hg38Maf (the region-sampling --chromSrc)" >&2; exit 1; }
    V6_FLAGS=( -i "$CHAINED_TUI" -d "$DEPTH_BW" -b "$SUMMARY_BB"
        --chromSrc "$HG38_MAF" --refGenome "$HAL_GENOME" -o "$O6"
        --sizes "$DEPTH_SUMMARY_SIZES" --nSamples "$DEPTH_SUMMARY_N" --seed "$DEPTH_SUMMARY_SEED"
        -T "$T_TOTAL" --timeBudget "$TIME_BUDGET"
        --time "$SBATCH_TIME" --mem "$SBATCH_MEM" )
    run_driver env TAFFY="$TAFFY" bash "$DEPTH_SUMMARY_DRV" "${V6_FLAGS[@]}" "${COMMON_FLAGS[@]}"
    echo ">> #6: depth-summary (depth bigWig lift vs bigMafSummary; $DEPTH_SUMMARY_N regions/size, sizes $DEPTH_SUMMARY_SIZES) -> hg38 zoom-out" >&2
fi

echo
echo "================================================================"
echo " wrapper done ($([[ $SUBMIT -eq 1 ]] && echo 'submitted' || echo 'dry-run only -- nothing submitted'))."
echo " per-comparison outputs under: $OUTROOT"
echo "================================================================"
