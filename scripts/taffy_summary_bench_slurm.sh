#!/bin/bash
#
# taffy summary vs bigMafSummary benchmark -- SLURM
#
# Draws N random hg38 regions per size and, for each, runs BOTH:
#   ours    taffy summary --bins N --sampleCols S    (live per-pixel coverage from the .tui)
#   theirs  bigBedToBed   on the real bigMafSummary  (precomputed multiz-score summary)
# then compares: our coverage vs their score (mean |diff| + correlation, binned to a
# common grid), plus the wall time of each.  Aggregated mean+range over the samples.
#
# This validates `taffy summary` as a runtime bigMafSummary replacement at scale:
# does it match the precomputed summary (accuracy), and is it fast enough (timing)?
#
# Sampling: runtime, length-weighted over canonical hg38 chroms (chr1..22,X,Y),
# fixed seed, only chroms >= the size eligible (so 200Mb draws only from chr1/chr2).
# hg38 is named $REF_GENOME in the .tui (GCA_000001405.15.chrN); the bigMafSummary
# uses bare chrN; -r / bigBedToBed both take bare chrN.
#
# STAGING IS ON BY DEFAULT: taffy summary random-accesses ~pixel-count blocks, and
# the meaningful number is the LOCAL-NVMe latency the design targets -- reading those
# scattered blocks over the network measures network latency, not the index.  So the
# inputs are copied to local scratch (~486GB for 577: -u + .tui + summary.bb) first.
# --no-stage-local opts out (reads from the network; timing then network-bound).
#
# Usage:
#   taffy_summary_bench_slurm.sh -u UNI.taf.gz -b SUMMARY.bb -R GCA_000001405.15 \
#       -o OUTDIR [options]

set -euo pipefail

UNI=""                         # universal TAF/MAF (+.tui sibling) for taffy summary
BB=""                          # the bigMafSummary .bb (ground truth, bigBedToBed)
CHROM_SRC=""                   # hg38-anchored MAF/.tai (row-0 == hg38) for chrom lengths;
                               # the universal's stats -s only lists ancestor row-0 seqs
MAF_SRC=""                     # --maf: forward-extraction source if != UNI (e.g. local .tui + remote MAF)
OUTDIR=""
REF_GENOME="GCA_000001405.15"  # hg38 label in the .tui; bigMafSummary uses bare chrN
SIZES_CSV="1000000,10000000,50000000,100000000,200000000"   # 1,10,50,100,200 Mb
N_SAMPLES=10
SEED=20260620
N_BINS=1000                    # pixel bins per query
SAMPLE_COLS=100                # cols/anchor for the sampled coverage estimate.  Block DECODE
                               # dominates (profile via -t) and scales with columns, while
                               # accuracy is anchor-density-limited (S barely moves it), so a
                               # small S is near-free speed.  Lower (e.g. 50) for more, sweep to taste.
THREADS=16                     # taffy summary -T (parallel anchor reads)
T_TOTAL=32                     # concurrency thread-slot budget across cells
TIME_BUDGET=3600               # per-cell wall cap (timeout SIGKILL)
SBATCH_TIME=24
SBATCH_MEM=128
TMP_GB=""
STAGE_LOCAL=1                  # ON by default: meaningful timing needs LOCAL-disk reads
                               # (network random-read latency would dominate, not the design's
                               # local-NVMe target).  --no-stage-local to opt out.
PARTITION=""
ACCOUNT=""
DRY_RUN=0
WAIT=1
TAFFY="${TAFFY:-$(command -v taffy || true)}"
BIGBED2BED="${BIGBED2BED:-$(command -v bigBedToBed || true)}"

usage() {
    cat >&2 <<EOF
taffy_summary_bench_slurm.sh -- taffy summary vs bigMafSummary, N random hg38
                               regions per size, accuracy + timing, mean+range

Required:
  -u FILE       Universal TAF/MAF (.taf.gz/.maf.gz with .tui sibling) for taffy summary
  -b FILE       bigMafSummary .bb (ground truth)
  -o DIR        Output directory

  --chromSrc FILE  hg38-anchored MAF (row-0 == hg38, with .tai) used ONLY for
                   sampling chrom lengths via \`taffy stats -s\` -- the universal
                   -u lists only ancestor row-0 seqs, so this is required when -u
                   is universal (e.g. the view-bench's vgp-577way-v1-hg38.maf.gz).

Optional:
  -R NAME       hg38 label in the .tui (default $REF_GENOME).  Region/-r and
                bigBedToBed use the bare chrN (prefix stripped).
  --maf FILE    Forward-extraction source if it differs from -u (e.g. a local
                .tui for the load + a remote MAF for the block reads).
  --sizes CSV   Query sizes in bp (default $SIZES_CSV)
  --nSamples N  Random regions per size (default $N_SAMPLES)
  --seed N      RNG seed (default $SEED)
  --bins N      Pixel bins per query (default $N_BINS)
  --sampleCols S  Columns/anchor for the sampled coverage estimate (default $SAMPLE_COLS;
                  0 = exact full-window coverage -- expensive at 577)
  --threads N   taffy summary -T parallel anchor reads (default $THREADS)
  -T INT        Concurrency thread-slot budget across cells (default $T_TOTAL)
  --timeBudget SEC  Per-cell wall cap (default $TIME_BUDGET)
  --time HRS    sbatch wall (default $SBATCH_TIME)
  --mem GB      sbatch mem (default $SBATCH_MEM)
  --tmp GB      Per-task local scratch (sbatch --tmp=N).  Default unset.
  --no-stage-local  Read inputs from network (default: STAGE -u/.tui + -b to
                    \$TMPDIR for meaningful LOCAL-disk timing).  ~486GB for 577.
  --partition X --account X
  --no-wait     Submit and detach
  --dry-run     Print sbatch; do not submit
  -h            Help

Override binary paths via env: TAFFY, BIGBED2BED
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)                 UNI="$2"; shift 2;;
        -b)                 BB="$2"; shift 2;;
        --chromSrc)         CHROM_SRC="$2"; shift 2;;
        --maf)              MAF_SRC="$2"; shift 2;;
        -o)                 OUTDIR="$2"; shift 2;;
        -R)                 REF_GENOME="$2"; shift 2;;
        --sizes)            SIZES_CSV="$2"; shift 2;;
        --nSamples)         N_SAMPLES="$2"; shift 2;;
        --seed)             SEED="$2"; shift 2;;
        --bins)             N_BINS="$2"; shift 2;;
        --sampleCols)       SAMPLE_COLS="$2"; shift 2;;
        --threads)          THREADS="$2"; shift 2;;
        -T)                 T_TOTAL="$2"; shift 2;;
        --timeBudget)       TIME_BUDGET="$2"; shift 2;;
        --time)             SBATCH_TIME="$2"; shift 2;;
        --mem)              SBATCH_MEM="$2"; shift 2;;
        --tmp)              TMP_GB="$2"; shift 2;;
        --stage-local)      STAGE_LOCAL=1; shift;;
        --no-stage-local)   STAGE_LOCAL=0; shift;;
        --partition)        PARTITION="$2"; shift 2;;
        --account)          ACCOUNT="$2"; shift 2;;
        --no-wait)          WAIT=0; shift;;
        --dry-run)          DRY_RUN=1; shift;;
        -h|--help)          usage 0;;
        *)                  echo "unknown arg: $1" >&2; usage 1;;
    esac
done

for v in UNI BB OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: missing required input \$$v" >&2; usage 1; }
done
[[ -n "$TAFFY"      ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$BIGBED2BED" ]] || { echo "ERROR: bigBedToBed not on PATH (set \$BIGBED2BED)" >&2; exit 1; }
[[ -f "$UNI"        ]] || { echo "ERROR: $UNI not found" >&2; exit 1; }
[[ -f "${UNI}.tui"  ]] || { echo "ERROR: $UNI has no .tui sibling" >&2; exit 1; }
[[ -f "$BB"         ]] || { echo "ERROR: $BB not found" >&2; exit 1; }
[[ -z "$MAF_SRC" || -f "$MAF_SRC" ]] || { echo "ERROR: --maf $MAF_SRC not found" >&2; exit 1; }
[[ -n "$CHROM_SRC" ]] || CHROM_SRC="$UNI"   # default; warns below if it has no hg38.chrN
[[ -f "$CHROM_SRC" ]] || { echo "ERROR: --chromSrc $CHROM_SRC not found" >&2; exit 1; }
[[ "$SIZES_CSV" =~ ^[0-9]+(,[0-9]+)*$ ]] || { echo "ERROR: --sizes must be a CSV of integers" >&2; exit 1; }
[[ "$N_SAMPLES" =~ ^[0-9]+$ && "$N_SAMPLES" -ge 1 ]] || { echo "ERROR: --nSamples must be a positive integer" >&2; exit 1; }
[[ "$SEED" =~ ^[0-9]+$ ]] || { echo "ERROR: --seed must be a non-negative integer" >&2; exit 1; }

IFS=',' read -r -a SIZES <<< "$SIZES_CSV"

mkdir -p "$OUTDIR" "$OUTDIR/logs"
echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> bigBedToBed:   $BIGBED2BED"
echo ">> uni (.tui):    $UNI (+.tui)${MAF_SRC:+   extract-src: $MAF_SRC}"
echo ">> bigMafSummary: $BB"
echo ">> ref genome:    $REF_GENOME (region/-r + bigBedToBed use bare chrN)"
echo ">> sizes:         $SIZES_CSV bp"
echo ">> samples/size:  $N_SAMPLES   seed: $SEED"
echo ">> bins/query:    $N_BINS   sampleCols: $SAMPLE_COLS   summary threads: $THREADS"
echo ">> local-stage:   $([[ $STAGE_LOCAL -eq 1 ]] && echo ON || echo 'OFF (random-access)')"

# --- The region sampler (separate file, quoted heredoc -> no shell escaping).
SAMPLER="$OUTDIR/sample_regions.py"
cat > "$SAMPLER" <<'PYEOF'
#!/usr/bin/env python3
"""Sample N random hg38 regions per size, length-weighted over canonical chroms,
deterministically from a seed.  Reads a `taffy stats -s` dump (name<ws>length),
keeps `<refGenome>.chrN` canonical chroms (chr1..22,X,Y), prints
`size  chrom  start  end` (0-based half-open, bare chrN)."""
import argparse, random, re
ap = argparse.ArgumentParser()
ap.add_argument("--stats", required=True)
ap.add_argument("--refGenome", required=True)
ap.add_argument("--sizes", required=True)
ap.add_argument("--nSamples", type=int, required=True)
ap.add_argument("--seed", type=int, required=True)
a = ap.parse_args()
prefix = a.refGenome + "."
canon = re.compile(r"^chr([0-9]+|[XY])$")
chroms = []
with open(a.stats) as f:
    for line in f:
        p = line.split()
        if len(p) < 2 or not p[1].isdigit():
            continue
        if p[0].startswith(prefix) and canon.match(p[0][len(prefix):]):
            chroms.append((p[0][len(prefix):], int(p[1])))
if not chroms:
    raise SystemExit("sample_regions: no canonical %s.chrN sequences in %s" % (a.refGenome, a.stats))
sizes = [int(x) for x in a.sizes.split(",") if x]
rng = random.Random(a.seed)
for size in sizes:
    elig = [(c, L) for c, L in chroms if L >= size]
    if not elig:
        continue
    weights = [L - size + 1 for _, L in elig]
    for _ in range(a.nSamples):
        c, L = rng.choices(elig, weights=weights, k=1)[0]
        start = rng.randint(0, L - size)
        print("%d\t%s\t%d\t%d" % (size, c, start, start + size))
PYEOF
chmod +x "$SAMPLER"

# --- The per-region comparison (ours vs theirs, binned to a common grid). -----
COMPARE="$OUTDIR/compare.py"
cat > "$COMPARE" <<'PYEOF'
#!/usr/bin/env python3
"""Compare taffy summary coverage (ours) vs bigMafSummary score (theirs) over one
region.  Both are bed-like (chrom start end src val [L R]); we bin each to N equal
bins over [start,end) by interval MIDPOINT (theirs intervals can be ~10-50kb and
offset), average per (src,bin), join, and print one TSV line:
  n_sp_ours  n_sp_theirs  n_matched  mean_abs_diff  pearson_corr
ours val = coverage fraction; theirs val = multiz identity-ish score (so ours is
expected systematically >= theirs, gap growing with divergence)."""
import argparse, math
ap = argparse.ArgumentParser()
ap.add_argument("--ours", required=True)
ap.add_argument("--theirs", required=True)
ap.add_argument("--bins", type=int, required=True)
ap.add_argument("--start", type=int, required=True)
ap.add_argument("--end", type=int, required=True)
a = ap.parse_args()
span = max(1, a.end - a.start)
def load(path):
    acc = {}   # (src,bin) -> [sum,count]
    sp = set()
    try:
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                try:
                    s, e, src, v = int(p[1]), int(p[2]), p[3], float(p[4])
                except ValueError:
                    continue
                mid = (s + e) // 2
                if mid < a.start or mid >= a.end:
                    continue
                b = (mid - a.start) * a.bins // span
                if b < 0 or b >= a.bins:
                    continue
                k = (src, b)
                t = acc.get(k)
                if t is None:
                    acc[k] = [v, 1]
                else:
                    t[0] += v; t[1] += 1
                sp.add(src)
    except FileNotFoundError:
        pass
    return {k: t[0] / t[1] for k, t in acc.items()}, sp
ours, sp_o = load(a.ours)
theirs, sp_t = load(a.theirs)
keys = set(ours) & set(theirs)
n = len(keys)
if n == 0:
    print("%d\t%d\t0\tNA\tNA" % (len(sp_o), len(sp_t)))
    raise SystemExit(0)
so = st = sd = soo = stt = sot = 0.0
for k in keys:
    o, t = ours[k], theirs[k]
    so += o; st += t; sd += abs(o - t); soo += o * o; stt += t * t; sot += o * t
mo, mt = so / n, st / n
denom = math.sqrt(max(0.0, soo / n - mo * mo)) * math.sqrt(max(0.0, stt / n - mt * mt))
corr = (sot / n - mo * mt) / denom if denom > 1e-12 else 0.0
print("%d\t%d\t%d\t%.4f\t%.4f" % (len(sp_o), len(sp_t), n, sd / n, corr))
PYEOF
chmod +x "$COMPARE"

# --- Generate the runner script (the thing sbatch executes). -----------
RUNNER="$OUTDIR/bench.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -uo pipefail

UNI="$UNI"
BB="$BB"
CHROM_SRC="$CHROM_SRC"
MAF_SRC="$MAF_SRC"
OUTDIR="$OUTDIR"
REF_GENOME="$REF_GENOME"
SIZES_CSV="$SIZES_CSV"
N_SAMPLES=$N_SAMPLES
SEED=$SEED
N_BINS=$N_BINS
SAMPLE_COLS=$SAMPLE_COLS
THREADS=$THREADS
T_TOTAL=$T_TOTAL
TIME_BUDGET=$TIME_BUDGET
STAGE_LOCAL=$STAGE_LOCAL
TAFFY="$TAFFY"
BIGBED2BED="$BIGBED2BED"
SIZES=( ${SIZES[*]} )

BENCH_TSV="\$OUTDIR/bench.tsv"
LOGDIR="\$OUTDIR/logs"
mkdir -p "\$LOGDIR"

# --- Optional stage-in (OFF by default; summary is random-access). ----
if [[ "\$STAGE_LOCAL" -eq 1 ]]; then
    SCRATCH="\${TMPDIR:-/tmp/taffy_summary_\${SLURM_JOB_ID:-\$\$}}"
    STAGE_DIR="\$SCRATCH/taffy_summary_stage_\${SLURM_JOB_ID:-\$\$}"
    mkdir -p "\$STAGE_DIR"
    trap 'rm -rf "\$STAGE_DIR" 2>/dev/null || true' EXIT
    stage_pids=()
    stage_bg() { local s="\$1"; [[ -n "\$s" && -f "\$s" ]] || return 0; cp "\$s" "\$STAGE_DIR/\$(basename "\$s")" & stage_pids+=( \$! ); }
    stage_srcs=( "\$UNI" "\$UNI.tui" "\$BB" )
    [[ -n "\$MAF_SRC" ]] && stage_srcs+=( "\$MAF_SRC" )
    for f in "\${stage_srcs[@]}"; do stage_bg "\$f"; done
    src_rc=0; for p in "\${stage_pids[@]}"; do wait "\$p" || src_rc=1; done
    [[ "\$src_rc" -eq 0 ]] || { echo "ERROR: stage-in cp returned non-zero" >&2; exit 1; }
    # VERIFY each staged copy matches its source byte count.  A cp to a full FS
    # can return 0 yet truncate (the write succeeds, the flush doesn't), and the
    # .tui/.tai-based reads then fail cryptically downstream -- catch it here.
    for f in "\${stage_srcs[@]}"; do
        [[ -f "\$f" ]] || continue
        ss=\$(stat -Lc %s "\$f" 2>/dev/null); ds=\$(stat -c %s "\$STAGE_DIR/\$(basename "\$f")" 2>/dev/null || echo -1)
        [[ "\$ss" == "\$ds" ]] || { echo "ERROR: staged \$(basename "\$f") is \$ds bytes, source is \$ss -- truncated (scratch full?)" >&2; exit 1; }
    done
    # NB: no index-mtime touch needed here -- taffy summary reads only the .tui
    # (tui_load does no staleness check), and the sampling stats -s runs on the
    # UNSTAGED --chromSrc.  The .tai staleness issue is the view-sample bench's.
    UNI="\$STAGE_DIR/\$(basename "\$UNI")"
    BB="\$STAGE_DIR/\$(basename "\$BB")"
    [[ -n "\$MAF_SRC" ]] && MAF_SRC="\$STAGE_DIR/\$(basename "\$MAF_SRC")"
    echo "stage: inputs staged + size-verified to \$STAGE_DIR" >&2
fi

# --- Sample the regions once. ------------------------------------------
CHROM_STATS="\$OUTDIR/hg38.stats.txt"
REGIONS="\$OUTDIR/regions.tsv"
"\$TAFFY" stats -s -i "\$CHROM_SRC" > "\$CHROM_STATS" 2> "\$OUTDIR/stats.err" || { echo "ERROR: taffy stats -s failed on \$CHROM_SRC:" >&2; cat "\$OUTDIR/stats.err" >&2; exit 1; }
python3 "\$OUTDIR/sample_regions.py" --stats "\$CHROM_STATS" --refGenome "\$REF_GENOME" \\
    --sizes "\$SIZES_CSV" --nSamples "\$N_SAMPLES" --seed "\$SEED" > "\$REGIONS" \\
    || { echo "ERROR: region sampling failed" >&2; exit 1; }
echo "sampled \$(wc -l < "\$REGIONS") regions (\$N_SAMPLES/size, seed \$SEED) over \$SIZES_CSV" >&2

printf "size_bp\tchrom\tstart\tsample\tour_wall\tour_rss_kb\tour_bytes\ttheir_wall\ttheir_bytes\tn_sp_ours\tn_sp_theirs\tn_matched\tmean_abs_diff\tcorr\tour_exit\ttheir_exit\n" > "\$BENCH_TSV"

MAF_FLAG=()
[[ -n "\$MAF_SRC" ]] && MAF_FLAG=( --maf "\$MAF_SRC" )

# Run one region: taffy summary (ours, timed) + bigBedToBed (theirs, timed) + compare.
run_region() {
    local N="\$1" chrom="\$2" s="\$3" e="\$4" sample="\$5"
    local stem="\${N}_\${sample}"
    local ours="\$LOGDIR/ours_\${stem}.tsv"  theirs="\$LOGDIR/theirs_\${stem}.bed"
    local ot="\$LOGDIR/ot_\${stem}.txt"      tt="\$LOGDIR/tt_\${stem}.txt"

    /usr/bin/time -q -f '%e %M' -o "\$ot" \\
        timeout --signal=KILL "\$TIME_BUDGET" \\
        "\$TAFFY" summary -t --bins "\$N_BINS" --sampleCols "\$SAMPLE_COLS" -T "\$THREADS" \\
            -R "\$REF_GENOME" -r "\${chrom}:\${s}-\${e}" -i "\$UNI" "\${MAF_FLAG[@]}" \\
        > "\$ours" 2> "\$LOGDIR/oerr_\${stem}.log"
    local orc=\$?
    /usr/bin/time -q -f '%e %M' -o "\$tt" \\
        timeout --signal=KILL "\$TIME_BUDGET" \\
        "\$BIGBED2BED" -chrom="\$chrom" -start="\$s" -end="\$e" "\$BB" stdout \\
        > "\$theirs" 2> "\$LOGDIR/terr_\${stem}.log"
    local trc=\$?

    local ow orss tw
    read -r ow orss < "\$ot" 2>/dev/null || { ow=NA; orss=NA; }
    read -r tw _   < "\$tt" 2>/dev/null || tw=NA
    [[ "\$ow"  =~ ^[0-9.]+\$ ]] || ow=NA
    [[ "\$tw"  =~ ^[0-9.]+\$ ]] || tw=NA
    [[ "\$orss" =~ ^[0-9]+\$  ]] || orss=NA
    local ob tb
    ob=\$(wc -c < "\$ours" 2>/dev/null || echo 0)
    tb=\$(wc -c < "\$theirs" 2>/dev/null || echo 0)

    local cmp
    cmp=\$(python3 "\$OUTDIR/compare.py" --ours "\$ours" --theirs "\$theirs" --bins "\$N_BINS" --start "\$s" --end "\$e" 2> "\$LOGDIR/cmperr_\${stem}.log")
    [[ -n "\$cmp" ]] || cmp=\$'0\t0\t0\tNA\tNA'

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n" \\
        "\$N" "\$chrom" "\$s" "\$sample" "\$ow" "\$orss" "\$ob" "\$tw" "\$tb" "\$cmp" "\$orc" "\$trc"
    rm -f "\$ours" "\$theirs"   # keep bench.tsv compact; raw outputs are large
}

# --- Concurrency throttle (each region uses THREADS slots). ------------
launched=0; declare -A pid_threads
acquire_slot() {
    local th=\$1; (( th > T_TOTAL )) && th=\$T_TOTAL
    while (( launched > 0 && launched + th > T_TOTAL )); do
        wait -n 2>/dev/null || true
        local p; for p in "\${!pid_threads[@]}"; do
            kill -0 \$p 2>/dev/null || { launched=\$(( launched - pid_threads[\$p] )); unset pid_threads[\$p]; }
        done
    done
    launched=\$(( launched + th ))
}

for N in "\${SIZES[@]}"; do
    echo "=== size wave N=\$N ===" >&2
    t_wave=\$SECONDS
    unset pids rowfiles; declare -A pids rowfiles   # CLEAR stale keys: a declare-A on an
                                                    # existing assoc array does NOT reset it, so a
                                                    # smaller/empty later wave would re-append a
                                                    # prior wave's rowfiles.  Match pid_threads below.
    launched=0; unset pid_threads; declare -A pid_threads
    sample=0
    while IFS=\$'\t' read -r rsize chrom rs re; do
        [[ "\$rsize" == "\$N" ]] || continue
        sample=\$(( sample + 1 ))
        rf="\$LOGDIR/row_\${N}_\${sample}.tsv"
        rowfiles["\$sample"]="\$rf"
        acquire_slot "\$THREADS"
        ( run_region "\$N" "\$chrom" "\$rs" "\$re" "\$sample" ) > "\$rf" &
        pids["\$sample"]=\$!; pid_threads[\$!]=\$THREADS
    done < "\$REGIONS"
    for k in "\${!pids[@]}"; do
        wait "\${pids[\$k]}" || true
        [[ -s "\${rowfiles[\$k]}" ]] && cat "\${rowfiles[\$k]}" >> "\$BENCH_TSV"
    done
    echo "=== size wave N=\$N took \$((SECONDS - t_wave)) s ===" >&2
done
echo "bench done.  TSV: \$BENCH_TSV" >&2
EOF
chmod +x "$RUNNER"

# --- Companion plot: timing (ours vs theirs) + accuracy + correlation. -----
PLOT="$OUTDIR/plot.py"
cat > "$PLOT" <<'PY'
#!/usr/bin/env python3
"""taffy summary vs bigMafSummary: mean+range over the sampled regions vs query
size.  Panel 1: wall time (ours = live coverage; theirs = bigBedToBed lookup).
Panel 2: accuracy = mean |our_coverage - their_score| (lower is closer; a floor
is expected since ours is coverage, theirs multiz identity).  Panel 3: per-bin
Pearson correlation."""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

bd = os.path.dirname(os.path.abspath(__file__))
def fnum(x):
    try: return float(x)
    except (ValueError, TypeError): return None
rows = []
with open(os.path.join(bd, "bench.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        r["size_bp"] = int(r["size_bp"])
        for k in ("our_wall", "their_wall", "mean_abs_diff", "corr"):
            r[k] = fnum(r.get(k))
        rows.append(r)

def agg(field, pred=lambda r: True):
    by = {}
    for r in rows:
        if r[field] is None or not pred(r):
            continue
        by.setdefault(r["size_bp"], []).append(r[field])
    return [(s, sum(v)/len(v), min(v), max(v)) for s, v in sorted(by.items())]

def draw(ax, series, title, ylab, logy=True):
    for pts, c, lab in series:
        if not pts: continue
        xs=[p[0] for p in pts]; m=[p[1] for p in pts]; lo=[p[2] for p in pts]; hi=[p[3] for p in pts]
        ax.plot(xs, m, "o-", color=c, label=lab); ax.fill_between(xs, lo, hi, color=c, alpha=0.15)
    ax.set_xscale("log")
    if logy: ax.set_yscale("log")
    ax.set_xlabel("query size (bp)"); ax.set_ylabel(ylab); ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y,_: f"{y:g}"))
    if series and any(s[0] for s in series): ax.legend(fontsize=8)

fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(21, 6))
fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.12, wspace=0.24)
draw(a1, [(agg("our_wall"), "#1f77b4", "taffy summary (live)"),
          (agg("their_wall"), "#ff7f0e", "bigBedToBed (precomputed)")],
     "wall time", "seconds")
draw(a2, [(agg("mean_abs_diff"), "#9467bd", "mean |coverage - score|")],
     "accuracy vs bigMafSummary", "mean abs diff", logy=False)
draw(a3, [(agg("corr"), "#2ca02c", "per-bin Pearson r")],
     "correlation vs bigMafSummary", "r", logy=False)
ns = max((sum(1 for r in rows if r["size_bp"] == s) for s in {r["size_bp"] for r in rows}), default=0)
fig.suptitle(f"taffy summary vs bigMafSummary -- 577-way hg38, random regions, "
             f"mean ± range over up to {ns} samples/size", fontsize=12, y=0.98)
out = os.path.join(bd, "bench.png"); fig.savefig(out, dpi=140); print(f"wrote {out}")
PY
chmod +x "$PLOT"

# --- Submit. -----------------------------------------------------------
SBATCH_ARGS=(
    --cpus-per-task="$T_TOTAL"
    --mem="${SBATCH_MEM}G"
    --time="${SBATCH_TIME}:00:00"
    --output="$OUTDIR/slurm_%j.log"
    --error="$OUTDIR/slurm_%j.err.log"
    -J taffy_summary_bench
)
[[ -n "$TMP_GB"    ]] && SBATCH_ARGS+=( --tmp="${TMP_GB}G" )
[[ -n "$PARTITION" ]] && SBATCH_ARGS+=( --partition="$PARTITION" )
[[ -n "$ACCOUNT"   ]] && SBATCH_ARGS+=( --account="$ACCOUNT" )

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo ">> DRY RUN -- would submit:"
    echo "sbatch ${SBATCH_ARGS[*]} --parsable $RUNNER"
else
    echo ">> submitting..."
    JOB=$(sbatch "${SBATCH_ARGS[@]}" --parsable "$RUNNER")
    echo ">> job id: $JOB"
    echo ">> stdout: $OUTDIR/slurm_${JOB}.log   stderr: $OUTDIR/slurm_${JOB}.err.log"
fi

if [[ "$DRY_RUN" -ne 1 && "$WAIT" -eq 1 ]]; then
    echo ">> waiting for job $JOB ..."
    while squeue -j "$JOB" -h -o "%T" 2>/dev/null | grep -qE "PENDING|RUNNING|CONFIGURING|COMPLETING|RESIZING|SUSPENDED|REQUEUED"; do
        sleep 60
    done
    FINAL_STATE=$(sacct -j "$JOB" -X -n -o State 2>/dev/null | head -1 | tr -d ' ')
    echo ">> job $JOB final state: ${FINAL_STATE:-UNKNOWN}"
    case "$FINAL_STATE" in COMPLETED) ;; *) echo ">> NON-SUCCESS -- check $OUTDIR/slurm_${JOB}.err.log"; exit 1;; esac
fi

echo ">> done.  after completion: python3 $OUTDIR/plot.py   # writes bench.png"
