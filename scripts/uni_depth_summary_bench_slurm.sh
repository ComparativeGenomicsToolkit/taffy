#!/bin/bash
#
# Universal-depth SUMMARY (zoom-out) benchmark driver -- SLURM
#
# The high-zoom companion to taffy_view_sample_bench_slurm.sh.  Where that one
# benches FULL-MAF browser queries at 100bp..1Mb (taffy view / hal2maf /
# bigMafToMaf), this one benches the SUMMARY / depth-track fallback that a MAF
# browser serves at chromosome scale (1Mb..200Mb): our universal-depth bigWig
# lifted to the reference via the chained .tui, versus UCSC bigMafSummary.
#
# Two tools, each answering the SAME sampled hg38 regions, WALL TIME ONLY:
#   ours    taffy lift --bigwig   (PER-SPECIES vector bigWig -> all-N per-species
#                                  coverage in ref coords, via the chained .tui)
#   theirs  bigBedToBed           (precomputed hg38 bigMafSummary .bb)
#
# Both now serve PER-SPECIES coverage at chromosome scale (ours all-N in one
# matrix, theirs one row per species), so this is the apples-to-apples zoom-out
# comparison.  It is a LATENCY comparison: the VALUES use different metrics
# (ours = coverage fraction, bigMafSummary = scaled-multiz identity), so it is
# not an accuracy/correlation panel.  (A scalar total-depth bigWig also works but
# lifts as a single track.)
#
# Sampling happens AT RUNTIME on the cluster: N random genome-wide hg38 regions
# per size, length-weighted over the canonical chroms (chr1..22,X,Y), drawn once
# from a fixed seed so both tools see the SAME regions.  Per cell we record wall
# seconds + exit + timed-out flag, tagged with (tool,size,chrom,start,sample).
#
# This dataset names hg38 as GCA_000001405.15: the .tui coordinate space is
# `GCA_000001405.15.chrN`, the lift region is GCA_000001405.15.chrN:s-e
# (--genome GCA_000001405.15), and the bigMafSummary uses bare `chrN`.
#
# NB on -i: `taffy lift --bigwig` opens `<-i>.tui` (it appends `.tui`), and in
# bigWig mode it NEVER touches the .taf.gz itself -- only the .tui + the .bw.  To
# serve the ZOOM-OUT track we point it at the CHAINED .tui
# (<uni>.taf.gz.chained_gN.tui).  Pass -i as that .tui file (or its stem); the
# driver strips a trailing `.tui` before handing it to taffy, since taffy
# re-appends it.
#
# Usage:
#   uni_depth_summary_bench_slurm.sh \
#       -i UNI.taf.gz.chained_g10000.tui  -d DEPTH.bw  -b SUMMARY.bb \
#       --chromSrc HG38.maf.gz  -o OUTDIR  [options]

set -euo pipefail

ITUI=""                       # chained .tui (zoom-out) OR its stem -> `taffy lift --bigwig -i <stem>`
DEPTH_BW=""                   # PER-SPECIES vector bigWig (taffy depth --perSpecies --bigwig;
                              # one component/leaf, single uni0 axis) + a <bw>.names sidecar.
                              # A scalar total-depth bigWig also works but lifts as 1 track.
BB=""                         # precomputed hg38 bigMafSummary .bb (bigBedToBed)
CHROM_SRC=""                  # hg38-anchored MAF/TAF for `taffy stats -s` region sampling
OUTDIR=""
T_TOTAL=32                    # concurrency = max simultaneous regions (1 slot/region)
REF_GENOME="GCA_000001405.15" # the hg38 label in this dataset (.tui genome, lift --genome)
# Default sizes (bp): 11 LINEAR points 1Mb..200Mb (the zoom-out range).
SIZES_CSV="1000000,20000000,40000000,60000000,80000000,100000000,120000000,140000000,160000000,180000000,200000000"
N_SAMPLES=10                  # random regions drawn per size
SEED=20260620                 # fixed RNG seed -> reproducible, same regions for both tools
TIME_BUDGET=3600              # per-cell wall seconds (timeout sends SIGKILL)
SBATCH_TIME=24
SBATCH_MEM=64
TMP_GB=""
STAGE_LOCAL=1                 # stage the .tui/.bw/.bb to local scratch (random-access seeks)
PARTITION=""
ACCOUNT=""
DRY_RUN=0
WAIT=1
TAFFY="${TAFFY:-$(command -v taffy || true)}"
BIGBEDTOBED="${BIGBEDTOBED:-$(command -v bigBedToBed || true)}"

usage() {
    cat >&2 <<EOF
uni_depth_summary_bench_slurm.sh -- bench zoom-out (summary) wall / RSS / output:
   universal-depth bigWig lift  vs  UCSC bigMafSummary,  N random regions/size

Required:
  -i FILE       Chained .tui (zoom-out) OR its stem.  taffy lift opens <stem>.tui;
                pass e.g. UNI.taf.gz.chained_g10000.tui (the .tui is stripped).
  -d FILE       Per-species vector bigWig (taffy depth --perSpecies --bigwig; needs
                its <bw>.names sidecar -> lift emits all-N).  A scalar total-depth
                bigWig also works but lifts as 1 track (NOT apples-to-apples).
  -b FILE       Precomputed hg38 bigMafSummary .bb (for bigBedToBed)
  --chromSrc FILE  hg38-anchored MAF/TAF for region sampling (taffy stats -s)
  -o DIR        Output directory

Optional:
  -T INT        Max simultaneous regions (1 slot/region; default $T_TOTAL)
  --refGenome NAME  hg38 label in the .tui (default $REF_GENOME).  lift region =
                    \`NAME.chrN:s-e\` --genome NAME; bigMafSummary uses bare chrN.
  --sizes CSV   Query sizes in bp (default: 11 linear points 1Mb..200Mb)
  --nSamples INT  Random regions drawn per size (default $N_SAMPLES)
  --seed INT    RNG seed for the draw (default $SEED).  Same seed -> same regions;
                both tools always see the same regions.
  --timeBudget SEC  Per-cell wall cap (default $TIME_BUDGET)
  --time HRS    sbatch wall (default $SBATCH_TIME)
  --mem GB      sbatch mem (default $SBATCH_MEM)
  --tmp GB      Per-task local scratch (sbatch --tmp=N).  Default unset.
  --no-stage-local  Read inputs from network paths (no copy to \$TMPDIR)
  --partition X --account X
  --no-wait     Submit and detach (default: block until SLURM finishes)
  --dry-run     Print sbatch; do not submit (missing inputs warn, don't error)
  -h            Help

Override binary paths via env: TAFFY, BIGBEDTOBED

Bin width per region is auto: BIN = max(1000, size/1000) -> ~1000 bins/region.
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)                 ITUI="$2"; shift 2;;
        -d)                 DEPTH_BW="$2"; shift 2;;
        -b)                 BB="$2"; shift 2;;
        --chromSrc)         CHROM_SRC="$2"; shift 2;;
        -o)                 OUTDIR="$2"; shift 2;;
        -T)                 T_TOTAL="$2"; shift 2;;
        --refGenome)        REF_GENOME="$2"; shift 2;;
        --sizes)            SIZES_CSV="$2"; shift 2;;
        --nSamples)         N_SAMPLES="$2"; shift 2;;
        --seed)             SEED="$2"; shift 2;;
        --timeBudget)       TIME_BUDGET="$2"; shift 2;;
        --time)             SBATCH_TIME="$2"; shift 2;;
        --mem)              SBATCH_MEM="$2"; shift 2;;
        --tmp)              TMP_GB="$2"; shift 2;;
        --no-stage-local)   STAGE_LOCAL=0; shift;;
        --partition)        PARTITION="$2"; shift 2;;
        --account)          ACCOUNT="$2"; shift 2;;
        --no-wait)          WAIT=0; shift;;
        --dry-run)          DRY_RUN=1; shift;;
        -h|--help)          usage 0;;
        *)                  echo "unknown arg: $1" >&2; usage 1;;
    esac
done

for v in ITUI DEPTH_BW BB CHROM_SRC OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: missing required input \$$v" >&2; usage 1; }
done
[[ -n "$TAFFY"       ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$BIGBEDTOBED" ]] || { echo "WARN: bigBedToBed not found (set \$BIGBEDTOBED); the 'theirs' tool will fail" >&2; }

# Resolve the .tui that lift will open.  `taffy lift -i X` opens `X.tui`, so the
# value we pass must be the .tui MINUS its extension.  Accept either form from
# the user: the .tui file itself, or its stem.
if [[ "$ITUI" == *.tui ]]; then
    ITUI_FILE="$ITUI"
else
    ITUI_FILE="${ITUI}.tui"
fi

# Input existence: hard error on a real run, soft warn under --dry-run (so the
# sbatch line can be previewed anywhere, including with cluster-only inputs).
check_input() {
    local path="$1" what="$2"
    if [[ ! -f "$path" ]]; then
        if [[ "$DRY_RUN" -eq 1 ]]; then
            echo "WARN (dry-run): $what not found: $path" >&2
        else
            echo "ERROR: $what not found: $path" >&2; exit 1
        fi
    fi
}
check_input "$ITUI_FILE" "chained .tui"
check_input "$DEPTH_BW"  "depth bigWig"
check_input "$BB"        "bigMafSummary .bb"
check_input "$CHROM_SRC" "--chromSrc (sampling source)"

# The per-species vector bigWig (taffy depth --perSpecies --bigwig) carries a
# <bw>.names sidecar (one leaf name per line, gerp order) that the vector lift
# reads for the species legend.  Require it for the vector case so a missing
# sidecar fails here, not mid-run.  A scalar total-depth bigWig (BIGWIG64_MAGIC,
# no .names) lifts as a single track and is allowed (but is NOT apples-to-apples).
if [[ -f "$DEPTH_BW" ]]; then
    case "$(od -An -tx1 -N4 "$DEPTH_BW" 2>/dev/null | tr -d ' \n')" in
        65fc8f88)  # BIGWIG64VEC_MAGIC -> per-species vector
            [[ -f "$DEPTH_BW.names" ]] || { echo "ERROR: per-species vector bigWig is missing its sidecar: $DEPTH_BW.names" >&2; exit 1; }
            echo ">> depth bigWig: PER-SPECIES vector, $(wc -l < "$DEPTH_BW.names") components (lift emits all-N -> apples-to-apples vs bigMafSummary)";;
        64fc8f88)  # BIGWIG64_MAGIC -> scalar total-depth
            echo ">> depth bigWig: SCALAR total-depth (1 track); NOT the apples-to-apples per-species comparison";;
        *) echo ">> WARN: $DEPTH_BW magic unrecognized; expected a 64-bit (vector or scalar) bigWig" >&2;;
    esac
fi

[[ "$SIZES_CSV" =~ ^[0-9]+(,[0-9]+)*$ ]] || { echo "ERROR: --sizes must be a CSV of integers (got '$SIZES_CSV')" >&2; exit 1; }
[[ "$N_SAMPLES" =~ ^[0-9]+$ && "$N_SAMPLES" -ge 1 ]] || { echo "ERROR: --nSamples must be a positive integer" >&2; exit 1; }
[[ "$SEED" =~ ^[0-9]+$ ]] || { echo "ERROR: --seed must be a non-negative integer" >&2; exit 1; }
[[ "$T_TOTAL" =~ ^[0-9]+$ && "$T_TOTAL" -ge 1 ]] || { echo "ERROR: -T must be a positive integer" >&2; exit 1; }

# Bash size array from the CSV (drives the per-size waves; the sampler reads the
# CSV form directly).
IFS=',' read -r -a SIZES <<< "$SIZES_CSV"

mkdir -p "$OUTDIR" "$OUTDIR/logs"
echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> bigBedToBed:   ${BIGBEDTOBED:-(missing)}"
echo ">> chained .tui:  $ITUI_FILE  (lift -i ${ITUI_FILE%.tui})"
echo ">> depth bigWig:  $DEPTH_BW"
echo ">> bigMafSummary: $BB"
echo ">> chrom source:  $CHROM_SRC  (taffy stats -s for region draw)"
echo ">> ref genome:    $REF_GENOME (.tui NAME.chrN, lift --genome; bb bare chrN)"
echo ">> sizes:         $SIZES_CSV bp"
echo ">> samples/size:  $N_SAMPLES   seed: $SEED"
echo ">> bin width:     auto = max(1000, size/1000) bp (~1000 bins/region)"
echo ">> concurrency:   <= $T_TOTAL regions at once (1 slot/region)"
echo ">> time budget:   $TIME_BUDGET s per cell"
echo ">> local-stage:   $([[ $STAGE_LOCAL -eq 1 ]] && echo "ON (copies .tui/.bw/.bb to \$TMPDIR)" || echo "OFF (reads from network)")"

if [[ "$STAGE_LOCAL" -eq 1 ]]; then
    STAGE_BYTES=0
    for f in "$ITUI_FILE" "$DEPTH_BW" "$DEPTH_BW.names" "$BB"; do
        [[ -f "$f" ]] && STAGE_BYTES=$(( STAGE_BYTES + $(stat -Lc %s "$f" 2>/dev/null || echo 0) ))
    done
    STAGE_GB=$(( STAGE_BYTES / (1024**3) ))
    echo ">> stage-in size: ~${STAGE_GB} GB total"
    [[ -n "$TMP_GB" ]] || echo ">> --tmp:         not requested (pass --tmp $(( STAGE_GB + 20 )) only if your cluster enforces it)"
fi

# --- The region sampler (separate file, quoted heredoc -> no shell escaping).
#     Identical to the view-sample bench sampler: same draw given the same seed.
SAMPLER="$OUTDIR/sample_regions.py"
cat > "$SAMPLER" <<'PYEOF'
#!/usr/bin/env python3
"""Sample N random genome-wide hg38 regions per size, length-weighted over the
canonical chroms, deterministically from a seed.  Reads a `taffy stats -s` dump
(name<ws>length per line), keeps `<refGenome>.chrN` canonical chroms (chr1..22,
X,Y -- no alt/unplaced contigs, which carry underscores and may be absent from
the bigMafSummary), and prints `size  chrom  start  end` (0-based half-open, bare chrN).
"""
import argparse, random, re

ap = argparse.ArgumentParser()
ap.add_argument("--stats", required=True)
ap.add_argument("--refGenome", required=True)
ap.add_argument("--sizes", required=True)        # CSV
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
        name = p[0]
        if not name.startswith(prefix):
            continue
        bare = name[len(prefix):]
        if canon.match(bare):
            chroms.append((bare, int(p[1])))

if not chroms:
    raise SystemExit("sample_regions: no canonical %s.chrN sequences in %s" % (a.refGenome, a.stats))

sizes = [int(x) for x in a.sizes.split(",") if x]
rng = random.Random(a.seed)
for size in sizes:
    elig = [(c, L) for c, L in chroms if L >= size]
    if not elig:
        continue
    # weight by number of valid start positions -> uniform over the genome
    weights = [L - size + 1 for _, L in elig]
    for _ in range(a.nSamples):
        c, L = rng.choices(elig, weights=weights, k=1)[0]
        start = rng.randint(0, L - size)
        print("%d\t%s\t%d\t%d" % (size, c, start, start + size))
PYEOF
chmod +x "$SAMPLER"

# --- Generate the runner script (the thing sbatch executes). -----------
RUNNER="$OUTDIR/bench.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -uo pipefail
# Don't 'set -e': we want per-cell exits captured, not the job aborted.

ITUI_FILE="$ITUI_FILE"
DEPTH_BW="$DEPTH_BW"
BB="$BB"
CHROM_SRC="$CHROM_SRC"
OUTDIR="$OUTDIR"
T_TOTAL=$T_TOTAL
REF_GENOME="$REF_GENOME"
SIZES_CSV="$SIZES_CSV"
N_SAMPLES=$N_SAMPLES
SEED=$SEED
TIME_BUDGET=$TIME_BUDGET
STAGE_LOCAL=$STAGE_LOCAL
TAFFY="$TAFFY"
BIGBEDTOBED="$BIGBEDTOBED"
SIZES=( ${SIZES[*]} )

BENCH_TSV="\$OUTDIR/bench.tsv"
LOGDIR="\$OUTDIR/logs"
mkdir -p "\$LOGDIR"

# --- Stage inputs to local scratch (\$TMPDIR or /tmp fallback). -----
# The .tui is seeked at random per query and the .bw / .bb likewise; on a
# network FS that random access dominates.  PARALLEL cp + trap-cleanup so an
# aborted job doesn't leave scratch behind.  tui_load does no index-staleness
# check, so no mtime touch is needed (unlike the .tai-based view bench).
if [[ "\$STAGE_LOCAL" -eq 1 ]]; then
    SCRATCH="\${TMPDIR:-/tmp/uni_depth_summary_\${SLURM_JOB_ID:-\$\$}}"
    STAGE_DIR="\$SCRATCH/uni_depth_summary_stage_\${SLURM_JOB_ID:-\$\$}"
    mkdir -p "\$STAGE_DIR"
    trap 'rm -rf "\$STAGE_DIR" 2>/dev/null || true' EXIT
    stage_pids=()
    stage_bg() {
        local src="\$1"
        [[ -n "\$src" && -f "\$src" ]] || return 0
        local dst="\$STAGE_DIR/\$(basename "\$src")"
        echo "stage: \$src -> \$dst (\$(stat -Lc %s "\$src" 2>/dev/null || echo ?) bytes)" >&2
        ( t0=\$SECONDS; cp "\$src" "\$dst"; \\
          echo "       done: \$(basename "\$src") in \$((SECONDS - t0)) s" >&2 ) &
        stage_pids+=( \$! )
    }
    stage_bg "\$ITUI_FILE"; stage_bg "\$DEPTH_BW"; stage_bg "\$DEPTH_BW.names"; stage_bg "\$BB"
    stage_rc=0
    for p in "\${stage_pids[@]}"; do wait "\$p" || stage_rc=1; done
    [[ "\$stage_rc" -eq 0 ]] || { echo "ERROR: a stage-in cp returned non-zero" >&2; exit 1; }
    # VERIFY each staged copy matches its source byte count.  A cp to a full FS
    # can return 0 yet truncate (write succeeds, flush doesn't); the .tui-based
    # reads then fail cryptically downstream -- catch truncation here.
    for f in "\$ITUI_FILE" "\$DEPTH_BW" "\$DEPTH_BW.names" "\$BB"; do
        [[ -f "\$f" ]] || continue
        ss=\$(stat -Lc %s "\$f" 2>/dev/null); ds=\$(stat -c %s "\$STAGE_DIR/\$(basename "\$f")" 2>/dev/null || echo -1)
        [[ "\$ss" == "\$ds" ]] || { echo "ERROR: staged \$(basename "\$f") is \$ds bytes, source is \$ss -- truncated (scratch full?)" >&2; exit 1; }
    done
    ITUI_FILE="\$STAGE_DIR/\$(basename "\$ITUI_FILE")"
    DEPTH_BW="\$STAGE_DIR/\$(basename "\$DEPTH_BW")"
    BB="\$STAGE_DIR/\$(basename "\$BB")"
    echo "stage: all inputs staged + size-verified to \$STAGE_DIR" >&2
fi

# taffy lift opens <stem>.tui -- hand it the .tui path minus its extension.
ITUI_STEM="\${ITUI_FILE%.tui}"

# --- Sample the regions once (same regions for both tools). ------------
# taffy stats -s on the hg38-anchored --chromSrc gives every ref-seq name +
# length; the python sampler filters to canonical \$REF_GENOME.chrN and draws.
CHROM_STATS="\$OUTDIR/hg38.stats.txt"
REGIONS="\$OUTDIR/regions.tsv"
"\$TAFFY" stats -s -i "\$CHROM_SRC" > "\$CHROM_STATS" 2> "\$OUTDIR/stats.err" || { echo "ERROR: taffy stats -s failed on \$CHROM_SRC:" >&2; cat "\$OUTDIR/stats.err" >&2; exit 1; }
python3 "\$OUTDIR/sample_regions.py" --stats "\$CHROM_STATS" --refGenome "\$REF_GENOME" \\
    --sizes "\$SIZES_CSV" --nSamples "\$N_SAMPLES" --seed "\$SEED" > "\$REGIONS" \\
    || { echo "ERROR: region sampling failed" >&2; exit 1; }
echo "sampled \$(wc -l < "\$REGIONS") regions (\$N_SAMPLES per size, seed \$SEED) over sizes \$SIZES_CSV" >&2

# Truncate + write the header fresh (one job -> one bench.tsv).
printf "tool\tsize_bp\tchrom\tstart\tsample\twall_s\tpeak_rss_kb\tout_bytes\texit\ttimed_out\n" > "\$BENCH_TSV"

# Run one cell.  Args: tool size chrom start sample budget cmd...
# Pipe the tool's stdout through wc -c so we never write a big bedGraph/BED to
# (network) OUTDIR -- but capture the byte COUNT for the output-size panel.
# /usr/bin/time -f '%e %M' gives wall + peak RSS (KB); PIPESTATUS[0] is the
# tool's real exit (124 = timeout expiry, 137 = SIGKILL).
run_cell() {
    local tool="\$1" N="\$2" chrom="\$3" rstart="\$4" sample="\$5" budget="\$6"
    shift 6
    local stem="\${tool}_\${N}_\${sample}"
    local time_file="\$LOGDIR/time_\${stem}.txt"
    local err_file="\$LOGDIR/err_\${stem}.log"
    local bytes_file="\$LOGDIR/bytes_\${stem}.txt"

    /usr/bin/time -q -f '%e %M' -o "\$time_file" \\
        timeout --signal=KILL "\$budget" "\$@" \\
        2> "\$err_file" | wc -c > "\$bytes_file"
    local rc=\${PIPESTATUS[0]}

    local wall rss timed_out=0
    if [[ -s "\$time_file" ]]; then
        read -r wall rss < "\$time_file"
        if ! [[ "\$wall" =~ ^[0-9.]+\$ ]]; then
            read -r wall rss < <(awk '/^[0-9.]+[ \t][0-9]+\$/ {l=\$0} END{print l}' "\$time_file")
            [[ -z "\$wall" ]] && wall="NA"
        fi
    else
        wall="NA"
    fi
    [[ "\$wall" =~ ^[0-9.]+\$ ]] || wall="NA"
    [[ "\$rss"  =~ ^[0-9]+\$  ]] || rss="NA"
    local out_bytes; read -r out_bytes < "\$bytes_file" 2>/dev/null || out_bytes=0
    [[ "\$out_bytes" =~ ^[0-9]+\$ ]] || out_bytes=0
    if (( rc == 137 || rc == 124 )); then timed_out=1; fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n" \\
        "\$tool" "\$N" "\$chrom" "\$rstart" "\$sample" "\$wall" "\$rss" "\$out_bytes" "\$rc" "\$timed_out"
}

# Run one region: ours (lift --bigwig) then theirs (bigBedToBed), SEQUENTIALLY,
# so the two tools never contend within a datapoint.  Regions run concurrently
# (throttle below).  BIN = max(1000, size/1000) -> ~1000 bins/region.
run_region() {
    local N="\$1" chrom="\$2" rstart="\$3" rend="\$4" sample="\$5"
    local bin=\$(( N / 1000 )); (( bin < 1000 )) && bin=1000
    local region="\${REF_GENOME}.\${chrom}:\${rstart}-\${rend}"

    run_cell ours "\$N" "\$chrom" "\$rstart" "\$sample" "\$TIME_BUDGET" \\
        "\$TAFFY" lift -i "\$ITUI_STEM" --bigwig "\$DEPTH_BW" \\
                 -g "\$REF_GENOME" -r "\$region" -B "\$bin"

    if [[ -n "\$BIGBEDTOBED" ]]; then
        run_cell theirs "\$N" "\$chrom" "\$rstart" "\$sample" "\$TIME_BUDGET" \\
            "\$BIGBEDTOBED" -chrom="\$chrom" -start="\$rstart" -end="\$rend" "\$BB" stdout
    fi
}

# --- Concurrency throttle: at most T_TOTAL regions running at once. -----
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

# --- Per-size wave: every sampled region of this size, both tools, throttled.
for N in "\${SIZES[@]}"; do
    echo "=== size wave N=\$N ===" >&2
    t_wave=\$SECONDS
    unset pids rowfiles; declare -A pids rowfiles   # CLEAR stale keys: declare -A on an
                                                    # existing assoc array does NOT reset it.
    launched=0; unset pid_threads; declare -A pid_threads
    sample=0
    while IFS=\$'\t' read -r rsize chrom rs re; do
        [[ "\$rsize" == "\$N" ]] || continue
        sample=\$(( sample + 1 ))
        rf="\$LOGDIR/row_\${N}_\${sample}.tsv"
        rowfiles["\$sample"]="\$rf"
        acquire_slot 1
        ( run_region "\$N" "\$chrom" "\$rs" "\$re" "\$sample" ) > "\$rf" &
        pids["\$sample"]=\$!; pid_threads[\$!]=1
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

# --- Companion plot: single LINEAR panel, mean wall + min/max band, 2 lines. --
PLOT="$OUTDIR/plot.py"
cat > "$PLOT" <<'PY'
#!/usr/bin/env python3
"""Zoom-out (summary) bench: for each backend, the MEAN over the sampled regions
(shaded min/max band) vs query size, in THREE linear panels -- wall time, peak
RSS, output size -- two lines (universal-depth bigWig lift vs UCSC bigMafSummary)
to mirror the view-sample (#5) figure.

Timed-out / errored samples are excluded from the aggregate; a stderr note
reports how many were dropped per (tool,size).  If EVERY sample of a tool at a
size dropped, that point is simply absent -- the line stops.
"""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

bench_dir = os.path.dirname(os.path.abspath(__file__))
rows = []
with open(os.path.join(bench_dir, "bench.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            r["size_bp"]     = int(r["size_bp"])
            r["wall_s"]      = float(r["wall_s"]) if r["wall_s"] != "NA" else None
            r["peak_rss_kb"] = float(r["peak_rss_kb"]) if r["peak_rss_kb"] != "NA" else None
            r["out_bytes"]   = int(r["out_bytes"])
            r["timed_out"]   = int(r["timed_out"])
            r["exit"]        = int(r["exit"])
            rows.append(r)
        except (ValueError, KeyError):
            continue

label = {
    "ours":   "taffy lift --bigwig (depth)",
    "theirs": "bigMafSummary",
}
colors = {"ours": "#9467bd", "theirs": "#ff7f0e"}
order = ["ours", "theirs"]
present = {r["tool"] for r in rows}
tools = [t for t in order if t in present] + sorted(present - set(order))

def valid(r):
    return (not r["timed_out"] and r["exit"] == 0 and r["out_bytes"] > 0
            and r["wall_s"] is not None and r["peak_rss_kb"] is not None)

def agg(tool, field, scale):
    by = {}
    for r in rows:
        if r["tool"] != tool or not valid(r):
            continue
        by.setdefault(r["size_bp"], []).append(r[field] * scale)
    out = []
    for s in sorted(by):
        vs = by[s]
        out.append((s, sum(vs) / len(vs), min(vs), max(vs), len(vs)))
    return out

# Report dropped (timed-out / errored) samples per (tool,size).
drop = {}
for r in rows:
    if not valid(r):
        drop.setdefault((r["tool"], r["size_bp"]), 0)
        drop[(r["tool"], r["size_bp"])] += 1
for (t, s), n in sorted(drop.items()):
    print(f"note: {t} size={s}: {n} sample(s) dropped (timeout/error)", file=sys.stderr)

fig, axes = plt.subplots(1, 3, figsize=(21, 6))
fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.12, wspace=0.24)
panels = [("wall_s", 1.0, "wall time", "seconds"),
          ("peak_rss_kb", 1 / 1024.0, "peak RSS", "MB"),
          ("out_bytes", 1 / 1e6, "output size", "MB")]

nsamples = 0
for ax, (field, scale, title, ylab) in zip(axes, panels):
    for tool in tools:
        pts = agg(tool, field, scale)
        if not pts:
            continue
        xs   = [p[0] for p in pts]
        mean = [p[1] for p in pts]
        lo   = [p[2] for p in pts]
        hi   = [p[3] for p in pts]
        nsamples = max(nsamples, max(p[4] for p in pts))
        c = colors.get(tool)
        ax.plot(xs, mean, "o-", color=c, label=label.get(tool, tool))
        ax.fill_between(xs, lo, hi, color=c, alpha=0.15)
    # Linear axes (browser zoom-out).
    ax.set_xlabel("query size (bp)")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1e6:g}M" if x >= 1e6 else (f"{x/1e3:g}k" if x >= 1e3 else f"{x:g}")))
    ax.legend(fontsize=9)

fig.suptitle(f"universal-depth bigWig lift vs bigMafSummary: 577-way hg38 zoom-out, "
             f"random genome-wide regions, mean ± range over up to {nsamples} samples/size", fontsize=12, y=0.98)
out = os.path.join(bench_dir, "bench.png")
fig.savefig(out, dpi=140)
print(f"wrote {out}")
PY
chmod +x "$PLOT"

# --- Submit the SLURM job. --------------------------------------------
SBATCH_ARGS=(
    --cpus-per-task="$T_TOTAL"
    --mem="${SBATCH_MEM}G"
    --time="${SBATCH_TIME}:00:00"
    --output="$OUTDIR/slurm_%j.log"
    --error="$OUTDIR/slurm_%j.err.log"
    -J uni_depth_summary_bench
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
    echo ">> stdout: $OUTDIR/slurm_${JOB}.log"
    echo ">> stderr: $OUTDIR/slurm_${JOB}.err.log"
fi

if [[ "$DRY_RUN" -ne 1 && "$WAIT" -eq 1 ]]; then
    echo ">> waiting for job $JOB ..."
    while squeue -j "$JOB" -h -o "%T" 2>/dev/null | grep -qE "PENDING|RUNNING|CONFIGURING|COMPLETING|RESIZING|SUSPENDED|REQUEUED"; do
        sleep 60
    done
    FINAL_STATE=$(sacct -j "$JOB" -X -n -o State 2>/dev/null | head -1 | tr -d ' ')
    echo ">> job $JOB final state: ${FINAL_STATE:-UNKNOWN}"
    case "$FINAL_STATE" in
        COMPLETED) ;;
        *)         echo ">> NON-SUCCESS state -- check $OUTDIR/slurm_${JOB}.err.log"; exit 1;;
    esac
fi

echo ">> done."
echo ">> after the job completes:"
echo "     cat $OUTDIR/bench.tsv"
echo "     python3 $OUTDIR/plot.py   # writes bench.png"
