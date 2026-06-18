#!/bin/bash
#
# taffyBlockViz vs halBlockViz browser-query benchmark -- SLURM
#
# Mimics what a genome browser does when it draws a "snake" track: for a
# display window [tStart, tEnd) on a reference (display) chromosome, ask
# the alignment "what does target species Q look like here?".  The two
# tools answer the SAME 6-arg query:
#
#   taffyBlockVizTest <CHAINED_TUI> <qSp> hg38 <tChrom> <tStart> <tEnd>
#   blockVizBed       <HAL>         <qSp> hg38 <tChrom> <tStart> <tEnd>
#
# tSpecies = hg38 (the browser display genome / reference), qSpecies =
# the target snake genome, tChrom = an hg38 chromosome, [tStart,tEnd) =
# the query window.  Same convention for both, so a cell is a clean
# parallel pair.
#
# THE SLIDE'S WHOLE POINT
# -----------------------
# The .tui carries a chained-LOD so it answers full-chromosome queries in
# a split second.  The HAL has NO LOD file (blockVizBed gets the raw .hal),
# so it is only viable on SMALL windows -- a full-chrom hal query can burn
# hours.  We therefore CAP the hal cells at --halMaxSize bp and simply do
# NOT launch a hal cell above it (we never run-then-timeout a doomed hal
# query).  The taffy cells run the FULL size ladder including whole-chrom.
#
# Matrix: target-genome panel (-L) x size ladder (-S CSV).  Per (species,
# size) the taffy cell always runs; the hal cell runs iff size <=
# --halMaxSize.  Each cell records wall + peak-RSS + exit + timed-out +
# out_bytes to bench.tsv.
#
# Run model: one SLURM job; all cells fire in parallel under a concurrency
# throttle bounded by T_TOTAL.  Each cell is single-threaded.
#
# Usage:
#   taffy_blockviz_bench_slurm.sh \
#       -u CHAINED.tui  -H HAL.hal  -L PANEL.tsv  -S 1000,10000,...  \
#       -r chr20  --halMaxSize 1000000  -o OUTDIR  [options]

set -euo pipefail

CHAINED_TUI=""
HAL=""
PANEL=""
SIZES_CSV=""
HAL_MAX_SIZE=1000000                 # --halMaxSize: hal cells skipped above this (bp)
REF_CHROM="chr20"                    # -r: tChrom passed to BOTH tools (bare hg38 seq name)
TSPECIES="hg38"                      # tSpecies (display/reference genome) -- both tools
START=0                              # window start; each ladder cell N -> [START, START+N)
OUTDIR=""
T_TOTAL=24                           # cpus-per-task; one core per concurrent cell
MAX_OUTPUT_BLOCKS=500                # --maxOutputBlocks for both tools (taffy native; we
                                     # mirror it onto hal's doDupes path by capping via budget)
TIME_BUDGET=1800                     # per-cell wall cap (timeout sends SIGKILL)
SBATCH_TIME=24
SBATCH_MEM=64
TMP_GB=""
STAGE_LOCAL=1
PARTITION=""
ACCOUNT=""
DRY_RUN=0
WAIT=1
# Resolve from env, else PATH, else next to the `taffy` binary
# (taffyBlockVizTest is built in the same bin/).  No hardcoded absolute path.
TAFFYBLOCKVIZ="${TAFFYBLOCKVIZ:-$(command -v taffyBlockVizTest 2>/dev/null || { _t=$(command -v taffy 2>/dev/null) && echo "${_t%/*}/taffyBlockVizTest"; } || true)}"
BLOCKVIZBED="${BLOCKVIZBED:-$(command -v blockVizBed 2>/dev/null || true)}"

# Default panel, baked here so the script is one-file portable.  3 cols
# tab-separated: <genome_id>\t<scientific>\t<english>.  Override with -L.
DEFAULT_PANEL_TSV=$(cat <<'EOF'
GCA_009914755.4	Homo_sapiens	Human_T2T_CHM13
GCA_944039275.1	Danio_rerio	Zebrafish
GCF_011100685.1	Canis_lupus_familiaris	Dog
GCF_016700215.2	Gallus_gallus	Chicken
GCF_037176765.1	Anolis_sagrei	Brown_Anole
EOF
)

usage() {
    cat >&2 <<EOF
taffy_blockviz_bench_slurm.sh -- taffyBlockViz vs halBlockViz browser-query
                                 bench across a target panel x size ladder

Required:
  -u FILE       CHAINED universal .tui (e.g. *.uni.taf.gz.chained_g10000.tui).
                Passed verbatim to taffyBlockVizTest as <tuiPath>; the tool
                opens this file directly (it is NOT a MAF sibling -- pass the
                .tui itself, not a .uni.taf.gz).
  -H FILE       HAL file (raw .hal -- there is NO LOD; blockVizBed reads it
                directly, which is why hal cells must be capped, see
                --halMaxSize).
  -L FILE       Target-genome panel.  Format: 3 whitespace-separated cols per
                line: <genome_id> <scientific> <english>.  '#' = comment.
                Default: a baked 5-species panel.
  -S CSV        Query size ladder in bp (e.g. 1000,10000,100000,1000000,...).
                Required (no default -- you choose the whole-chrom cap).
  -o DIR        Output directory

Optional:
  -r NAME       tChrom passed to BOTH tools (default $REF_CHROM).  This is a
                BARE hg38 sequence name (e.g. 'chr20'), NOT 'hg38.chr20':
                both taffyBlockVizTest and blockVizBed take qSpecies /
                tSpecies separately and want the raw chrom for tChrom.
  --tSpecies N  Display / reference genome name for tSpecies (default
                $TSPECIES).  Both tools receive it as the 3rd positional arg.
  --start INT   Window start (default $START).  Each ladder cell N becomes
                the window [START, START+N).
  --halMaxSize N  Cap (bp) above which the hal (blockVizBed) cell is SKIPPED,
                not launched (default $HAL_MAX_SIZE).  taffy cells always run
                the full ladder.  This is the core of the slide: the HAL has
                no LOD, so a whole-chrom blockVizBed query would burn hours --
                we never even launch it past the cap.
  --maxOutputBlocks N  --maxOutputBlocks for taffyBlockVizTest (default
                $MAX_OUTPUT_BLOCKS).  blockVizBed has no equivalent cap flag;
                its runaway protection is --timeBudget + --halMaxSize.
  -T INT        cpus-per-task (default $T_TOTAL; one core per concurrent cell)
  --timeBudget SEC  Per-cell wall cap (timeout) (default $TIME_BUDGET)
  --time HRS    sbatch wall (default $SBATCH_TIME)
  --mem GB      sbatch mem (default $SBATCH_MEM)
  --tmp GB      Per-task local scratch requirement (sbatch --tmp=N).
                Default = (staged .tui + HAL size) + 50 GB when local-stage
                is on; unset under --no-stage-local.
  --no-stage-local
                Skip the cp of .tui + HAL to \$TMPDIR (read from network).
  --partition X --account X
  --no-wait     Submit and detach (default: driver blocks until SLURM done)
  --dry-run     Print sbatch; do not submit
  -h            Help

Override binary paths via env: TAFFYBLOCKVIZ (default
$TAFFYBLOCKVIZ), BLOCKVIZBED (default $BLOCKVIZBED).
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)                 CHAINED_TUI="$2"; shift 2;;
        -H)                 HAL="$2"; shift 2;;
        -L)                 PANEL="$2"; shift 2;;
        -S)                 SIZES_CSV="$2"; shift 2;;
        -o)                 OUTDIR="$2"; shift 2;;
        -r)                 REF_CHROM="$2"; shift 2;;
        --tSpecies)         TSPECIES="$2"; shift 2;;
        --start)            START="$2"; shift 2;;
        --halMaxSize)       HAL_MAX_SIZE="$2"; shift 2;;
        --maxOutputBlocks)  MAX_OUTPUT_BLOCKS="$2"; shift 2;;
        -T)                 T_TOTAL="$2"; shift 2;;
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

# --- Required-arg validation.  Let stderr through (never 2>/dev/null a
# helper under set -euo pipefail -- a silent error aborts the whole driver
# with no diagnostic).
for v in CHAINED_TUI HAL OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: missing required arg for \$$v" >&2; usage 1; }
done
[[ -n "$SIZES_CSV" ]] || { echo "ERROR: -S (size ladder CSV) required" >&2; usage 1; }
[[ -n "$TAFFYBLOCKVIZ" && -x "$TAFFYBLOCKVIZ" ]] || { echo "ERROR: taffyBlockVizTest not executable: $TAFFYBLOCKVIZ (set \$TAFFYBLOCKVIZ)" >&2; exit 1; }
[[ -n "$BLOCKVIZBED"  && -x "$BLOCKVIZBED"  ]] || { echo "ERROR: blockVizBed not executable: $BLOCKVIZBED (set \$BLOCKVIZBED)" >&2; exit 1; }
[[ -f "$CHAINED_TUI" ]] || { echo "ERROR: chained .tui not found: $CHAINED_TUI" >&2; exit 1; }
[[ -f "$HAL"         ]] || { echo "ERROR: HAL not found: $HAL" >&2; exit 1; }
[[ "$HAL_MAX_SIZE" =~ ^[0-9]+$ ]] || { echo "ERROR: --halMaxSize must be a non-negative integer (got '$HAL_MAX_SIZE')" >&2; exit 1; }
[[ "$START" =~ ^[0-9]+$ ]] || { echo "ERROR: --start must be a non-negative integer (got '$START')" >&2; exit 1; }

mkdir -p "$OUTDIR" "$OUTDIR/logs"
echo ">> driver starting (output dir: $OUTDIR)" >&2

# --- Resolve panel ------------------------------------------------------
# Whitespace-tolerant parse (default awk FS) so tab- AND space-separated
# user-edited panel files both work; strip '#' comments + short lines.
PANEL_TSV="$OUTDIR/species.tsv"
if [[ -n "$PANEL" ]]; then
    [[ -f "$PANEL" ]] || { echo "ERROR: panel file not found: $PANEL" >&2; exit 1; }
    awk '!/^#/ && NF >= 3 {print $1"\t"$2"\t"$3}' "$PANEL" > "$PANEL_TSV"
else
    printf "%s\n" "$DEFAULT_PANEL_TSV" > "$PANEL_TSV"
fi
N_SPECIES=$(wc -l < "$PANEL_TSV")
[[ "$N_SPECIES" -gt 0 ]] || { echo "ERROR: empty species panel" >&2; exit 1; }
echo ">> species panel: $N_SPECIES entries" >&2

# --- Parse + sort the size ladder (ascending integers). -----------------
IFS=',' read -r -a SIZE_ARR <<< "$SIZES_CSV"
for S in "${SIZE_ARR[@]}"; do
    [[ "$S" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -S sizes must be positive integers (got '$S')" >&2; exit 1; }
done
# Sort numerically + de-dup so the ladder + the plot read cleanly.
readarray -t SIZE_ARR < <(printf '%s\n' "${SIZE_ARR[@]}" | sort -n -u)
N_SIZES=${#SIZE_ARR[@]}

# Count hal cells (sizes <= cap) up front -- pure reporting / sizing.
N_HAL_SIZES=0
for S in "${SIZE_ARR[@]}"; do (( S <= HAL_MAX_SIZE )) && N_HAL_SIZES=$((N_HAL_SIZES + 1)); done

echo ">> output dir:    $OUTDIR"
echo ">> taffyBlockViz: $TAFFYBLOCKVIZ"
echo ">> blockVizBed:   $BLOCKVIZBED"
echo ">> chained .tui:  $CHAINED_TUI"
echo ">> hal:           $HAL"
echo ">> tSpecies:      $TSPECIES  (display/reference genome)"
echo ">> ref chrom:     $REF_CHROM  (bare tChrom for both tools)"
echo ">> start:         $START bp  (window = [start, start+N))"
echo ">> sizes:         ${SIZE_ARR[*]}"
echo ">> halMaxSize:    $HAL_MAX_SIZE bp  (hal cells run for $N_HAL_SIZES of $N_SIZES sizes; larger SKIPPED)"
echo ">> maxOutBlocks:  $MAX_OUTPUT_BLOCKS  (taffyBlockVizTest --maxOutputBlocks)"
echo ">> species:       $N_SPECIES"
awk -F'\t' '{printf("                  %-18s %-32s %s\n", $1, $2, $3)}' "$PANEL_TSV"
# Cells: taffy = N_SPECIES * N_SIZES; hal = N_SPECIES * N_HAL_SIZES.
TOTAL_CELLS=$(( N_SPECIES * (N_SIZES + N_HAL_SIZES) ))
echo ">> cpus/task:     $T_TOTAL  (one core per concurrent cell; ~$TOTAL_CELLS cells total)"
(( T_TOTAL >= 1 )) || { echo "ERROR: -T must be >= 1" >&2; exit 1; }
echo ">> time budget:   $TIME_BUDGET s per cell"
echo ">> local-stage:   $([[ $STAGE_LOCAL -eq 1 ]] && echo "ON (copies to \$TMPDIR)" || echo "OFF (reads from network)")"

# --- Per-job --tmp sizing hint (when local-stage is on). ---------------
# Roll up the bytes the runner actually stages: the chained .tui + the HAL.
# Both blockViz tools open their input file directly, so both are copied to
# local scratch -- request scratch sized to them.
if [[ "$STAGE_LOCAL" -eq 1 ]]; then
    STAGE_BYTES=0
    for f in "$CHAINED_TUI" "$HAL"; do
        if [[ -e "$f" ]]; then
            STAGE_BYTES=$(( STAGE_BYTES + $(stat -Lc %s "$f" 2>/dev/null || stat -f %z "$f" 2>/dev/null || echo 0) ))
        fi
    done
    STAGE_GB=$(( STAGE_BYTES / (1024**3) ))
    echo ">> stage-in size: ~${STAGE_GB} GB total"
    # Promote the hint to the actual --tmp default: scratch sized to what we
    # stage (+50 GB headroom).  An explicit --tmp still overrides.
    TMP_GB=${TMP_GB:-$(( STAGE_GB + 50 ))}
    echo ">> --tmp default: ${TMP_GB} GB per task (stage ~${STAGE_GB} GB + 50 GB headroom; override with --tmp)"
elif [[ -n "$TMP_GB" ]]; then
    echo ">> --tmp request: ${TMP_GB} GB per task"
fi

# --- Generate the runner script (the thing sbatch executes). -----------
RUNNER="$OUTDIR/bench.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -uo pipefail
# Don't 'set -e': we want per-cell exits captured, not the job aborted.

CHAINED_TUI="$CHAINED_TUI"
HAL="$HAL"
OUTDIR="$OUTDIR"
TSPECIES="$TSPECIES"
REF_CHROM="$REF_CHROM"
START=$START
HAL_MAX_SIZE=$HAL_MAX_SIZE
MAX_OUTPUT_BLOCKS=$MAX_OUTPUT_BLOCKS
TIME_BUDGET=$TIME_BUDGET
STAGE_LOCAL=$STAGE_LOCAL
T_TOTAL=$T_TOTAL
TAFFYBLOCKVIZ="$TAFFYBLOCKVIZ"
BLOCKVIZBED="$BLOCKVIZBED"
SPECIES_TSV="$OUTDIR/species.tsv"
SIZES=( ${SIZE_ARR[*]} )

BENCH_TSV="\$OUTDIR/bench.tsv"
LOGDIR="\$OUTDIR/logs"
mkdir -p "\$LOGDIR"

# --- Stage inputs to local scratch (\$TMPDIR or /tmp fallback). -----
# .tui + HAL are both large so this is the dominant pre-bench cost;
# staging once amortises it across the N_SPECIES * sizes cells that
# follow.  Both blockViz tools open their input file directly (unlike
# taffy lift, which only needs the .tui sibling) -- so we stage the
# files we actually pass.  Trap-cleanup so an aborted job doesn't leak.
if [[ "\$STAGE_LOCAL" -eq 1 ]]; then
    SCRATCH="\${TMPDIR:-/tmp/taffy_blockviz_bench_\${SLURM_JOB_ID:-\$\$}}"
    STAGE_DIR="\$SCRATCH/taffy_blockviz_stage_\${SLURM_JOB_ID:-\$\$}"
    mkdir -p "\$STAGE_DIR"
    trap 'rm -rf "\$STAGE_DIR" 2>/dev/null || true' EXIT
    stage_one() {
        local src="\$1"
        local dst="\$STAGE_DIR/\$(basename "\$src")"
        echo "stage: \$src -> \$dst (\$(stat -Lc %s "\$src" 2>/dev/null || echo ?) bytes)" >&2
        local t0=\$SECONDS
        cp "\$src" "\$dst"
        echo "       done in \$((SECONDS - t0)) s" >&2
        echo "\$dst"
    }
    CHAINED_TUI=\$(stage_one "\$CHAINED_TUI")
    HAL=\$(stage_one "\$HAL")
    echo "stage: all inputs staged to \$STAGE_DIR" >&2
fi

# Write header if file is empty / fresh.
if [[ ! -s "\$BENCH_TSV" ]]; then
    printf "tool\tgenome_id\tsci_name\tcommon_name\tsize_bp\twall_s\tpeak_rss_kb\texit\ttimed_out\tout_bytes\n" > "\$BENCH_TSV"
fi

# run_cell tool species_id sci common size cmd...
# Writes a single TSV row to stdout.  out_bytes = size of the tool's
# stdout (the BED blocks it emitted); we keep the err log for debugging
# and discard the output dump after measuring.
run_cell() {
    local tool="\$1" sid="\$2" sci="\$3" common="\$4" N="\$5"
    shift 5
    local stem="\${tool}_\${sid}_\${N}"
    local time_file="\$LOGDIR/time_\${stem}.txt"
    local out_file="\$LOGDIR/out_\${stem}"     # discarded after wc
    local err_file="\$LOGDIR/err_\${stem}.log"
    # Per-cell scratch CWD.  blockVizBed writes <tSp>.bed + <qSp>.bed into the
    # CWD with HARDCODED names, so without an isolated dir concurrent hal cells
    # would all clobber a shared hg38.bed.  taffyBlockVizTest writes only to
    # stdout, so its scratch dir just stays empty.  The stdout/err/time fds are
    # opened by THIS shell BEFORE the cd, so they're unaffected by it.
    local cell_dir="\$LOGDIR/cwd_\${stem}"
    rm -rf "\$cell_dir"; mkdir -p "\$cell_dir"

    # -q suppresses /usr/bin/time's "Command exited with non-zero status N"
    # line, which would otherwise occupy the time_file's first row and
    # break our read.  The bash -c cd-wrapper runs the tool in cell_dir;
    # 'exec' so timeout's SIGKILL lands on the real tool, not the wrapper.
    /usr/bin/time -q -f '%e %M' -o "\$time_file" \\
        timeout --signal=KILL "\$TIME_BUDGET" \\
        bash -c 'cd "\$0" && exec "\$@"' "\$cell_dir" "\$@" \\
        > "\$out_file" 2> "\$err_file"
    local rc=\$?

    local wall rss out_bytes timed_out=0
    if [[ -s "\$time_file" ]]; then
        read -r wall rss < "\$time_file"
        # Belt-and-suspenders: if the first line isn't "<float> <int>",
        # scan for the last matching line instead.
        if ! [[ "\$wall" =~ ^[0-9.]+\$ && "\$rss" =~ ^[0-9]+\$ ]]; then
            read -r wall rss < <(awk '/^[0-9.]+[ \t][0-9]+\$/ {l=\$0} END{print l}' "\$time_file")
            [[ -z "\$wall" ]] && { wall="NA"; rss="NA"; }
        fi
    else
        wall="NA"; rss="NA"
    fi
    (( rc == 137 || rc == 124 )) && timed_out=1
    out_bytes=\$(stat -c %s "\$out_file" 2>/dev/null || echo 0)
    # Fold in any files the tool dropped in its scratch CWD (blockVizBed's
    # .bed output; taffyBlockVizTest leaves none) so out_bytes reflects real
    # output, not blockVizBed's ~70-byte stdout status line.  find always
    # exits 0 (empty dir -> awk prints 0), so this is pipefail/set -e safe.
    out_bytes=\$(( out_bytes + \$(find "\$cell_dir" -type f -printf '%s\\n' 2>/dev/null | awk '{s+=\$1} END{print s+0}') ))
    rm -rf "\$cell_dir"
    rm -f "\$out_file"

    printf "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%s\n" \\
        "\$tool" "\$sid" "\$sci" "\$common" "\$N" "\$wall" "\$rss" "\$rc" "\$timed_out" "\$out_bytes"
}

# --- Concurrency throttle ----------------------------------------------
# Bound the number of running cells by T_TOTAL so we never oversubscribe
# the SLURM allocation -- every cell here is single-threaded, so the
# budget is a plain count.  Without it the N_SPECIES * (sizes + hal-sizes)
# cells all launch at once and fight for cores.
THREAD_BUDGET=\$T_TOTAL
launched=0
declare -A pid_threads

acquire_slot() {
    local threads=\$1
    while (( launched + threads > THREAD_BUDGET )); do
        wait -n 2>/dev/null || true   # block until any child finishes
        local p
        for p in "\${!pid_threads[@]}"; do
            if ! kill -0 \$p 2>/dev/null; then
                launched=\$(( launched - pid_threads[\$p] ))
                unset pid_threads[\$p]
            fi
        done
    done
    launched=\$(( launched + threads ))
}
register_pid() {
    pid_threads[\$1]=\$2
}

# --- Fire all cells in parallel.  No size waves -- the throttle keeps
# concurrency bounded, and we want hal + taffy of the SAME (species,size)
# to overlap so they see the same FS load. -----------------------------
echo "=== launching cells: \$(wc -l < "\$SPECIES_TSV") species x \${#SIZES[@]} sizes (taffy) + capped hal ==="
t0=\$SECONDS
declare -A pids rowfiles

while IFS=\$'\t' read -r sid sci common; do
    [[ -n "\$sid" ]] || continue
    for N in "\${SIZES[@]}"; do
        END=\$((START + N))

        # ---- taffyBlockViz cell (FULL ladder, always runs) ----
        # taffyBlockVizTest <tuiPath> <qSp> <tSp> <tChrom> <tStart> <tEnd>
        stem_tv="taffyBlockViz_\${sid}_\${N}"
        rowfiles[\$stem_tv]="\$LOGDIR/row_\${stem_tv}.tsv"
        acquire_slot 1
        ( run_cell taffyBlockViz "\$sid" "\$sci" "\$common" "\$N" \\
            "\$TAFFYBLOCKVIZ" --maxOutputBlocks "\$MAX_OUTPUT_BLOCKS" \\
                "\$CHAINED_TUI" "\$sid" "\$TSPECIES" "\$REF_CHROM" "\$START" "\$END" \\
          ) > "\${rowfiles[\$stem_tv]}" &
        pids[\$stem_tv]=\$!; register_pid \$! 1

        # ---- blockVizBed cell (CAPPED: only when N <= halMaxSize) ----
        # We never launch a doomed whole-chrom hal query: the raw .hal has
        # no LOD, so above the cap it would burn hours.  Skip = no cell.
        # blockVizBed <halLodPath> <qSp> <tSp> <tChrom> <tStart> <tEnd>
        #             [doSeq=0] [doDupes=0] [udcPath].  We pass the raw
        # .hal as halLodPath (there is no LOD) + defaults for the rest.
        if (( N <= HAL_MAX_SIZE )); then
            stem_bv="halBlockViz_\${sid}_\${N}"
            rowfiles[\$stem_bv]="\$LOGDIR/row_\${stem_bv}.tsv"
            acquire_slot 1
            ( run_cell halBlockViz "\$sid" "\$sci" "\$common" "\$N" \\
                "\$BLOCKVIZBED" "\$HAL" "\$sid" "\$TSPECIES" "\$REF_CHROM" "\$START" "\$END" \\
              ) > "\${rowfiles[\$stem_bv]}" &
            pids[\$stem_bv]=\$!; register_pid \$! 1
        fi
    done
done < "\$SPECIES_TSV"

# Wait for all cells, append rows.
for stem in "\${!pids[@]}"; do
    wait "\${pids[\$stem]}" || true
    [[ -s "\${rowfiles[\$stem]}" ]] && cat "\${rowfiles[\$stem]}" >> "\$BENCH_TSV"
done

echo "=== all cells done in \$((SECONDS - t0)) s ==="
echo "bench done.  TSV: \$BENCH_TSV"
EOF
chmod +x "$RUNNER"

# --- Companion plot script.  Best-effort / optional: a missing
# matplotlib must NOT fail the job, so the wrapper guards the import. ----
PLOT="$OUTDIR/plot.py"
cat > "$PLOT" <<'PY'
#!/usr/bin/env python3
"""Plot wall + peak-RSS vs query size for the blockViz bench (log-x).

One line per tool; taffyBlockViz runs the full ladder, halBlockViz stops
at the --halMaxSize cap.  We draw a vertical dashed line at the largest
size halBlockViz actually ran (the visible "hal stops here" marker) and
annotate it.

Timed-out cells get a hollow 'X' at the budget value -- "killed at
budget" rather than silently dropped.  Best-effort: a missing matplotlib
prints a hint and exits 0 (never fails the SLURM job)."""
import csv, os, sys
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as e:
    print(f"plot skipped (matplotlib unavailable: {e})")
    sys.exit(0)

bench_dir = os.path.dirname(os.path.abspath(__file__))
rows = []
with open(os.path.join(bench_dir, "bench.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            r["size_bp"]     = int(r["size_bp"])
            r["wall_s"]      = float(r["wall_s"]) if r["wall_s"] != "NA" else None
            r["peak_rss_kb"] = float(r["peak_rss_kb"]) if r["peak_rss_kb"] != "NA" else None
            r["timed_out"]   = int(r["timed_out"])
            r["exit"]        = int(r["exit"])
            r["out_bytes"]   = int(r["out_bytes"])
            rows.append(r)
        except (ValueError, KeyError):
            continue

# One series per tool, aggregating across the species panel: at each size
# we plot the MEAN wall / RSS over species (a representative browser-query
# cost), with the min/max as a light band.
colors = {"taffyBlockViz": "#1f77b4", "halBlockViz": "#d62728"}
order  = ["taffyBlockViz", "halBlockViz"]
present = {r["tool"] for r in rows}
tools = [t for t in order if t in present] + sorted(present - set(order))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.subplots_adjust(left=0.06, right=0.97, top=0.90, bottom=0.12, wspace=0.22)

def agg(tool, field):
    """size -> (mean, lo, hi) over species for completed (non-timed-out,
    exit==0) cells."""
    by_size = {}
    for r in rows:
        if r["tool"] != tool or r["timed_out"] or r["exit"] != 0:
            continue
        v = r[field]
        if v is None:
            continue
        by_size.setdefault(r["size_bp"], []).append(v)
    out = []
    for s in sorted(by_size):
        vs = by_size[s]
        out.append((s, sum(vs) / len(vs), min(vs), max(vs)))
    return out

hal_max_run = 0
for r in rows:
    if r["tool"] == "halBlockViz" and not r["timed_out"] and r["exit"] == 0:
        hal_max_run = max(hal_max_run, r["size_bp"])

for tool in tools:
    color = colors.get(tool)
    for ax, field, scale in ((ax1, "wall_s", 1.0), (ax2, "peak_rss_kb", 1 / 1024.0)):
        pts = agg(tool, field)
        if not pts:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] * scale for p in pts]
        lo = [p[2] * scale for p in pts]
        hi = [p[3] * scale for p in pts]
        ax.plot(xs, ys, "o-", label=tool, color=color)
        ax.fill_between(xs, lo, hi, color=color, alpha=0.15)
    # Timed-out markers (wall panel only): hollow X at the budget.
    to = [(r["size_bp"], r["wall_s"]) for r in rows
          if r["tool"] == tool and r["timed_out"] and r["wall_s"] is not None]
    if to:
        ax1.plot([p[0] for p in to], [p[1] for p in to], "X", color=color,
                 markerfacecolor="none", markeredgewidth=2, markersize=11,
                 label=f"{tool} (timed out)")

# "hal stops here" marker: vertical line at the largest size halBlockViz ran.
if hal_max_run > 0:
    for ax in (ax1, ax2):
        ax.axvline(hal_max_run, color="#d62728", linestyle=":", alpha=0.7)
        ymin, ymax = ax.get_ylim()
        ax.text(hal_max_run, ymax, " halBlockViz cap\n (no LOD)",
                color="#d62728", fontsize=8, va="top", ha="left")

for ax, title, ylab in [(ax1, "browser-query wall time (mean over panel)", "seconds"),
                        (ax2, "peak RSS (mean over panel)", "MB")]:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("query window size (bp)")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

out = os.path.join(bench_dir, "bench.png")
fig.savefig(out, dpi=140)
print(f"wrote {out}")
PY
chmod +x "$PLOT"

# --- Submit. -----------------------------------------------------------
SBATCH_ARGS=(
    --cpus-per-task="$T_TOTAL"
    --mem="${SBATCH_MEM}G"
    --time="${SBATCH_TIME}:00:00"
    --output="$OUTDIR/slurm_%j.log"
    --error="$OUTDIR/slurm_%j.err.log"
    -J taffy_blockviz_bench
)
[[ -n "$TMP_GB"    ]] && SBATCH_ARGS+=( --tmp="${TMP_GB}G" )
[[ -n "$PARTITION" ]] && SBATCH_ARGS+=( --partition="$PARTITION" )
[[ -n "$ACCOUNT"   ]] && SBATCH_ARGS+=( --account="$ACCOUNT" )

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo ">> DRY RUN -- would submit:"
    echo "sbatch ${SBATCH_ARGS[*]} --parsable $RUNNER"
    exit 0
fi

echo ">> submitting..."
JOB=$(sbatch "${SBATCH_ARGS[@]}" --parsable "$RUNNER")
echo ">> job id: $JOB"
echo ">> stdout: $OUTDIR/slurm_${JOB}.log"
echo ">> stderr: $OUTDIR/slurm_${JOB}.err.log"

if [[ "$WAIT" -eq 1 ]]; then
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
