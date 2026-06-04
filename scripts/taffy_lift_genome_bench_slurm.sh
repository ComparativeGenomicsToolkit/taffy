#!/bin/bash
#
# taffy lift / halLiftover whole-genome benchmark -- SLURM
#
# Lifts the WHOLE reference genome (one bed line per chromosome covering
# 0..chrom_length) from REF (default: hg38) to N target species spread
# across the tree.  No size-sweep, no chains, no liftover -- this bench is
# specifically about whole-genome lifts where chain-based tools fail
# trivially and the contest is taffy lift vs halLiftover.
#
# Run model: one SLURM job; all (n_species x 2 tools) cells fire in
# parallel.  Cells are independent (different output paths, no shared
# writes), and each cell is single-threaded.
#
# Default 9-species panel (override with -L FILE; same format as the
# chicken lift-bench).  Spread chosen for x-axis coverage from human:
#   GCF_028858775.2   Pan_troglodytes              (0.007)  primate
#   GCF_011100685.1   Canis_lupus_familiaris       (0.251)  carnivore
#   GCA_000001635.9   mm39                         (0.328)  rodent
#   GCF_027887165.1   Monodelphis_domestica        (0.362)  marsupial
#   GCF_004115215.2   Ornithorhynchus_anatinus     (0.468)  monotreme
#   GCF_015237465.2   Chelonia_mydas               (0.588)  reptile
#   GCF_016700215.2   Gallus_gallus                (0.746)  bird
#   GCF_905171765.1   Bufo_bufo                    (0.958)  amphibian
#   GCA_944039275.1   Danio_rerio                  (1.449)  fish
#
# Usage:
#   taffy_lift_genome_bench_slurm.sh \
#       -u UNI.maf.gz  -H HAL.hal  -t TREE.nwk \
#       -o OUTDIR  [options]

set -euo pipefail

UNI=""
UNI_TAF=""
HAL=""
TREE=""
FAST_MODE=0                          # --fast: add chunk-walk variant cells alongside default
BIN_SIZES_CSV=""                     # --bin CSV: add bedGraph variants per bin size; empty = off.
THREADS_PER_CELL_CSV="1"             # --threadsPerCell CSV: OMP_NUM_THREADS values to bench; tool gets _t<N> suffix when N>1.
OUTDIR=""
REF="GCA_000001405.15"              # hg38
T_TOTAL=32                          # cpus-per-task; ≥ 2 × n_species
TIME_BUDGET=14400                   # per-cell wall cap (4 h)
SBATCH_TIME=24
SBATCH_MEM=128
TMP_GB=""
STAGE_LOCAL=1
SPECIES_FILE=""
PARTITION=""
ACCOUNT=""
DRY_RUN=0
WAIT=1
TAFFY="${TAFFY:-$(command -v taffy || true)}"
HALLIFTOVER="${HALLIFTOVER:-$(command -v halLiftover || true)}"

DEFAULT_SPECIES_TSV=$(cat <<'EOF'
GCF_028858775.2	Pan_troglodytes	Chimpanzee
GCF_011100685.1	Canis_lupus_familiaris	Dog
GCA_000001635.9	mm39	House_Mouse
GCF_027887165.1	Monodelphis_domestica	Opossum
GCF_004115215.2	Ornithorhynchus_anatinus	Platypus
GCF_015237465.2	Chelonia_mydas	Green_Sea_Turtle
GCF_016700215.2	Gallus_gallus	Chicken
GCF_905171765.1	Bufo_bufo	Common_Toad
GCA_944039275.1	Danio_rerio	Zebrafish
EOF
)

usage() {
    cat >&2 <<EOF
taffy_lift_genome_bench_slurm.sh -- whole-genome lift bench, taffy vs halLiftover

Required (at least one of -u / --uniTaf must be set):
  -u FILE       Universal MAF (.uni.maf.gz; only its .tui sibling is read).
                If set, adds a "maf.tui" cell per species (taffy lift
                against the MAF-anchored .tui).
  --uniTaf FILE Universal TAF (.uni.taf.gz; only its .tui sibling is read).
                If set, adds a "taf.tui" cell per species.  Both flags
                can be set together to bench both .tui formats at once.
  -H FILE       HAL file (for halLiftover)
  -t FILE       Species tree (.nwk) -- used by the plot script for divergence
  -o DIR        Output directory

Optional:
  -r ID         Reference genome ID (default $REF; must match the genome
                prefix in the universal MAF / HAL)
  -T INT        cpus-per-task (default $T_TOTAL).  Needs ≥ N_SPECIES *
                (1 + N_SOURCES * N_MODES) for full parallelism, where
                N_SOURCES = (-u set) + (--uniTaf set) and N_MODES = 2 if
                --fast else 1.  E.g. 9 species + both sources + --fast =
                45 cells; lower T_TOTAL queues them.
  -L FILE       Override species panel.  Format: 3 whitespace-separated
                cols per line: <genome_id> <scientific> <english>.
                '#' = comment.
  --fast        For each taffy cell, also launch a --fast variant using
                the chunk-iteration lift path (10-50x faster on whole-
                chromosome / whole-genome queries).  Tool name in
                bench.tsv gets a '_fast' suffix: 'maf.tui_fast' /
                'taf.tui_fast'.  Default cells are still emitted so the
                bench can compare default-vs-fast wall + verify parity.
  --bin CSV     taffy lift --bin sizes (bedGraph, coarse-grained browser
                tracks) to bench side-by-side.  Each value adds a
                --fast --bin <N> variant cell per species per source.
                Auto-enables --fast.  Tool name suffix '_bin<S>' (human-
                readable when divisible: 1M, 100k, 1G; raw int otherwise).
                Example: --bin 100000,1000000 adds 'maf.tui_fast_bin100k'
                and 'maf.tui_fast_bin1M' per species.
  --threadsPerCell CSV
                OMP_NUM_THREADS values to bench (default "1").  Each
                value becomes a variant cell per (species, source, mode,
                bin) -- A/B compare parallel chunk_decode wall.  Tool
                name suffix '_t<N>' when N>1; N=1 keeps the bare name.
                halLiftover ignores OMP_NUM_THREADS so doesn't get _t<N>
                variants.  T_TOTAL should be >= concurrent-cells *
                max-N to avoid oversubscription; driver warns at submit.
  --timeBudget SEC  Per-cell wall cap (timeout) (default $TIME_BUDGET)
  --time HRS    sbatch wall (default $SBATCH_TIME)
  --mem GB      sbatch mem (default $SBATCH_MEM)
  --tmp GB      Per-task local scratch requirement (sbatch --tmp=N).
                Default unset.  Size to roughly (.tui + HAL + 10%) --
                ~1100 GB for a 92 GB .tui + 965 GB HAL.
  --no-stage-local
                Skip the cp of .tui + HAL to \$TMPDIR.  Cells read from
                the network paths instead.
  --partition X --account X
  --no-wait     Submit and detach (default: driver blocks until SLURM done)
  --dry-run     Print sbatch; do not submit
  -h            Help

Override binary paths via env: TAFFY, HALLIFTOVER
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)             UNI="$2"; shift 2;;
        --uniTaf)       UNI_TAF="$2"; shift 2;;
        --fast)         FAST_MODE=1; shift;;
        --bin)          BIN_SIZES_CSV="$2"; shift 2;;
        --threadsPerCell) THREADS_PER_CELL_CSV="$2"; shift 2;;
        -H)             HAL="$2"; shift 2;;
        -t)             TREE="$2"; shift 2;;
        -o)             OUTDIR="$2"; shift 2;;
        -r)             REF="$2"; shift 2;;
        -T)             T_TOTAL="$2"; shift 2;;
        -L)             SPECIES_FILE="$2"; shift 2;;
        --timeBudget)   TIME_BUDGET="$2"; shift 2;;
        --time)         SBATCH_TIME="$2"; shift 2;;
        --mem)          SBATCH_MEM="$2"; shift 2;;
        --tmp)          TMP_GB="$2"; shift 2;;
        --no-stage-local)   STAGE_LOCAL=0; shift;;
        --partition)    PARTITION="$2"; shift 2;;
        --account)      ACCOUNT="$2"; shift 2;;
        --no-wait)      WAIT=0; shift;;
        --dry-run)      DRY_RUN=1; shift;;
        -h|--help)      usage 0;;
        *)              echo "unknown arg: $1" >&2; usage 1;;
    esac
done

for v in HAL TREE OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: -$(echo $v | cut -c1) required" >&2; usage 1; }
done
[[ -n "$UNI" || -n "$UNI_TAF" ]] || {
    echo "ERROR: at least one of -u / --uniTaf must be set" >&2; usage 1;
}
[[ -n "$TAFFY"       ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$HALLIFTOVER" ]] || { echo "ERROR: halLiftover not on PATH (set \$HALLIFTOVER)" >&2; exit 1; }
if [[ -n "$UNI" ]]; then
    [[ -f "$UNI"        ]] || { echo "ERROR: $UNI not found" >&2; exit 1; }
    [[ -f "${UNI}.tui"  ]] || { echo "ERROR: $UNI has no .tui sibling" >&2; exit 1; }
fi
if [[ -n "$UNI_TAF" ]]; then
    [[ -f "$UNI_TAF"        ]] || { echo "ERROR: $UNI_TAF not found" >&2; exit 1; }
    [[ -f "${UNI_TAF}.tui"  ]] || { echo "ERROR: $UNI_TAF has no .tui sibling" >&2; exit 1; }
fi
[[ -f "$HAL"         ]] || { echo "ERROR: $HAL not found" >&2; exit 1; }
[[ -f "$TREE"        ]] || { echo "ERROR: $TREE not found" >&2; exit 1; }

mkdir -p "$OUTDIR" "$OUTDIR/logs" "$OUTDIR/beds"
echo ">> driver starting (output dir: $OUTDIR)" >&2

# --- Resolve species panel ----------------------------------------------
SPECIES_TSV="$OUTDIR/species.tsv"
if [[ -n "$SPECIES_FILE" ]]; then
    awk '!/^#/ && NF >= 3 {print $1"\t"$2"\t"$3}' "$SPECIES_FILE" > "$SPECIES_TSV"
else
    printf "%s\n" "$DEFAULT_SPECIES_TSV" > "$SPECIES_TSV"
fi
N_SPECIES=$(wc -l < "$SPECIES_TSV")
[[ "$N_SPECIES" -gt 0 ]] || { echo "ERROR: empty species panel" >&2; exit 1; }
echo ">> species panel: $N_SPECIES entries" >&2

# --- Pull REF chrom sizes from the .tui --------------------------------
# Both taffy and halLiftover get the SAME bed.  The .tui lists every
# sequence in the universal MAF; filter to REF.* and emit one bed line
# per chrom covering 0..seqLen.
# Pull chrom list from whichever .tui is set (both have the same chroms,
# they're built from the same alignment).
UNI_STATS_SRC="${UNI:-$UNI_TAF}"
echo ">> querying .tui ($UNI_STATS_SRC) for ${REF}.* chroms via taffy stats -s ..." >&2
NATIVE_BED="$OUTDIR/beds/genome.native.bed"
PREFIXED_BED="$OUTDIR/beds/genome.prefixed.bed"
"$TAFFY" stats -i "$UNI_STATS_SRC" -s \
    | awk -v p="${REF}." 'index($1, p) == 1 {
        sub(p, "", $1);
        printf "%s\t0\t%d\n", $1, $2;
      }' \
    | sort -k1,1 > "$NATIVE_BED"
[[ -s "$NATIVE_BED" ]] || {
    echo "ERROR: no sequences matching ${REF}.* in $UNI_STATS_SRC's .tui." >&2
    exit 1
}
# Prefixed bed for taffy lift (its --bed expects "<genome>.<chrom>"):
awk -v p="${REF}." '{print p$0}' "$NATIVE_BED" > "$PREFIXED_BED"

N_CHROMS=$(wc -l < "$NATIVE_BED")
REF_BP=$(awk '{s += $3 - $2} END {print s}' "$NATIVE_BED")
echo ">>   ${REF} input bed: $N_CHROMS chroms, $REF_BP bp total" >&2

# --- Resolve --bin CSV.  Default empty -> no bin cells.  Non-empty
# implies --fast (auto-enable) -- the taffy CLI rejects --bin without
# it, and we always want a plain-fast cell alongside for comparison.
if [[ -z "$BIN_SIZES_CSV" ]]; then
    BIN_SIZES_ARR=( )
else
    IFS=',' read -r -a BIN_SIZES_ARR <<< "$BIN_SIZES_CSV"
    for B in "${BIN_SIZES_ARR[@]}"; do
        [[ "$B" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --bin values must be positive integers (got '$B')" >&2; exit 1; }
    done
    FAST_MODE=1
fi
N_BIN=${#BIN_SIZES_ARR[@]}

# --- Resolve --threadsPerCell CSV.  Default "1" -> one variant with
# OMP_NUM_THREADS=1 (preserves prior tool names).  Each value fans out
# a separate taffy cell with OMP_NUM_THREADS=N; suffix '_t<N>' when N>1.
IFS=',' read -r -a THREADS_PER_CELL_ARR <<< "$THREADS_PER_CELL_CSV"
for T in "${THREADS_PER_CELL_ARR[@]}"; do
    [[ "$T" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threadsPerCell values must be positive integers (got '$T')" >&2; exit 1; }
done
N_TPC=${#THREADS_PER_CELL_ARR[@]}

echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> halLiftover:   $HALLIFTOVER"
echo ">> uni:           $UNI  (using \$UNI.tui only)"
echo ">> hal:           $HAL"
echo ">> tree:          $TREE"
echo ">> reference:     $REF"
echo ">> species:       $N_SPECIES"
awk -F'\t' '{printf("                  %-18s %-32s %s\n", $1, $2, $3)}' "$SPECIES_TSV"
# Cells per species: 1 halLiftover + (one taffy cell per .tui source) x
# (1 mode if --fast off, 2 modes if --fast on).  This drives the T_TOTAL
# sizing recommendation.
N_SOURCES=0
[[ -n "$UNI"     ]] && N_SOURCES=$((N_SOURCES + 1))
[[ -n "$UNI_TAF" ]] && N_SOURCES=$((N_SOURCES + 1))
N_MODES=$([[ "$FAST_MODE" -eq 1 ]] && echo 2 || echo 1)
# Per source: N_MODES standard cells (default + fast) + N_BIN bin cells
# (always fast-based, only added when N_BIN > 0).  Each taffy cell is
# itself replicated N_TPC times (one per threadsPerCell value).
CELLS_PER_SPECIES=$(( 1 + N_SOURCES * (N_MODES + N_BIN) * N_TPC ))
TOTAL_CELLS=$(( N_SPECIES * CELLS_PER_SPECIES ))
# Thread demand: sum of all OMP_NUM_THREADS values across concurrent
# taffy cells (halLiftover gets 1).  Used for the oversubscription warn.
_TPC_SUM=0
for T in "${THREADS_PER_CELL_ARR[@]}"; do _TPC_SUM=$((_TPC_SUM + T)); done
TOTAL_THREAD_DEMAND=$(( N_SPECIES * (1 + N_SOURCES * (N_MODES + N_BIN) * _TPC_SUM) ))
echo ">> cpus/task:     $T_TOTAL  (need ≥ $TOTAL_CELLS cells; ~$TOTAL_THREAD_DEMAND total threads at peak)"
echo ">> --fast mode:   $([[ $FAST_MODE -eq 1 ]] && echo "ON (adds _fast variant per taffy cell)" || echo "OFF (default column-walk only)")"
echo ">> --bin sizes:   $([[ $N_BIN -gt 0 ]] && echo "${BIN_SIZES_ARR[*]}  (adds 1 fast+bin cell per source per size)" || echo "OFF")"
echo ">> --threadsPerCell: ${THREADS_PER_CELL_ARR[*]}  (OMP_NUM_THREADS per taffy cell; N>1 cells get _t<N> suffix)"
(( T_TOTAL >= TOTAL_THREAD_DEMAND )) || \
    echo ">> WARN: T_TOTAL=$T_TOTAL < TOTAL_THREAD_DEMAND=$TOTAL_THREAD_DEMAND; cells will oversubscribe / queue" >&2
echo ">> time budget:   $TIME_BUDGET s per cell"
echo ">> local-stage:   $([[ $STAGE_LOCAL -eq 1 ]] && echo "ON (copies to \$TMPDIR)" || echo "OFF (reads from network)")"
[[ -n "$TMP_GB" ]] && echo ">> --tmp request: ${TMP_GB} GB per task"

# --- Generate the runner script (the thing sbatch executes). -----------
RUNNER="$OUTDIR/bench.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -uo pipefail

UNI="$UNI"
UNI_TAF="$UNI_TAF"
HAL="$HAL"
OUTDIR="$OUTDIR"
REF="$REF"
TREE="$TREE"     # for plot.py: divergence-from-ref x-axis
TIME_BUDGET=$TIME_BUDGET
STAGE_LOCAL=$STAGE_LOCAL
FAST_MODE=$FAST_MODE
T_TOTAL=$T_TOTAL
BIN_SIZES=( ${BIN_SIZES_ARR[*]:-} )
THREADS_PER_CELL=( ${THREADS_PER_CELL_ARR[*]} )

# Format an integer bin size as the human-readable suffix used in tool
# names ('1M', '100k'); divisibility-exact only, no rounding.
fmt_bin() {
    local n=\$1
    if   (( n % 1000000000 == 0 )); then echo "\$((n/1000000000))G"
    elif (( n % 1000000    == 0 )); then echo "\$((n/1000000))M"
    elif (( n % 1000       == 0 )); then echo "\$((n/1000))k"
    else                                  echo "\$n"
    fi
}
TAFFY="$TAFFY"
HALLIFTOVER="$HALLIFTOVER"
SPECIES_TSV="$OUTDIR/species.tsv"
NATIVE_BED="$OUTDIR/beds/genome.native.bed"
PREFIXED_BED="$OUTDIR/beds/genome.prefixed.bed"
N_CHROMS_INNER=$N_CHROMS
REF_BP_INNER=$REF_BP

BENCH_TSV="\$OUTDIR/bench.tsv"
LOGDIR="\$OUTDIR/logs"
mkdir -p "\$LOGDIR"

# --- Stage inputs to local scratch (\$TMPDIR or /tmp fallback). -----
# .tui (taffy never opens the MAF itself) + HAL.  Both are large so
# this is the dominant pre-bench cost; staging just once amortises it
# across the N_SPECIES * cells_per_species cells that follow.
if [[ "\$STAGE_LOCAL" -eq 1 ]]; then
    SCRATCH="\${TMPDIR:-/tmp/taffy_lift_genome_bench_\${SLURM_JOB_ID:-\$\$}}"
    STAGE_DIR="\$SCRATCH/taffy_lift_genome_stage_\${SLURM_JOB_ID:-\$\$}"
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
    # .tui (taffy lift never opens the source MAF/TAF itself).
    # Each provided source (-u, --uniTaf) gets its .tui staged + a 0-byte
    # stub at <input> so taffy's tui_path() resolves to the staged file.
    if [[ -n "\$UNI" ]]; then
        stage_one "\$UNI.tui" > /dev/null
        LOCAL_UNI="\$STAGE_DIR/\$(basename "\$UNI")"
        : > "\$LOCAL_UNI"
        UNI="\$LOCAL_UNI"
    fi
    if [[ -n "\$UNI_TAF" ]]; then
        stage_one "\$UNI_TAF.tui" > /dev/null
        LOCAL_UNI_TAF="\$STAGE_DIR/\$(basename "\$UNI_TAF")"
        : > "\$LOCAL_UNI_TAF"
        UNI_TAF="\$LOCAL_UNI_TAF"
    fi

    # HAL for halLiftover
    HAL=\$(stage_one "\$HAL")

    echo "stage: all inputs staged to \$STAGE_DIR" >&2
fi

# Write header if file is empty.
if [[ ! -s "\$BENCH_TSV" ]]; then
    printf "tool\tgenome_id\tsci_name\tcommon_name\tn_chroms_in\tref_bp\twall_s\tpeak_rss_kb\texit\ttimed_out\tn_mapped\tn_mapped_bp\n" > "\$BENCH_TSV"
fi

# run_cell tool species_id sci common cmd...
# n_mapped = output bed line count; n_mapped_bp = sum of (end - start)
# over all output rows.
run_cell() {
    local tool="\$1" sid="\$2" sci="\$3" common="\$4"
    shift 4
    local stem="\${tool}_\${sid}"
    local time_file="\$LOGDIR/time_\${stem}.txt"
    local err_file="\$LOGDIR/err_\${stem}.log"
    local out_bed="\$LOGDIR/mapped_\${stem}.bed"

    /usr/bin/time -q -f '%e %M' -o "\$time_file" \\
        timeout --signal=KILL "\$TIME_BUDGET" "\$@" \\
        > /dev/null 2> "\$err_file"
    local rc=\$?

    local wall rss timed_out=0
    if [[ -s "\$time_file" ]]; then
        read -r wall rss < "\$time_file"
        if ! [[ "\$wall" =~ ^[0-9.]+\$ && "\$rss" =~ ^[0-9]+\$ ]]; then
            read -r wall rss < <(awk '/^[0-9.]+[ \t][0-9]+\$/ {l=\$0} END{print l}' "\$time_file")
            [[ -z "\$wall" ]] && { wall="NA"; rss="NA"; }
        fi
    else
        wall="NA"; rss="NA"
    fi
    (( rc == 137 || rc == 124 )) && timed_out=1

    local n_mapped=0 n_mapped_bp=0
    if [[ -s "\$out_bed" ]]; then
        # For BED output bp = sum(end-start) (target bp covered); for
        # bedGraph from --bin variants \$3-\$2 is the constant bin width
        # and bin_size * n_rows is meaningless -- the actual source-bp
        # count is in the value column (\$4).  Branch on tool name.
        if [[ "\$tool" == *_bin* ]]; then
            read -r n_mapped n_mapped_bp < <(awk '{n++; bp += \$4} END {print n+0, bp+0}' "\$out_bed")
        else
            read -r n_mapped n_mapped_bp < <(awk '{n++; bp += \$3 - \$2} END {print n+0, bp+0}' "\$out_bed")
        fi
    fi

    printf "%s\t%s\t%s\t%s\t\$N_CHROMS_INNER\t\$REF_BP_INNER\t%s\t%s\t%d\t%d\t%d\t%d\n" \\
        "\$tool" "\$sid" "\$sci" "\$common" "\$wall" "\$rss" "\$rc" "\$timed_out" "\$n_mapped" "\$n_mapped_bp"
}

# --- Concurrency throttle ----------------------------------------------
# Bound the SUM of running-cell OMP_NUM_THREADS by T_TOTAL so we never
# oversubscribe the SLURM allocation -- same as lift-bench.  Without
# this, the 9 species x (halLiftover + N modes x M thread counts x
# K bins) cells all launch at once and oversubscribe the cpus-per-task
# allocation.
THREAD_BUDGET=\$T_TOTAL
launched=0
declare -A pid_threads

acquire_slot() {
    local threads=\$1
    while (( launched + threads > THREAD_BUDGET )); do
        wait -n 2>/dev/null || true
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

# --- Fire all cells in parallel.  No waves. ----------------------------
# Cells per species: 1 halLiftover + 1 per provided .tui source (1-2).
echo "=== launching cells: \$(wc -l < "\$SPECIES_TSV") species x (1 halLiftover + .tui sources) ==="
t0=\$SECONDS
declare -A pids rowfiles

while IFS=\$'\t' read -r sid sci common; do
    [[ -n "\$sid" ]] || continue

    # ---- halLiftover cell ----
    # halLiftover HAL SRC_GENOME SRC_BED TGT_GENOME TGT_BED
    # Native (un-prefixed) bed -- halLiftover wants bare chrom names.
    stem_hl="halLiftover_\${sid}"
    rowfiles[\$stem_hl]="\$LOGDIR/row_\${stem_hl}.tsv"
    acquire_slot 1
    ( run_cell halLiftover "\$sid" "\$sci" "\$common" \\
        "\$HALLIFTOVER" "\$HAL" "\$REF" "\$NATIVE_BED" "\$sid" \\
        "\$LOGDIR/mapped_\${stem_hl}.bed" \\
      ) > "\${rowfiles[\$stem_hl]}" &
    pids[\$stem_hl]=\$!; register_pid \$! 1

    # ---- taffy lift cells: one per (.tui source) x (mode) ----
    # Mode loop: 'default' = column-walk (existing); 'fast' = chunk-walk
    # via --fast.  FAST_MODE=1 adds the _fast variant alongside; default
    # cells are still emitted so the bench compares both.  Tool name:
    #   default -> 'maf.tui' / 'taf.tui'           (existing names)
    #   fast    -> 'maf.tui_fast' / 'taf.tui_fast'
    # Prefixed bed -- taffy lift -b expects "<genome>.<chrom>".
    # ttag(): tool-name suffix for OMP_NUM_THREADS=N (empty when N==1
    # so the bare "maf.tui_fast" name still appears in single-thread
    # runs -- matters for legacy bench.tsv readers).
    ttag() { [[ "\$1" == "1" ]] && echo "" || echo "_t\$1"; }

    modes=( default )
    [[ "\$FAST_MODE" -eq 1 ]] && modes+=( fast )
    for mode in "\${modes[@]}"; do
        if [[ "\$mode" == "fast" ]]; then
            ftag="_fast"; ffl=( --fast )
        else
            ftag=""; ffl=()
        fi
        for THREADS in "\${THREADS_PER_CELL[@]}"; do
            tt="\$(ttag \$THREADS)"
            if [[ -n "\$UNI" ]]; then
                tool="maf.tui\${ftag}\${tt}"
                stem="\${tool}_\${sid}"
                rowfiles[\$stem]="\$LOGDIR/row_\${stem}.tsv"
                acquire_slot \$THREADS
                ( export OMP_NUM_THREADS=\$THREADS; \\
                  run_cell "\$tool" "\$sid" "\$sci" "\$common" \\
                    "\$TAFFY" lift -i "\$UNI" -b "\$PREFIXED_BED" -g "\$sid" \\
                                  "\${ffl[@]}" \\
                                  -o "\$LOGDIR/mapped_\${stem}.bed" \\
                ) > "\${rowfiles[\$stem]}" &
                pids[\$stem]=\$!; register_pid \$! \$THREADS
            fi
            if [[ -n "\$UNI_TAF" ]]; then
                tool="taf.tui\${ftag}\${tt}"
                stem="\${tool}_\${sid}"
                rowfiles[\$stem]="\$LOGDIR/row_\${stem}.tsv"
                acquire_slot \$THREADS
                ( export OMP_NUM_THREADS=\$THREADS; \\
                  run_cell "\$tool" "\$sid" "\$sci" "\$common" \\
                    "\$TAFFY" lift -i "\$UNI_TAF" -b "\$PREFIXED_BED" -g "\$sid" \\
                                  "\${ffl[@]}" \\
                                  -o "\$LOGDIR/mapped_\${stem}.bed" \\
                ) > "\${rowfiles[\$stem]}" &
                pids[\$stem]=\$!; register_pid \$! \$THREADS
            fi
        done
    done

    # --bin variants: --fast --bin <B>, one cell per (source, bin size,
    # threads-per-cell).
    if [[ \${#BIN_SIZES[@]} -gt 0 ]]; then
        for B in "\${BIN_SIZES[@]}"; do
            btag="_bin\$(fmt_bin \$B)"
            for THREADS in "\${THREADS_PER_CELL[@]}"; do
                tt="\$(ttag \$THREADS)"
                if [[ -n "\$UNI" ]]; then
                    tool="maf.tui_fast\${btag}\${tt}"
                    stem="\${tool}_\${sid}"
                    rowfiles[\$stem]="\$LOGDIR/row_\${stem}.tsv"
                    acquire_slot \$THREADS
                    ( export OMP_NUM_THREADS=\$THREADS; \\
                      run_cell "\$tool" "\$sid" "\$sci" "\$common" \\
                        "\$TAFFY" lift -i "\$UNI" -b "\$PREFIXED_BED" -g "\$sid" \\
                                      --fast --bin "\$B" \\
                                      -o "\$LOGDIR/mapped_\${stem}.bed" \\
                    ) > "\${rowfiles[\$stem]}" &
                    pids[\$stem]=\$!; register_pid \$! \$THREADS
                fi
                if [[ -n "\$UNI_TAF" ]]; then
                    tool="taf.tui_fast\${btag}\${tt}"
                    stem="\${tool}_\${sid}"
                    rowfiles[\$stem]="\$LOGDIR/row_\${stem}.tsv"
                    acquire_slot \$THREADS
                    ( export OMP_NUM_THREADS=\$THREADS; \\
                      run_cell "\$tool" "\$sid" "\$sci" "\$common" \\
                        "\$TAFFY" lift -i "\$UNI_TAF" -b "\$PREFIXED_BED" -g "\$sid" \\
                                      --fast --bin "\$B" \\
                                      -o "\$LOGDIR/mapped_\${stem}.bed" \\
                    ) > "\${rowfiles[\$stem]}" &
                    pids[\$stem]=\$!; register_pid \$! \$THREADS
                fi
            done
        done
    fi
done < "\$SPECIES_TSV"

# Wait for all cells, append rows in stable order.
for stem in "\${!pids[@]}"; do
    wait "\${pids[\$stem]}" || true
    [[ -s "\${rowfiles[\$stem]}" ]] && cat "\${rowfiles[\$stem]}" >> "\$BENCH_TSV"
done

echo "=== all cells done in \$((SECONDS - t0)) s ==="
echo "bench done.  TSV: \$BENCH_TSV"
EOF
chmod +x "$RUNNER"

# --- Companion plot script. -------------------------------------------
PLOT="$OUTDIR/plot.py"
cat > "$PLOT" <<'PY'
#!/usr/bin/env python3
"""Plot wall + RSS + bp-lifted for the whole-genome bench.  X-axis is
species rank by divergence from REF (evenly spaced)."""
import csv, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

bench_dir = os.path.dirname(os.path.abspath(__file__))

# Tree from bench.sh
tree_path = None
ref = None
with open(os.path.join(bench_dir, "bench.sh")) as f:
    for line in f:
        m = re.match(r'^TREE="([^"]+)"', line.strip())
        if m: tree_path = m.group(1)
        m = re.match(r'^REF="([^"]+)"', line.strip())
        if m: ref = m.group(1)
if not tree_path or not os.path.exists(tree_path):
    sys.exit(f"ERROR: tree not found: {tree_path}")
if not ref:
    sys.exit("ERROR: REF not found in bench.sh")

# Parse newick + compute leaf distances.
text = open(tree_path).read().strip()
if text.endswith(';'): text = text[:-1]
i = [0]
def parse_clade():
    children = []
    if text[i[0]] == '(':
        i[0] += 1
        while True:
            children.append(parse_clade())
            if text[i[0]] == ',': i[0] += 1
            elif text[i[0]] == ')': i[0] += 1; break
    j = i[0]
    while j < len(text) and text[j] not in '(),:;': j += 1
    name = text[i[0]:j]; i[0] = j
    bl = 0.0
    if i[0] < len(text) and text[i[0]] == ':':
        i[0] += 1; k = i[0]
        while k < len(text) and text[k] not in '(),;': k += 1
        bl = float(text[i[0]:k]); i[0] = k
    return (name, bl, children)
root = parse_clade()
parent = {}; nm = {}; nid = [0]
def annotate(n, p=None):
    x = nid[0]; nid[0] += 1
    name, bl, ch = n
    parent[x] = (p, name, bl, ch)
    if not ch: nm[name] = x
    for c in ch: annotate(c, x)
annotate(root)
def root_path(li):
    out = []; cur = li; cum = 0.0
    while True:
        out.append((cur, cum))
        p, _, bl, _ = parent[cur]
        if p is None: break
        cur = p; cum += bl
    return out
def dist(a, b):
    pa = {n: c for n, c in root_path(a)}
    for n, c in root_path(b):
        if n in pa: return pa[n] + c
    return None
ref_id = nm[ref]

# Load bench.tsv
rows = []
with open(os.path.join(bench_dir, "bench.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            r["wall_s"] = float(r["wall_s"]) if r["wall_s"] != "NA" else None
            r["peak_rss_kb"] = float(r["peak_rss_kb"]) if r["peak_rss_kb"] != "NA" else None
            r["timed_out"] = int(r["timed_out"])
            r["n_mapped"] = int(r["n_mapped"])
            r["n_mapped_bp"] = int(r["n_mapped_bp"])
            r["ref_bp"] = int(r["ref_bp"])
            g = r["genome_id"]
            r["dist"] = dist(ref_id, nm[g]) if g in nm else None
            if r["dist"] is not None:
                rows.append(r)
        except (ValueError, KeyError):
            continue

species_order = sorted({(r["dist"], r["genome_id"], r["common_name"]) for r in rows})
rank = {gid: i for i, (_, gid, _) in enumerate(species_order)}
ref_bp = rows[0]["ref_bp"] if rows else 0

# 3-panel: wall, RSS, bp.  Bar chart with two bars per species (tool).
fig, axes = plt.subplots(3, 1, figsize=(14, 14))
ax_wall, ax_rss, ax_bp = axes
fig.subplots_adjust(left=0.10, right=0.97, top=0.95, bottom=0.10, hspace=0.40)

BASE_COLOR = {
    "halLiftover": "#2ca02c",
    "maf.tui":     "#1f77b4",
    "taf.tui":     "#9467bd",
}
def _parse_bin_suffix(s):
    m = re.match(r'^(\d+)([GMk]?)$', s)
    if not m: return None
    return int(m.group(1)) * {'G': 10**9, 'M': 10**6, 'k': 10**3, '': 1}[m.group(2)]
def parse_tool(tool):
    """Return (family, is_fast, bin_size, threads).  Strips suffixes
    outermost-first: _t<N>, _bin<S>, _fast.  Bin variants stay in the
    same colour family as the column-walk sibling."""
    threads = 1
    m = re.search(r'_t(\d+)$', tool)
    if m:
        threads = int(m.group(1))
        tool = tool[:m.start()]
    bin_size = None
    m = re.search(r'_bin(\d+[GMk]?)$', tool)
    if m:
        bin_size = _parse_bin_suffix(m.group(1))
        tool = tool[:m.start()]
    is_fast = tool.endswith("_fast")
    base = tool[:-len("_fast")] if is_fast else tool
    return base, is_fast, bin_size, threads
def tool_color(tool):
    fam, _, _, _ = parse_tool(tool)
    return BASE_COLOR.get(fam, "#888888")
# Order: halLiftover, maf.tui (column-walk + fast paired + bin variants),
# taf.tui (same), anything else.  Within family: default-before-fast,
# no-bin-before-bin (ascending bin size), ascending thread count.
def tool_sort_key(t):
    fam, is_fast, bin_size, threads = parse_tool(t)
    fam_order = {"halLiftover": 0, "maf.tui": 1, "taf.tui": 2}.get(fam, 9)
    return (fam_order, int(is_fast),
            -1 if bin_size is None else bin_size, threads)
tools = sorted(set(r["tool"] for r in rows), key=tool_sort_key)

xs = np.arange(len(species_order))
n_tools = max(len(tools), 1)
width = 0.8 / n_tools

for ti, tool in enumerate(tools):
    by_rank = {rank[r["genome_id"]]: r for r in rows
               if r["tool"] == tool and r["wall_s"] is not None}
    wt_done = [None]*len(xs); wt_to = [None]*len(xs)
    rs_done = [None]*len(xs); rs_to = [None]*len(xs)
    bp_done = [None]*len(xs); bp_to = [None]*len(xs)
    for i in xs:
        r = by_rank.get(i)
        if not r: continue
        bp_clamped = max(r["n_mapped_bp"], 0)
        if r["timed_out"]:
            wt_to[i] = r["wall_s"]; rs_to[i] = r["peak_rss_kb"]/1024.0; bp_to[i] = bp_clamped
        else:
            wt_done[i] = r["wall_s"]; rs_done[i] = r["peak_rss_kb"]/1024.0; bp_done[i] = bp_clamped
    z = lambda lst: [v if v is not None else 0 for v in lst]
    off = (ti - (n_tools - 1) / 2) * width
    color = tool_color(tool)
    _, is_fast, bin_size, _ = parse_tool(tool)
    # Three visual styles:
    #   default        -> solid filled bar (family colour)
    #   --fast         -> hollow / outlined bar (same colour)
    #   --fast --bin N -> hollow + diagonal stripes ('xx' hatch), same colour
    # Timed-out variant of each keeps its base style + a '//' hatch.
    if bin_size:
        bar_kw = dict(facecolor="none", edgecolor=color,
                      hatch="xx", linewidth=1.0)
        to_kw  = dict(facecolor="none", edgecolor=color,
                      hatch="xx//", linewidth=1.0, alpha=0.85)
    elif is_fast:
        bar_kw = dict(facecolor="none", edgecolor=color, linewidth=1.5)
        to_kw  = dict(facecolor="none", edgecolor=color,
                      hatch="//", linewidth=1.5, alpha=0.85)
    else:
        bar_kw = dict(color=color)
        to_kw  = dict(color=color, hatch="//", edgecolor="black",
                      linewidth=0.5, alpha=0.55)
    # Suppress the "(timed out)" legend entry when no cell actually
    # timed out for this tool -- otherwise the legend doubles in length.
    # Use the same label on all three axes so per-panel legends agree.
    has_to = any(v is not None for v in wt_to)
    to_label = f"{tool} (timed out)" if has_to else "_nolegend_"
    ax_wall.bar(xs + off, z(wt_done), width=width, label=tool, **bar_kw)
    ax_wall.bar(xs + off, z(wt_to),   width=width, label=to_label, **to_kw)
    ax_rss.bar (xs + off, z(rs_done), width=width, label=tool, **bar_kw)
    ax_rss.bar (xs + off, z(rs_to),   width=width, label=to_label, **to_kw)
    ax_bp.bar  (xs + off, z(bp_done), width=width, label=tool, **bar_kw)
    ax_bp.bar  (xs + off, z(bp_to),   width=width, label=to_label, **to_kw)

xlabels = [f"{c}\n({d:.2f})" for d, _, c in species_order]
for ax, title, ylab, yscale in [
    (ax_wall, f"wall time -- whole {ref} genome ({ref_bp/1e9:.2f} Gb input)", "seconds", "log"),
    (ax_rss,  "peak RSS", "MB", "linear"),
    (ax_bp,   "total target bp lifted", "bp mapped", "linear"),
]:
    ax.set_yscale(yscale)
    ax.set_xticks(xs)
    ax.set_xticklabels(xlabels, rotation=30, ha="right", fontsize=10)
    ax.set_xlabel(f"target species (ranked by divergence from {ref})")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True, which="both", axis="y", alpha=0.3)
    ax.legend(fontsize=10, loc="upper left")

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
    -J taffy_lift_genome_bench
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
