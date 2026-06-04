#!/bin/bash
#
# taffy lift benchmark driver -- SLURM
#
# Benchmarks two interval-lift tools against a panel of species at
# varying divergence from the reference (default: chicken):
#
#   liftover  -- UCSC liftOver with a per-pair .chain.gz
#   taffy     -- taffy lift -b with the universal MAF's .tui index
#
# For each species S and interval-size N we lift the SAME randomly
# generated bed of K intervals (length N each) from REF to S and record
# wall + peak-RSS of the lift call.  Mean per-interval time = wall / K.
#
# X-axis on the plot is divergence from REF (sum of branch lengths on
# the species tree); Y is wall (panel 1) and peak RSS (panel 2).  Each
# size N is a separate line / marker shape.
#
# Inputs to stage (default ON, copies to $TMPDIR):
#   UNI.tui                  -- only the .tui; taffy lift does not read the MAF
#   <REF>_vs_<S>.chain.gz    -- one per species in the panel
#
# Default species panel (override with -L FILE; format: GCA<TAB>scientific<TAB>common):
#   GCA_039878825.1  Coturnix_chinensis         King_Quail
#   GCF_963932015.1  Anas_acuta                 Pintail_Duck
#   GCF_028018845.1  Melospiza_georgiana        Swamp_Sparrow
#   GCF_040807025.1  Struthio_camelus_australis Ostrich
#   GCF_028858775.2  Pan_troglodytes            Chimpanzee
#   GCA_000001405.15 hg38                       Human
#   GCA_000001635.9  mm39                       House_Mouse
#   GCF_037176765.1  Anolis_sagrei              Brown_Anole
#   GCA_944039275.1  Danio_rerio                Zebrafish
#   GCF_905171765.1  Bufo_bufo                  Common_Toad
#
# Run model: one SLURM job; outer loop = size waves (default 3); within
# each wave all species x 2 tools = up to 20 cells run in parallel.  Each
# cell is single-threaded.
#
# Usage:
#   taffy_lift_bench_slurm.sh \
#       -u UNI.maf.gz   -c CHAINS_DIR  -t TREE.nwk \
#       -o OUTDIR  [options]

set -euo pipefail

UNI=""
UNI_TAF=""
HAL=""
CHAINS_DIR=""
TREE=""
OUTDIR=""
T_TOTAL=24                          # cpus-per-task; one core per concurrent cell
REF="GCF_016700215.2"               # chicken
MAX_GAPS_CSV=""                     # --maxGap CSV: e.g. "0,1000,10000".  Empty = single taffy cell (K=0).
FAST_MODE=0                         # --fast: add a --fast variant cell next to each taffy cell.
BIN_SIZES_CSV=""                    # --bin CSV: e.g. "100000,1000000".  Empty = no binned variants.
                                    # Non-empty implies --fast (the binned mode requires it at the taffy CLI).
THREADS_PER_CELL_CSV="1"            # OMP_NUM_THREADS values for taffy cells, comma-separated.
                                    # Each value becomes a variant cell with OMP_NUM_THREADS=N
                                    # exported -- lets you A/B compare parallel-decode wall.
                                    # Default "1" = single-threaded (preserves prior behaviour).
                                    # Tool name suffix '_t<N>' when N>1; N=1 keeps the bare name.
SIZES="1000,100000,1000000"
N_INTERVALS=100                     # random intervals per (species, size) cell
SEED=42                             # RNG seed for the random-interval generator
TIME_BUDGET=3600                    # per-cell wall cap (timeout)
SBATCH_TIME=24
SBATCH_MEM=64
TMP_GB=""
STAGE_LOCAL=1
SPECIES_FILE=""                     # optional override (see format above)
PARTITION=""
ACCOUNT=""
DRY_RUN=0
WAIT=1
LIFTOVER="${LIFTOVER:-$(command -v liftOver || true)}"
HALLIFTOVER="${HALLIFTOVER:-$(command -v halLiftover || true)}"
TAFFY="${TAFFY:-$(command -v taffy || true)}"

# Default 10-species panel, baked here so the script is one-file portable.
# Tab-separated: <genome_id>\t<scientific>\t<english>.
DEFAULT_SPECIES_TSV=$(cat <<'EOF'
GCA_039878825.1	Coturnix_chinensis	King_Quail
GCF_963932015.1	Anas_acuta	Pintail_Duck
GCF_028018845.1	Melospiza_georgiana	Swamp_Sparrow
GCF_040807025.1	Struthio_camelus_australis	Ostrich
GCF_028858775.2	Pan_troglodytes	Chimpanzee
GCA_000001405.15	hg38	Human
GCA_000001635.9	mm39	House_Mouse
GCF_037176765.1	Anolis_sagrei	Brown_Anole
GCA_944039275.1	Danio_rerio	Zebrafish
GCF_905171765.1	Bufo_bufo	Common_Toad
EOF
)

usage() {
    cat >&2 <<EOF
taffy_lift_bench_slurm.sh -- liftover vs taffy lift benchmark across a
                             panel of species at varying divergence

Required (at least one of -u / --uniTaf must be set):
  -u FILE       Universal MAF (.uni.maf.gz; only its .tui sibling is read).
                If set, generates a "maf.tui" cell per species (taffy lift
                against the MAF-anchored .tui).
  --uniTaf FILE Universal TAF (.uni.taf.gz; only its .tui sibling is read).
                If set, generates a "taf.tui" cell per species (taffy lift
                against the TAF-anchored .tui).  Both flags can be set
                together to bench both .tui formats side by side.
  -H FILE       HAL file (for halLiftover)
  -c DIR        Directory of chain files named <REF>_vs_<GENOME_ID>.chain.gz
  -t FILE       Species tree (.nwk) -- used by the plot script for divergence
  -o DIR        Output directory

Optional:
  -r ID         Reference genome ID (default $REF; must match chain tName prefix)
  -T INT        cpus-per-task (default $T_TOTAL; one core per concurrent cell)
  -S CSV        Interval sizes in bp (default $SIZES)
  -N INT        # random intervals per (species, size) cell (default $N_INTERVALS)
  --seed INT    RNG seed for interval generation (default $SEED)
  -L FILE       Override species panel.  Format: 3 tab-separated cols per line:
                <genome_id>\\t<scientific>\\t<english>.  '#' = comment.
                Useful for smoke tests: pass a file with just 1-2 species.
  --maxGap CSV  taffy lift --maxGap values to bench side-by-side.  Default
                empty = one taffy cell at K=0 (current behaviour).  Each K
                in the CSV becomes its own taffy cell per (species, size,
                source).  Tool name in bench.tsv: 'maf.tui' / 'taf.tui' for
                K=0 (preserves existing names); 'maf.tui_g<K>' /
                'taf.tui_g<K>' for K > 0.  Example: --maxGap 0,1000,10000
                runs 3 taffy variants per cell.
  --fast        For each taffy cell (default + per-K variant if --maxGap is
                set), also launch a --fast variant using the chunk-iteration
                lift path (10-50x faster on multi-Mb queries).  Tool name
                in bench.tsv gets a '_fast' suffix: e.g. 'maf.tui_fast',
                'maf.tui_g1000_fast'.  Lets you compare default-vs-fast
                wall + verify output parity in one job.
  --bin CSV     taffy lift --bin sizes (bedGraph, coarse-grained browser
                tracks) to bench side-by-side.  Each value adds a
                --fast --bin <N> variant cell PER SPECIES PER SIZE (only
                at K=0 -- --bin is mutex with --maxGap at the CLI).
                Auto-enables --fast.  Tool name suffix: '_bin<S>' where
                <S> is human-readable when divisible (1M, 100k, 1G) and
                raw int otherwise.  Example: --bin 100000,1000000 adds
                'maf.tui_fast_bin100k' and 'maf.tui_fast_bin1M' per cell.
  --threadsPerCell CSV
                OMP_NUM_THREADS values for taffy cells, comma-separated
                (e.g. "1,4,8").  Each value becomes its own variant cell
                per (species, size, source, mode) -- lets you A/B compare
                parallel-decode wall in a single bench.  Tool name suffix
                '_t<N>' when N>1; N=1 keeps the bare name for compat.
                Default "1" = single-threaded.  halLiftover / liftover
                ignore OMP_NUM_THREADS so don't get _t<N> variants.
                T_TOTAL should be >= concurrent-cells * max-N to avoid
                oversubscription; the driver warns at submit time.
  --timeBudget SEC  Per-cell wall cap (timeout) (default $TIME_BUDGET)
  --time HRS    sbatch wall (default $SBATCH_TIME)
  --mem GB      sbatch mem (default $SBATCH_MEM)
  --tmp GB      Per-task local scratch requirement (sbatch --tmp=N).
                Default unset.
  --no-stage-local
                Skip the copy of .tui + chains to \$TMPDIR (read straight
                from the network paths).  Only sensible for small tests.
  --partition X --account X
  --no-wait     Submit and detach (default: driver blocks until SLURM done)
  --dry-run     Print sbatch; do not submit
  -h            Help

Override binary paths via env: TAFFY, LIFTOVER, HALLIFTOVER
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)             UNI="$2"; shift 2;;
        --uniTaf)       UNI_TAF="$2"; shift 2;;
        -H)             HAL="$2"; shift 2;;
        -c)             CHAINS_DIR="$2"; shift 2;;
        -t)             TREE="$2"; shift 2;;
        -o)             OUTDIR="$2"; shift 2;;
        -r)             REF="$2"; shift 2;;
        -T)             T_TOTAL="$2"; shift 2;;
        -S)             SIZES="$2"; shift 2;;
        -N)             N_INTERVALS="$2"; shift 2;;
        --seed)         SEED="$2"; shift 2;;
        -L)             SPECIES_FILE="$2"; shift 2;;
        --maxGap)       MAX_GAPS_CSV="$2"; shift 2;;
        --fast)         FAST_MODE=1; shift;;
        --bin)          BIN_SIZES_CSV="$2"; shift 2;;
        --threadsPerCell) THREADS_PER_CELL_CSV="$2"; shift 2;;
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

for v in HAL CHAINS_DIR TREE OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: -$(echo $v | cut -c1) required" >&2; usage 1; }
done
[[ -n "$UNI" || -n "$UNI_TAF" ]] || {
    echo "ERROR: at least one of -u / --uniTaf must be set" >&2; usage 1;
}
[[ -n "$TAFFY"       ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$LIFTOVER"    ]] || { echo "ERROR: liftOver not on PATH (set \$LIFTOVER)" >&2; exit 1; }
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
[[ -d "$CHAINS_DIR" ]] || { echo "ERROR: $CHAINS_DIR not found" >&2; exit 1; }
[[ -f "$TREE"     ]] || { echo "ERROR: $TREE not found" >&2; exit 1; }

mkdir -p "$OUTDIR" "$OUTDIR/logs" "$OUTDIR/beds"
echo ">> driver starting (output dir: $OUTDIR)" >&2

# --- Resolve species panel ----------------------------------------------
SPECIES_TSV="$OUTDIR/species.tsv"
if [[ -n "$SPECIES_FILE" ]]; then
    # Strip '#' comments + blank lines; canonicalize to tab-separated.
    # Default awk FS = whitespace, so tab- AND space-separated input both work.
    awk '!/^#/ && NF >= 3 {print $1"\t"$2"\t"$3}' "$SPECIES_FILE" > "$SPECIES_TSV"
else
    printf "%s\n" "$DEFAULT_SPECIES_TSV" > "$SPECIES_TSV"
fi
N_SPECIES=$(wc -l < "$SPECIES_TSV")
[[ "$N_SPECIES" -gt 0 ]] || { echo "ERROR: empty species panel" >&2; exit 1; }
echo ">> species panel: $N_SPECIES entries" >&2

# --- Resolve & validate chains for every species ------------------------
declare -A CHAIN_OF       # species_id -> chain.gz path
MISSING=0
while IFS=$'\t' read -r sid sci common; do
    [[ -n "$sid" ]] || continue
    ch="$CHAINS_DIR/${REF}_vs_${sid}.chain.gz"
    if [[ -f "$ch" ]]; then
        CHAIN_OF[$sid]="$ch"
    else
        echo "ERROR: missing chain for $sid ($sci): $ch" >&2
        MISSING=1
    fi
done < "$SPECIES_TSV"
(( MISSING == 0 )) || exit 1

# --- Reference chrom sizes (sourced from the .tui) ----------------------
# Source of truth for chicken chrom names is the .tui itself (those are
# what taffy lift will look up); we ASSUME the chains' tName uses the
# same names (cactus.chain inherits chrom naming from the HAL).  If the
# chain disagrees, the liftOver cells will exit nonzero -- check the
# logs/err_liftover_* files.  Different leaf assemblies use different
# naming conventions in VGP cactus (UCSC chr1 for mm39, NCBI accessions
# like NC_051216.1 for VGP species), so we cannot guess.
REF_SIZES="$OUTDIR/ref.sizes"
# Prefer the MAF.tui if available; the chrom list is identical from either
# (both .tui files were built from the same alignment), but we only need
# to query one.
UNI_STATS_SRC="${UNI:-$UNI_TAF}"
echo ">> querying .tui ($UNI_STATS_SRC) for ${REF}.* chroms via taffy stats -s ..." >&2
"$TAFFY" stats -i "$UNI_STATS_SRC" -s \
    | awk -v p="${REF}." 'index($1, p) == 1 {sub(p, "", $1); print $1"\t"$2}' \
    | sort -k2,2nr > "$REF_SIZES"
[[ -s "$REF_SIZES" ]] || {
    echo "ERROR: no sequences matching ${REF}.* in $UNI_STATS_SRC's .tui." >&2
    echo "       Check that -r matches the genome prefix used in the universal MAF." >&2
    exit 1
}
echo ">>   found $(wc -l < "$REF_SIZES") chroms matching ${REF}.*" >&2

# Quick cross-check: does at least one chain header tName match one of
# our REF_SIZES chroms?  Catches the bench-killing case where chain
# names disagree with .tui names (e.g., chr1 vs NC_006088.5) BEFORE the
# job runs.  Loop over chains since the first one may legitimately have
# 0 chrom intersection (an alt-only chrom pair) if we're unlucky.
#
# Subshell with pipefail off + `|| true`: `zcat | awk | head -200` early-
# closes the pipe, which SIGPIPEs zcat (exit 141).  With outer pipefail+
# set -e, that would silently abort the whole driver.  Disabling pipefail
# inside the substitution scopes the relaxation to this one check.
echo ">> validating chain tName ↔ .tui chrom intersection ..." >&2
INTERSECT_OK=0
for ch in "${CHAIN_OF[@]}"; do
    CHAIN_TNAMES=$(set +o pipefail; \
        zcat "$ch" 2>/dev/null | awk '/^chain/ {print $3}' | head -200 | sort -u; \
        true)
    if [[ -n "$CHAIN_TNAMES" ]] \
       && awk 'NR==FNR{t[$0]=1; next} ($1 in t)' \
              <(echo "$CHAIN_TNAMES") "$REF_SIZES" \
              | grep -q .; then
        INTERSECT_OK=1
        break
    fi
done
if [[ "$INTERSECT_OK" -ne 1 ]]; then
    echo "WARN: no overlap between chain tNames and ${REF}.* in .tui." >&2
    echo "      First chain head chroms (sample):" >&2
    first_chain=$(printf '%s\n' "${CHAIN_OF[@]}" | head -1)
    ( set +o pipefail; \
        zcat "$first_chain" 2>/dev/null \
            | awk '/^chain/ {print "       "$3}' | head -5 >&2; \
        true )
    echo "      .tui ${REF}.* chroms (sample):" >&2
    head -5 "$REF_SIZES" | awk '{print "       "$1}' >&2
    echo "      liftOver cells will likely fail; taffy cells should still work." >&2
else
    echo ">>   chain ↔ .tui chrom naming OK" >&2
fi

# --- Generate the bed files (one per size; shared across species) -------
# Each line: <chrom>\t<start>\t<end>.  Chrom names come from the chain
# tName (bare accession; e.g. NC_006088.5).  For taffy lift we prefix
# with "<REF>." on-the-fly inside the runner.  Picks chroms weighted by
# size, then offsets uniformly inside.
IFS=',' read -r -a SIZE_ARR <<< "$SIZES"

# --- Resolve --maxGap CSV.  Default empty -> a single taffy cell at K=0
# named 'maf.tui' / 'taf.tui' (existing behaviour).  Non-empty CSV: each
# K becomes its own taffy cell suffixed '_g<K>' (K=0 keeps the
# unsuffixed name for compat).
if [[ -z "$MAX_GAPS_CSV" ]]; then
    MAX_GAPS_ARR=( 0 )
else
    IFS=',' read -r -a MAX_GAPS_ARR <<< "$MAX_GAPS_CSV"
    for K in "${MAX_GAPS_ARR[@]}"; do
        [[ "$K" =~ ^[0-9]+$ ]] || { echo "ERROR: --maxGap values must be non-negative integers (got '$K')" >&2; exit 1; }
    done
fi

# --- Resolve --bin CSV.  Default empty -> no bin cells.  Non-empty
# implies --fast (auto-enable).  Bin cells launch only at K=0 since the
# taffy CLI rejects --bin with --maxGap; this mirrors that constraint
# without bothering to enumerate the cartesian product.
if [[ -z "$BIN_SIZES_CSV" ]]; then
    BIN_SIZES_ARR=( )
else
    IFS=',' read -r -a BIN_SIZES_ARR <<< "$BIN_SIZES_CSV"
    for B in "${BIN_SIZES_ARR[@]}"; do
        [[ "$B" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --bin values must be positive integers (got '$B')" >&2; exit 1; }
    done
    FAST_MODE=1
fi

# --- Resolve --threadsPerCell CSV.  Default "1" -> a single variant
# with OMP_NUM_THREADS=1 (preserves existing tool names).  Each value
# fans out a separate taffy cell with OMP_NUM_THREADS=N exported;
# tool name suffix '_t<N>' when N > 1, no suffix at N == 1.
IFS=',' read -r -a THREADS_PER_CELL_ARR <<< "$THREADS_PER_CELL_CSV"
for T in "${THREADS_PER_CELL_ARR[@]}"; do
    [[ "$T" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threadsPerCell values must be positive integers (got '$T')" >&2; exit 1; }
done

python3 - "$REF_SIZES" "$OUTDIR/beds" "$SEED" "$N_INTERVALS" "${SIZE_ARR[@]}" <<'PY'
import os, random, sys
sizes_path, beds_dir, seed, n_intervals, *sizes = sys.argv[1:]
seed = int(seed); n_intervals = int(n_intervals)
sizes = [int(s) for s in sizes]
chroms = []
with open(sizes_path) as f:
    for line in f:
        c, sz = line.rstrip().split('\t')
        chroms.append((c, int(sz)))
total = sum(sz for _, sz in chroms)
rng = random.Random(seed)
for S in sizes:
    # filter to chroms that can hold an interval of length S
    cand = [(c, sz) for c, sz in chroms if sz >= S]
    if not cand:
        print(f"WARN: no chrom big enough for size {S}", file=sys.stderr)
        open(os.path.join(beds_dir, f"intervals_{S}.bed"), 'w').close()
        continue
    cum = []
    acc = 0
    for c, sz in cand:
        acc += sz
        cum.append((acc, c, sz))
    with open(os.path.join(beds_dir, f"intervals_{S}.bed"), 'w') as out:
        for _ in range(n_intervals):
            r = rng.randrange(acc)
            for lim, c, sz in cum:
                if r < lim:
                    start = rng.randint(0, sz - S)
                    out.write(f"{c}\t{start}\t{start + S}\n")
                    break
    print(f"size {S}: wrote {n_intervals} intervals -> intervals_{S}.bed", file=sys.stderr)
PY

echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> liftOver:      $LIFTOVER"
echo ">> halLiftover:   $HALLIFTOVER"
echo ">> uni:           $UNI  (using \$UNI.tui only)"
echo ">> hal:           $HAL"
echo ">> chains dir:    $CHAINS_DIR"
echo ">> tree:          $TREE"
echo ">> reference:     $REF"
echo ">> species:       $N_SPECIES"
awk -F'\t' '{printf("                  %-18s %-30s %s\n", $1, $2, $3)}' "$SPECIES_TSV"
echo ">> sizes:         ${SIZE_ARR[*]}  ($N_INTERVALS intervals each, seed=$SEED)"
echo ">> --maxGap vals: ${MAX_GAPS_ARR[*]}  (taffy cells per source per (species, size))"
echo ">> --fast mode:   $([[ $FAST_MODE -eq 1 ]] && echo "ON (adds _fast variant per taffy cell)" || echo "OFF (default column-walk only)")"
echo ">> --bin sizes:   $([[ ${#BIN_SIZES_ARR[@]} -gt 0 ]] && echo "${BIN_SIZES_ARR[*]}  (adds 1 fast+bin cell per size at K=0)" || echo "OFF")"
echo ">> --threadsPerCell: ${THREADS_PER_CELL_ARR[*]}  (OMP_NUM_THREADS per taffy cell; N>1 cells get _t<N> suffix)"
# Oversubscription check: concurrent taffy cells = N_SPECIES * N_TAFFY_TOOLS;
# total CPU demand = sum(threads_per_cell) per (species, size, source, mode, bin).
# N_TAFFY_TOOLS = N_MAX_GAPS * (1 + FAST_MODE) + (FAST_MODE * N_BIN), times sources,
# all times threadsPerCell values length.  Approximate the order-of-magnitude check.
_T_SUM=0
for T in "${THREADS_PER_CELL_ARR[@]}"; do _T_SUM=$((_T_SUM + T)); done
_N_TAFFY_VARIANTS=$(( ${#MAX_GAPS_ARR[@]} * (1 + FAST_MODE) + FAST_MODE * ${#BIN_SIZES_ARR[@]} ))
_N_SRC=0; [[ -n "$UNI" ]] && _N_SRC=$((_N_SRC + 1)); [[ -n "$UNI_TAF" ]] && _N_SRC=$((_N_SRC + 1))
_THREADS_NEEDED=$(( N_SPECIES * _N_SRC * _N_TAFFY_VARIANTS * _T_SUM / ${#THREADS_PER_CELL_ARR[@]} ))
(( T_TOTAL >= _THREADS_NEEDED )) || echo ">> WARN: T_TOTAL=$T_TOTAL < estimated ${_THREADS_NEEDED} threads needed for full-parallel taffy fan-out; will queue / oversubscribe" >&2
echo ">> cpus/task:     $T_TOTAL  (one core per concurrent cell)"
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
CHAINS_DIR="$CHAINS_DIR"
OUTDIR="$OUTDIR"
REF="$REF"
TREE="$TREE"     # for plot.py: divergence-from-ref x-axis
T_TOTAL=$T_TOTAL
TIME_BUDGET=$TIME_BUDGET
STAGE_LOCAL=$STAGE_LOCAL
TAFFY="$TAFFY"
LIFTOVER="$LIFTOVER"
HALLIFTOVER="$HALLIFTOVER"
SPECIES_TSV="$OUTDIR/species.tsv"
SIZES=( ${SIZE_ARR[*]} )
MAX_GAPS=( ${MAX_GAPS_ARR[*]} )
FAST_MODE=$FAST_MODE
BIN_SIZES=( ${BIN_SIZES_ARR[*]:-} )
THREADS_PER_CELL=( ${THREADS_PER_CELL_ARR[*]} )

# Format an integer bin size as the human-readable suffix used in tool
# names ('1M', '100k', '500'); divisibility-exact only, no rounding.
fmt_bin() {
    local n=\$1
    if   (( n % 1000000000 == 0 )); then echo "\$((n/1000000000))G"
    elif (( n % 1000000    == 0 )); then echo "\$((n/1000000))M"
    elif (( n % 1000       == 0 )); then echo "\$((n/1000))k"
    else                                  echo "\$n"
    fi
}

BENCH_TSV="\$OUTDIR/bench.tsv"
LOGDIR="\$OUTDIR/logs"
BEDS="\$OUTDIR/beds"
mkdir -p "\$LOGDIR"

# --- Stage inputs to local scratch (\$TMPDIR or /tmp fallback). -----
# We stage the .tui (taffy lift never reads the MAF itself) + every
# chain.gz in the panel.  Trap-cleanup so an aborted job doesn't leak.
if [[ "\$STAGE_LOCAL" -eq 1 ]]; then
    SCRATCH="\${TMPDIR:-/tmp/taffy_lift_bench_\${SLURM_JOB_ID:-\$\$}}"
    STAGE_DIR="\$SCRATCH/taffy_lift_bench_stage_\${SLURM_JOB_ID:-\$\$}"
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
    # .tui only -- taffy lift never opens the source MAF/TAF itself.
    # Each provided source (-u, --uniTaf) gets its .tui staged.  The
    # input file itself isn't read; we just need a stub at <input> so
    # taffy lift's tui_path() resolves to the staged .tui.
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

    # HAL file for halLiftover.
    HAL=\$(stage_one "\$HAL")

    # Stage every chain in the panel into a local mirror directory.
    LOCAL_CHAINS="\$STAGE_DIR/chains"
    mkdir -p "\$LOCAL_CHAINS"
    while IFS=\$'\t' read -r sid sci common; do
        [[ -n "\$sid" ]] || continue
        src="\$CHAINS_DIR/\${REF}_vs_\${sid}.chain.gz"
        dst="\$LOCAL_CHAINS/\$(basename "\$src")"
        echo "stage: \$src -> \$dst (\$(stat -Lc %s "\$src" 2>/dev/null || echo ?) bytes)" >&2
        t0=\$SECONDS
        cp "\$src" "\$dst"
        echo "       done in \$((SECONDS - t0)) s" >&2
    done < "\$SPECIES_TSV"
    CHAINS_DIR="\$LOCAL_CHAINS"
    echo "stage: all inputs staged to \$STAGE_DIR" >&2
fi

# Write header if file is empty.
if [[ ! -s "\$BENCH_TSV" ]]; then
    printf "tool\tgenome_id\tsci_name\tcommon_name\tsize_bp\tn_intervals\twall_s\tpeak_rss_kb\texit\ttimed_out\tn_mapped\tn_mapped_bp\tn_unmapped\n" > "\$BENCH_TSV"
fi

# run_cell tool species_id sci common size budget cmd...
# Writes one TSV row to stdout.  n_mapped = output bed row count;
# n_mapped_bp = sum of (end - start) over all output rows; n_unmapped =
# count of input intervals that left no output (liftOver writes them to
# unmapped.bed; taffy/halLiftover never do, so for those tools n_unmapped
# is reported as 0 -- prefer n_mapped vs n_intervals to gauge coverage).
run_cell() {
    local tool="\$1" sid="\$2" sci="\$3" common="\$4" N="\$5" budget="\$6"
    shift 6
    local stem="\${tool}_\${sid}_\${N}"
    local time_file="\$LOGDIR/time_\${stem}.txt"
    local err_file="\$LOGDIR/err_\${stem}.log"
    local out_bed="\$LOGDIR/mapped_\${stem}.bed"
    local unm_bed="\$LOGDIR/unmapped_\${stem}.bed"

    /usr/bin/time -q -f '%e %M' -o "\$time_file" \\
        timeout --signal=KILL "\$budget" "\$@" \\
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

    local n_mapped=0 n_mapped_bp=0 n_unmapped=0
    if [[ -s "\$out_bed" ]]; then
        # Single awk pass for both row count and total bp.  For BED output
        # bp = sum(end-start) gives target bp covered; for bedGraph (the
        # --bin variants) the bin WIDTH would be constant and bin_size *
        # n_rows is meaningless -- the real source-bp count lives in the
        # value column (\$4).  Detect bin variants by tool-name suffix.
        if [[ "\$tool" == *_bin* ]]; then
            read -r n_mapped n_mapped_bp < <(awk '{n++; bp += \$4} END {print n+0, bp+0}' "\$out_bed")
        else
            read -r n_mapped n_mapped_bp < <(awk '{n++; bp += \$3 - \$2} END {print n+0, bp+0}' "\$out_bed")
        fi
    fi
    [[ -s "\$unm_bed" ]] && n_unmapped=\$(grep -v '^#' "\$unm_bed" | wc -l)

    printf "%s\t%s\t%s\t%s\t%d\t\$N_INTERVALS_INNER\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n" \\
        "\$tool" "\$sid" "\$sci" "\$common" "\$N" "\$wall" "\$rss" "\$rc" "\$timed_out" "\$n_mapped" "\$n_mapped_bp" "\$n_unmapped"
}

# Per cell N_intervals is fixed across the run (constant from driver).
N_INTERVALS_INNER=$N_INTERVALS

# --- Outer loop: size waves; inner = species x 3 tools in parallel. ----
for N in "\${SIZES[@]}"; do
    BED_NATIVE="\$BEDS/intervals_\${N}.bed"
    [[ -s "\$BED_NATIVE" ]] || { echo "skip wave N=\$N: empty bed"; continue; }
    BED_PREFIXED="\$LOGDIR/intervals_\${N}.prefixed.bed"
    awk -v p="\$REF." '{print p\$0}' "\$BED_NATIVE" > "\$BED_PREFIXED"

    echo "=== wave N=\$N  intervals=\$(wc -l < "\$BED_NATIVE") ==="
    t_wave=\$SECONDS
    declare -A pids rowfiles

    while IFS=\$'\t' read -r sid sci common; do
        [[ -n "\$sid" ]] || continue
        chain="\$CHAINS_DIR/\${REF}_vs_\${sid}.chain.gz"

        # ---- liftover cell ----
        stem_lo="liftover_\${sid}_\${N}"
        rowfiles[\$stem_lo]="\$LOGDIR/row_\${stem_lo}.tsv"
        out_lo="\$LOGDIR/mapped_\${stem_lo}.bed"
        unm_lo="\$LOGDIR/unmapped_\${stem_lo}.bed"
        # -minMatch=0: liftOver's default 0.95 rejects any interval where
        # fewer than 95% of bases can be mapped.  At cross-species + multi-
        # megabase scale that drops EVERY interval at sizes >= 100 kb.
        # taffy/halLiftover have no analogous threshold (they preserve any
        # mapping, even partial).
        # -multiple: keep all output regions when an interval maps to more
        # than one place in the target (paralogs / chain duplicates).
        # Default rejects the whole interval as "Duplicated in new" -- a
        # serious undercount at Mb scale, where most chicken regions have
        # at least one paralog hit in a close target.  taffy lift and
        # halLiftover both emit one row per paralog x per gap-free run, so
        # -multiple aligns liftover's semantics with theirs.
        ( run_cell liftover "\$sid" "\$sci" "\$common" "\$N" "\$TIME_BUDGET" \\
            "\$LIFTOVER" -minMatch=0 -multiple "\$BED_NATIVE" "\$chain" "\$out_lo" "\$unm_lo" \\
          ) > "\${rowfiles[\$stem_lo]}" &
        pids[\$stem_lo]=\$!

        # ---- halLiftover cell ----
        # halLiftover's input bed uses chain-native chrom names (no
        # <REF>. prefix), same as UCSC liftOver.  Argument order:
        # halLiftover HAL SRC_GENOME SRC_BED TGT_GENOME TGT_BED.
        stem_hl="halLiftover_\${sid}_\${N}"
        rowfiles[\$stem_hl]="\$LOGDIR/row_\${stem_hl}.tsv"
        out_hl="\$LOGDIR/mapped_\${stem_hl}.bed"
        ( run_cell halLiftover "\$sid" "\$sci" "\$common" "\$N" "\$TIME_BUDGET" \\
            "\$HALLIFTOVER" "\$HAL" "\$REF" "\$BED_NATIVE" "\$sid" "\$out_hl" \\
          ) > "\${rowfiles[\$stem_hl]}" &
        pids[\$stem_hl]=\$!

        # ---- taffy lift cells: one per (.tui source) x (--maxGap K) x (mode) -
        # Mode loop: 'default' = column-walk, 'fast' = chunk-walk (--fast).
        # When FAST_MODE=0 we only run 'default'; when 1 we run both so the
        # bench can compare side-by-side.  Tool name suffix:
        #   K=0, default      : "maf.tui"                  (existing name)
        #   K>0, default      : "maf.tui_g<K>"
        #   K=0, fast         : "maf.tui_fast"
        #   K>0, fast         : "maf.tui_g<K>_fast"
        # (Same scheme for taf.tui.)
        # ttag(): tool-name suffix for OMP_NUM_THREADS=N (empty for N==1
        # so the bare "maf.tui_fast" name still appears in single-thread
        # runs -- matters for legacy bench.tsv readers).
        ttag() { [[ "\$1" == "1" ]] && echo "" || echo "_t\$1"; }

        modes=( default )
        [[ "\$FAST_MODE" -eq 1 ]] && modes+=( fast )
        for K in "\${MAX_GAPS[@]}"; do
            if [[ "\$K" == "0" ]]; then gtag=""; else gtag="_g\${K}"; fi
            for mode in "\${modes[@]}"; do
                if [[ "\$mode" == "fast" ]]; then
                    ftag="_fast"; ffl=( --fast )
                else
                    ftag=""; ffl=()
                fi
                for THREADS in "\${THREADS_PER_CELL[@]}"; do
                    tt="\$(ttag \$THREADS)"
                    if [[ -n "\$UNI" ]]; then
                        tool="maf.tui\${gtag}\${ftag}\${tt}"
                        stem="\${tool}_\${sid}_\${N}"
                        rowfiles[\$stem]="\$LOGDIR/row_\${stem}.tsv"
                        ( export OMP_NUM_THREADS=\$THREADS; \\
                          run_cell "\$tool" "\$sid" "\$sci" "\$common" "\$N" "\$TIME_BUDGET" \\
                            "\$TAFFY" lift -i "\$UNI" -b "\$BED_PREFIXED" -g "\$sid" \\
                                          --maxGap "\$K" "\${ffl[@]}" \\
                                          -o "\$LOGDIR/mapped_\${stem}.bed" \\
                        ) > "\${rowfiles[\$stem]}" &
                        pids[\$stem]=\$!
                    fi
                    if [[ -n "\$UNI_TAF" ]]; then
                        tool="taf.tui\${gtag}\${ftag}\${tt}"
                        stem="\${tool}_\${sid}_\${N}"
                        rowfiles[\$stem]="\$LOGDIR/row_\${stem}.tsv"
                        ( export OMP_NUM_THREADS=\$THREADS; \\
                          run_cell "\$tool" "\$sid" "\$sci" "\$common" "\$N" "\$TIME_BUDGET" \\
                            "\$TAFFY" lift -i "\$UNI_TAF" -b "\$BED_PREFIXED" -g "\$sid" \\
                                          --maxGap "\$K" "\${ffl[@]}" \\
                                          -o "\$LOGDIR/mapped_\${stem}.bed" \\
                        ) > "\${rowfiles[\$stem]}" &
                        pids[\$stem]=\$!
                    fi
                done

                # --bin variants: only at K=0 + mode=fast (--bin is mutex
                # with --maxGap at the CLI, and requires --fast).  One
                # extra cell per (bin size, threads_per_cell).
                if [[ "\$mode" == "fast" && "\$K" == "0" && \${#BIN_SIZES[@]} -gt 0 ]]; then
                    for B in "\${BIN_SIZES[@]}"; do
                        btag="_bin\$(fmt_bin \$B)"
                        for THREADS in "\${THREADS_PER_CELL[@]}"; do
                            tt="\$(ttag \$THREADS)"
                            if [[ -n "\$UNI" ]]; then
                                tool="maf.tui\${ftag}\${btag}\${tt}"
                                stem="\${tool}_\${sid}_\${N}"
                                rowfiles[\$stem]="\$LOGDIR/row_\${stem}.tsv"
                                ( export OMP_NUM_THREADS=\$THREADS; \\
                                  run_cell "\$tool" "\$sid" "\$sci" "\$common" "\$N" "\$TIME_BUDGET" \\
                                    "\$TAFFY" lift -i "\$UNI" -b "\$BED_PREFIXED" -g "\$sid" \\
                                                  --fast --bin "\$B" \\
                                                  -o "\$LOGDIR/mapped_\${stem}.bed" \\
                                ) > "\${rowfiles[\$stem]}" &
                                pids[\$stem]=\$!
                            fi
                            if [[ -n "\$UNI_TAF" ]]; then
                                tool="taf.tui\${ftag}\${btag}\${tt}"
                                stem="\${tool}_\${sid}_\${N}"
                                rowfiles[\$stem]="\$LOGDIR/row_\${stem}.tsv"
                                ( export OMP_NUM_THREADS=\$THREADS; \\
                                  run_cell "\$tool" "\$sid" "\$sci" "\$common" "\$N" "\$TIME_BUDGET" \\
                                    "\$TAFFY" lift -i "\$UNI_TAF" -b "\$BED_PREFIXED" -g "\$sid" \\
                                                  --fast --bin "\$B" \\
                                                  -o "\$LOGDIR/mapped_\${stem}.bed" \\
                                ) > "\${rowfiles[\$stem]}" &
                                pids[\$stem]=\$!
                            fi
                        done
                    done
                fi
            done
        done
    done < "\$SPECIES_TSV"

    # Wait on all cells, append rows.
    for stem in "\${!pids[@]}"; do
        wait "\${pids[\$stem]}" || true
        [[ -s "\${rowfiles[\$stem]}" ]] && cat "\${rowfiles[\$stem]}" >> "\$BENCH_TSV"
    done
    unset pids rowfiles

    echo "=== wave N=\$N took \$((SECONDS - t_wave)) s ==="
done

echo "bench done.  TSV: \$BENCH_TSV"
EOF
chmod +x "$RUNNER"

# --- Companion plot script. -------------------------------------------
PLOT="$OUTDIR/plot.py"
cat > "$PLOT" <<'PY'
#!/usr/bin/env python3
"""Plot wall + peak RSS vs divergence-from-reference for each (tool, size)."""
import csv, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

bench_dir = os.path.dirname(os.path.abspath(__file__))

# Find the tree path: stamped into the runner script's TREE= line.
tree_path = None
with open(os.path.join(bench_dir, "bench.sh")) as f:
    for line in f:
        m = re.match(r'^TREE="([^"]+)"', line.strip())
        if m:
            tree_path = m.group(1); break
if not tree_path:
    tree_path = os.environ.get("TREE_NWK") or os.path.join(bench_dir, "tree.nwk")
if not os.path.exists(tree_path):
    sys.exit(f"ERROR: tree file not found: {tree_path}  (set TREE_NWK or copy as tree.nwk)")

# Compute leaf-to-leaf branch-length distances on the newick.
text = open(tree_path).read().strip()
if text.endswith(';'):
    text = text[:-1]
i = [0]
def parse_clade():
    children = []
    if text[i[0]] == '(':
        i[0] += 1
        while True:
            children.append(parse_clade())
            if text[i[0]] == ',':
                i[0] += 1
            elif text[i[0]] == ')':
                i[0] += 1
                break
    j = i[0]
    while j < len(text) and text[j] not in '(),:;':
        j += 1
    name = text[i[0]:j]; i[0] = j
    bl = 0.0
    if i[0] < len(text) and text[i[0]] == ':':
        i[0] += 1
        k = i[0]
        while k < len(text) and text[k] not in '(),;':
            k += 1
        bl = float(text[i[0]:k]); i[0] = k
    return (name, bl, children)
root = parse_clade()

parent = {}
node_of_name = {}
nid = [0]
def annotate(node, par=None):
    n = nid[0]; nid[0] += 1
    name, bl, ch = node
    parent[n] = (par, name, bl, ch)
    if not ch:
        node_of_name[name] = n
    for c in ch:
        annotate(c, n)
annotate(root)

def root_path(leaf_id):
    out = []; cur = leaf_id; cum = 0.0
    while True:
        out.append((cur, cum))
        par, _, bl, _ = parent[cur]
        if par is None: break
        cur = par; cum += bl
    return out

def dist(a, b):
    pa = {n: c for n, c in root_path(a)}
    for n, c in root_path(b):
        if n in pa:
            return pa[n] + c
    return None

# Load bench.tsv
rows = []
with open(os.path.join(bench_dir, "bench.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            r["size_bp"] = int(r["size_bp"])
            r["wall_s"]  = float(r["wall_s"]) if r["wall_s"] != "NA" else None
            r["peak_rss_kb"] = float(r["peak_rss_kb"]) if r["peak_rss_kb"] != "NA" else None
            r["timed_out"] = int(r["timed_out"])
            r["n_intervals"] = int(r["n_intervals"])
            rows.append(r)
        except (ValueError, KeyError):
            continue

# Resolve reference from bench.sh (REF= line).
ref = None
with open(os.path.join(bench_dir, "bench.sh")) as f:
    for line in f:
        m = re.match(r'^REF="([^"]+)"', line.strip())
        if m:
            ref = m.group(1); break
if not ref or ref not in node_of_name:
    sys.exit(f"ERROR: reference {ref!r} not found in tree leaves")
ref_id = node_of_name[ref]

# Augment each row with divergence + display label.
for r in rows:
    g = r["genome_id"]
    if g in node_of_name:
        r["dist"] = dist(ref_id, node_of_name[g])
    else:
        r["dist"] = None
rows = [r for r in rows if r["dist"] is not None]

# Distinct sizes & order species by distance.
sizes  = sorted({r["size_bp"] for r in rows})
species_order = sorted({(r["dist"], r["genome_id"], r["common_name"]) for r in rows})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
fig.subplots_adjust(left=0.08, right=0.97, top=0.90, bottom=0.20, wspace=0.27)

# Base palette per tool family.  --maxGap variants (maf.tui_gK / taf.tui_gK)
# inherit the family colour but shift hue/lightness by K so K=0 is the
# saturated base and larger K is paler / brighter.
BASE_COLOR = {
    "liftover":    "#d62728",
    "halLiftover": "#2ca02c",
    "maf.tui":     "#1f77b4",
    "taf.tui":     "#9467bd",
}
import matplotlib.colors as mcolors

def _parse_bin_suffix(s):
    # '1M' -> 1_000_000; '100k' -> 100_000; '500' -> 500; bad -> None.
    m = re.match(r'^(\d+)([GMk]?)$', s)
    if not m: return None
    return int(m.group(1)) * {'G': 10**9, 'M': 10**6, 'k': 10**3, '': 1}[m.group(2)]

def parse_tool(tool):
    """Return (family, gap, is_fast, bin, threads) for taffy variants:
    'maf.tui' / 'maf.tui_fast' / 'maf.tui_fast_bin1M_t8' / etc.
    Suffix strip order (outermost first): _t<N>, _bin<S>, _fast, _g<K>.
    Non-tui tools return (tool, None, False, None, 1)."""
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
    for fam in ("maf.tui", "taf.tui"):
        if base == fam:           return fam, 0, is_fast, bin_size, threads
        prefix = fam + "_g"
        if base.startswith(prefix):
            try: return fam, int(base[len(prefix):]), is_fast, bin_size, threads
            except ValueError: pass
    return tool, None, False, bin_size, threads

def tool_color(tool):
    fam, gap, _, _, _ = parse_tool(tool)
    base = BASE_COLOR.get(fam, "#666666")
    if gap is None or gap == 0:
        return base
    # Lighten toward white as log10(K) grows.  K=10 -> small shift,
    # K=10000 -> ~60% toward white.  Caps at 75% lightness shift.
    import math
    t = min(0.75, 0.20 * math.log10(max(gap, 1)))
    r, g, b = mcolors.to_rgb(base)
    return (r + (1 - r) * t, g + (1 - g) * t, b + (1 - b) * t)

size_marker = {sz: m for sz, m in zip(sizes, ["o", "s", "^", "D", "v", "*"])}

# Tools present in the data, ordered: non-tui first, then tui families
# sorted within by K.
def tool_sort_key(t):
    fam, gap, is_fast, bin_size, threads = parse_tool(t)
    fam_order = {"liftover": 0, "halLiftover": 1, "maf.tui": 2, "taf.tui": 3}.get(fam, 9)
    # Order within family: by gap, then default-before-fast, then no-bin
    # before binned variants (ascending bin size), then ascending thread count.
    return (fam_order, gap if gap is not None else -1, int(is_fast),
            -1 if bin_size is None else bin_size, threads)
tools = sorted({r["tool"] for r in rows}, key=tool_sort_key)

def by_tool_size(tool, size):
    # Wall-time data point = mean per-interval seconds (wall_s / n_intervals).
    out = [(r["dist"], r["wall_s"] / max(r["n_intervals"], 1), r["peak_rss_kb"])
           for r in rows
           if r["tool"] == tool and r["size_bp"] == size
           and not r["timed_out"] and r["wall_s"] is not None]
    out.sort()
    return out

for tool in tools:
    color = tool_color(tool)
    _, _, is_fast, bin_size, _ = parse_tool(tool)
    # --fast variants render with dashed line + open marker so they
    # overlay cleanly on top of their column-walk counterpart (same
    # color/marker family).  --bin variants are dotted (":") to stack
    # a third visually-distinct layer on top of (default, fast).
    linestyle = ":" if bin_size else ("--" if is_fast else "-")
    markerfacecolor = "none" if (is_fast or bin_size) else None
    for sz in sizes:
        pts = by_tool_size(tool, sz)
        if not pts: continue
        xs = [p[0] for p in pts]
        wt = [p[1] for p in pts]
        rs = [p[2] / 1024.0 for p in pts]  # MB
        label = f"{tool} N={sz:,}"
        kw = dict(marker=size_marker[sz], linestyle=linestyle, color=color, label=label)
        if markerfacecolor is not None:
            kw["markerfacecolor"] = markerfacecolor
        ax1.plot(xs, wt, **kw)
        ax2.plot(xs, rs, **kw)

# X-axis labels = species common names at their divergence positions.
xticks = [d for d, _, _ in species_order]
xlabels = [f"{c}\n({d:.2f})" for d, _, c in species_order]
for ax, title, ylab in [(ax1, "mean lift time per interval", "seconds / interval"),
                         (ax2, "peak RSS", "MB")]:
    ax.set_yscale("log")
    ax.set_xlabel(f"divergence from {ref} (subs/site)")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=8)
    ax.grid(True, which="both", axis="y", alpha=0.3)
    ax.legend(fontsize=8, loc="best")

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
    -J taffy_lift_bench
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
