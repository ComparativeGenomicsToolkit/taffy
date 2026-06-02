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
CHAINS_DIR=""
TREE=""
OUTDIR=""
T_TOTAL=24                          # cpus-per-task; one core per concurrent cell
REF="GCF_016700215.2"               # chicken
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

Required:
  -u FILE       Universal MAF (.uni.maf.gz; only its .tui sibling is read)
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

Override binary paths via env: TAFFY, LIFTOVER
EOF
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u)             UNI="$2"; shift 2;;
        -c)             CHAINS_DIR="$2"; shift 2;;
        -t)             TREE="$2"; shift 2;;
        -o)             OUTDIR="$2"; shift 2;;
        -r)             REF="$2"; shift 2;;
        -T)             T_TOTAL="$2"; shift 2;;
        -S)             SIZES="$2"; shift 2;;
        -N)             N_INTERVALS="$2"; shift 2;;
        --seed)         SEED="$2"; shift 2;;
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

for v in UNI CHAINS_DIR TREE OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: -$(echo $v | cut -c1) required" >&2; usage 1; }
done
[[ -n "$TAFFY"    ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$LIFTOVER" ]] || { echo "ERROR: liftOver not on PATH (set \$LIFTOVER)" >&2; exit 1; }
[[ -f "$UNI"      ]] || { echo "ERROR: $UNI not found" >&2; exit 1; }
[[ -f "${UNI}.tui" ]] || { echo "ERROR: $UNI has no .tui sibling" >&2; exit 1; }
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
echo ">> querying .tui for ${REF}.* chroms via taffy stats -s ..." >&2
# Don't pipe-redirect stderr -- we want any taffy errors to surface.
# `set -e` + `pipefail` will still abort the script if anything failed,
# but with a visible cause.
"$TAFFY" stats -i "$UNI" -s \
    | awk -v p="${REF}." 'index($1, p) == 1 {sub(p, "", $1); print $1"\t"$2}' \
    | sort -k2,2nr > "$REF_SIZES"
[[ -s "$REF_SIZES" ]] || {
    echo "ERROR: no sequences matching ${REF}.* in $UNI's .tui." >&2
    echo "       Check that -r matches the genome prefix used in the universal MAF." >&2
    exit 1
}
echo ">>   found $(wc -l < "$REF_SIZES") chroms matching ${REF}.*" >&2

# Quick cross-check: does at least one chain header tName match one of
# our REF_SIZES chroms?  Catches the bench-killing case where chain
# names disagree with .tui names (e.g., chr1 vs NC_006088.5) BEFORE the
# job runs.  Loop over chains since the first one may legitimately have
# 0 chrom intersection (an alt-only chrom pair) if we're unlucky.
INTERSECT_OK=0
for ch in "${CHAIN_OF[@]}"; do
    # Sample first 200 chain headers (cheap enough on a multi-GB chain.gz).
    CHAIN_TNAMES=$(zcat "$ch" | awk '/^chain/ {print $3}' | head -200 | sort -u)
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
    zcat "$(echo "${CHAIN_OF[@]}" | awk '{print $1}')" \
        | awk '/^chain/ {print "       "$3}' | head -5 >&2
    echo "      .tui ${REF}.* chroms (sample):" >&2
    head -5 "$REF_SIZES" | awk '{print "       "$1}' >&2
    echo "      liftOver cells will likely fail; taffy cells should still work." >&2
fi

# --- Generate the bed files (one per size; shared across species) -------
# Each line: <chrom>\t<start>\t<end>.  Chrom names come from the chain
# tName (bare accession; e.g. NC_006088.5).  For taffy lift we prefix
# with "<REF>." on-the-fly inside the runner.  Picks chroms weighted by
# size, then offsets uniformly inside.
IFS=',' read -r -a SIZE_ARR <<< "$SIZES"

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
echo ">> uni:           $UNI  (using \$UNI.tui only)"
echo ">> chains dir:    $CHAINS_DIR"
echo ">> tree:          $TREE"
echo ">> reference:     $REF"
echo ">> species:       $N_SPECIES"
awk -F'\t' '{printf("                  %-18s %-30s %s\n", $1, $2, $3)}' "$SPECIES_TSV"
echo ">> sizes:         ${SIZE_ARR[*]}  ($N_INTERVALS intervals each, seed=$SEED)"
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
CHAINS_DIR="$CHAINS_DIR"
OUTDIR="$OUTDIR"
REF="$REF"
TREE="$TREE"     # for plot.py: divergence-from-ref x-axis
T_TOTAL=$T_TOTAL
TIME_BUDGET=$TIME_BUDGET
STAGE_LOCAL=$STAGE_LOCAL
TAFFY="$TAFFY"
LIFTOVER="$LIFTOVER"
SPECIES_TSV="$OUTDIR/species.tsv"
SIZES=( ${SIZE_ARR[*]} )

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
    # .tui only -- taffy lift never opens the MAF itself.
    stage_one "\$UNI.tui" > /dev/null
    LOCAL_UNI="\$STAGE_DIR/\$(basename "\$UNI")"
    # taffy lift wants <input>.tui; create a stub <input> next to it so
    # it can derive the .tui path.  Stub is 0 bytes, never read.
    : > "\$LOCAL_UNI"
    UNI="\$LOCAL_UNI"

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
    printf "tool\tgenome_id\tsci_name\tcommon_name\tsize_bp\tn_intervals\twall_s\tpeak_rss_kb\texit\ttimed_out\tn_mapped\tn_unmapped\n" > "\$BENCH_TSV"
fi

# run_cell tool species_id sci common size budget cmd...
# Writes one TSV row to stdout.  n_mapped / n_unmapped come from line
# counts of the cell's output beds (kept around for sanity-checking).
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

    local n_mapped=0 n_unmapped=0
    [[ -s "\$out_bed" ]] && n_mapped=\$(wc -l < "\$out_bed")
    [[ -s "\$unm_bed" ]] && n_unmapped=\$(grep -v '^#' "\$unm_bed" | wc -l)

    printf "%s\t%s\t%s\t%s\t%d\t\$N_INTERVALS_INNER\t%s\t%s\t%d\t%d\t%d\t%d\n" \\
        "\$tool" "\$sid" "\$sci" "\$common" "\$N" "\$wall" "\$rss" "\$rc" "\$timed_out" "\$n_mapped" "\$n_unmapped"
}

# Per cell N_intervals is fixed across the run (constant from driver).
N_INTERVALS_INNER=$N_INTERVALS

# --- Outer loop: size waves; inner = species x 2 tools in parallel. ----
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
        ( run_cell liftover "\$sid" "\$sci" "\$common" "\$N" "\$TIME_BUDGET" \\
            "\$LIFTOVER" "\$BED_NATIVE" "\$chain" "\$out_lo" "\$unm_lo" \\
          ) > "\${rowfiles[\$stem_lo]}" &
        pids[\$stem_lo]=\$!

        # ---- taffy lift cell ----
        stem_tf="taffy_\${sid}_\${N}"
        rowfiles[\$stem_tf]="\$LOGDIR/row_\${stem_tf}.tsv"
        out_tf="\$LOGDIR/mapped_\${stem_tf}.bed"
        ( run_cell taffy "\$sid" "\$sci" "\$common" "\$N" "\$TIME_BUDGET" \\
            "\$TAFFY" lift -i "\$UNI" -b "\$BED_PREFIXED" -g "\$sid" -o "\$out_tf" \\
          ) > "\${rowfiles[\$stem_tf]}" &
        pids[\$stem_tf]=\$!
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

tool_color = {"liftover": "#d62728", "taffy": "#1f77b4"}
size_marker = {sz: m for sz, m in zip(sizes, ["o", "s", "^", "D", "v", "*"])}

def by_tool_size(tool, size):
    # Wall-time data point = mean per-interval seconds (wall_s / n_intervals).
    # Lets the three size lines be compared directly across an axis where the
    # cell-totals would otherwise be dominated by N.
    out = [(r["dist"], r["wall_s"] / max(r["n_intervals"], 1), r["peak_rss_kb"])
           for r in rows
           if r["tool"] == tool and r["size_bp"] == size
           and not r["timed_out"] and r["wall_s"] is not None]
    out.sort()
    return out

for tool in ("liftover", "taffy"):
    for sz in sizes:
        pts = by_tool_size(tool, sz)
        if not pts: continue
        xs = [p[0] for p in pts]
        wt = [p[1] for p in pts]
        rs = [p[2] / 1024.0 for p in pts]  # MB
        label = f"{tool} N={sz:,}"
        ax1.plot(xs, wt, marker=size_marker[sz], linestyle="-", color=tool_color[tool], label=label)
        ax2.plot(xs, rs, marker=size_marker[sz], linestyle="-", color=tool_color[tool], label=label)

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
