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
HAL=""
TREE=""
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

Required:
  -u FILE       Universal MAF (.uni.maf.gz; only its .tui sibling is read by taffy)
  -H FILE       HAL file (for halLiftover)
  -t FILE       Species tree (.nwk) -- used by the plot script for divergence
  -o DIR        Output directory

Optional:
  -r ID         Reference genome ID (default $REF; must match the genome
                prefix in the universal MAF / HAL)
  -T INT        cpus-per-task (default $T_TOTAL; needs ≥ 2 × n_species
                for parallel cells)
  -L FILE       Override species panel.  Format: 3 whitespace-separated
                cols per line: <genome_id> <scientific> <english>.
                '#' = comment.
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

for v in UNI HAL TREE OUTDIR; do
    [[ -n "${!v}" ]] || { echo "ERROR: -$(echo $v | cut -c1) required" >&2; usage 1; }
done
[[ -n "$TAFFY"       ]] || { echo "ERROR: taffy not on PATH (set \$TAFFY)" >&2; exit 1; }
[[ -n "$HALLIFTOVER" ]] || { echo "ERROR: halLiftover not on PATH (set \$HALLIFTOVER)" >&2; exit 1; }
[[ -f "$UNI"         ]] || { echo "ERROR: $UNI not found" >&2; exit 1; }
[[ -f "${UNI}.tui"   ]] || { echo "ERROR: $UNI has no .tui sibling" >&2; exit 1; }
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
echo ">> querying .tui for ${REF}.* chroms via taffy stats -s ..." >&2
NATIVE_BED="$OUTDIR/beds/genome.native.bed"
PREFIXED_BED="$OUTDIR/beds/genome.prefixed.bed"
"$TAFFY" stats -i "$UNI" -s \
    | awk -v p="${REF}." 'index($1, p) == 1 {
        sub(p, "", $1);
        printf "%s\t0\t%d\n", $1, $2;
      }' \
    | sort -k1,1 > "$NATIVE_BED"
[[ -s "$NATIVE_BED" ]] || {
    echo "ERROR: no sequences matching ${REF}.* in $UNI's .tui." >&2
    exit 1
}
# Prefixed bed for taffy lift (its --bed expects "<genome>.<chrom>"):
awk -v p="${REF}." '{print p$0}' "$NATIVE_BED" > "$PREFIXED_BED"

N_CHROMS=$(wc -l < "$NATIVE_BED")
REF_BP=$(awk '{s += $3 - $2} END {print s}' "$NATIVE_BED")
echo ">>   ${REF} input bed: $N_CHROMS chroms, $REF_BP bp total" >&2

echo ">> output dir:    $OUTDIR"
echo ">> taffy:         $TAFFY"
echo ">> halLiftover:   $HALLIFTOVER"
echo ">> uni:           $UNI  (using \$UNI.tui only)"
echo ">> hal:           $HAL"
echo ">> tree:          $TREE"
echo ">> reference:     $REF"
echo ">> species:       $N_SPECIES"
awk -F'\t' '{printf("                  %-18s %-32s %s\n", $1, $2, $3)}' "$SPECIES_TSV"
echo ">> cpus/task:     $T_TOTAL  (need ≥ $((N_SPECIES * 2)) for full parallel)"
(( T_TOTAL >= N_SPECIES * 2 )) || \
    echo ">> WARN: T_TOTAL=$T_TOTAL < 2 × N_SPECIES=$((N_SPECIES * 2)); cells will queue" >&2
echo ">> time budget:   $TIME_BUDGET s per cell"
echo ">> local-stage:   $([[ $STAGE_LOCAL -eq 1 ]] && echo "ON (copies to \$TMPDIR)" || echo "OFF (reads from network)")"
[[ -n "$TMP_GB" ]] && echo ">> --tmp request: ${TMP_GB} GB per task"

# --- Generate the runner script (the thing sbatch executes). -----------
RUNNER="$OUTDIR/bench.sh"
cat > "$RUNNER" <<EOF
#!/bin/bash
set -uo pipefail

UNI="$UNI"
HAL="$HAL"
OUTDIR="$OUTDIR"
REF="$REF"
TREE="$TREE"     # for plot.py: divergence-from-ref x-axis
TIME_BUDGET=$TIME_BUDGET
STAGE_LOCAL=$STAGE_LOCAL
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
# across the 18 cells that follow.
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
    # .tui (taffy lift never opens the MAF itself)
    stage_one "\$UNI.tui" > /dev/null
    LOCAL_UNI="\$STAGE_DIR/\$(basename "\$UNI")"
    : > "\$LOCAL_UNI"   # stub <input> so taffy lift's tui_path() resolves
    UNI="\$LOCAL_UNI"

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
        read -r n_mapped n_mapped_bp < <(awk '{n++; bp += \$3 - \$2} END {print n+0, bp+0}' "\$out_bed")
    fi

    printf "%s\t%s\t%s\t%s\t\$N_CHROMS_INNER\t\$REF_BP_INNER\t%s\t%s\t%d\t%d\t%d\t%d\n" \\
        "\$tool" "\$sid" "\$sci" "\$common" "\$wall" "\$rss" "\$rc" "\$timed_out" "\$n_mapped" "\$n_mapped_bp"
}

# --- Fire all (species x 2 tools) cells in parallel.  No waves. --------
echo "=== launching \$(wc -l < "\$SPECIES_TSV") species x 2 tools in parallel ==="
t0=\$SECONDS
declare -A pids rowfiles

while IFS=\$'\t' read -r sid sci common; do
    [[ -n "\$sid" ]] || continue

    # ---- halLiftover cell ----
    # halLiftover HAL SRC_GENOME SRC_BED TGT_GENOME TGT_BED
    # Native (un-prefixed) bed -- halLiftover wants bare chrom names.
    stem_hl="halLiftover_\${sid}"
    rowfiles[\$stem_hl]="\$LOGDIR/row_\${stem_hl}.tsv"
    ( run_cell halLiftover "\$sid" "\$sci" "\$common" \\
        "\$HALLIFTOVER" "\$HAL" "\$REF" "\$NATIVE_BED" "\$sid" \\
        "\$LOGDIR/mapped_\${stem_hl}.bed" \\
      ) > "\${rowfiles[\$stem_hl]}" &
    pids[\$stem_hl]=\$!

    # ---- taffy lift cell ----
    # Prefixed bed -- taffy lift -b expects "<genome>.<chrom>".
    stem_tf="taffy_\${sid}"
    rowfiles[\$stem_tf]="\$LOGDIR/row_\${stem_tf}.tsv"
    ( run_cell taffy "\$sid" "\$sci" "\$common" \\
        "\$TAFFY" lift -i "\$UNI" -b "\$PREFIXED_BED" -g "\$sid" \\
        -o "\$LOGDIR/mapped_\${stem_tf}.bed" \\
      ) > "\${rowfiles[\$stem_tf]}" &
    pids[\$stem_tf]=\$!
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

tool_color = {"halLiftover": "#2ca02c", "taffy": "#1f77b4"}
tools = ["halLiftover", "taffy"]

xs = np.arange(len(species_order))
width = 0.4

for ti, tool in enumerate(tools):
    by_rank = {rank[r["genome_id"]]: r for r in rows
               if r["tool"] == tool and not r["timed_out"] and r["wall_s"] is not None}
    wt = [by_rank.get(i, {}).get("wall_s") for i in xs]
    rs = [by_rank.get(i, {}).get("peak_rss_kb") for i in xs]
    bp = [by_rank.get(i, {}).get("n_mapped_bp", 0) for i in xs]
    wt = [v if v is not None else 0 for v in wt]
    rs = [(v/1024.0) if v is not None else 0 for v in rs]
    off = (ti - 0.5) * width
    ax_wall.bar(xs + off, wt, width=width, color=tool_color[tool], label=tool)
    ax_rss.bar (xs + off, rs, width=width, color=tool_color[tool], label=tool)
    ax_bp.bar  (xs + off, bp, width=width, color=tool_color[tool], label=tool)

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
