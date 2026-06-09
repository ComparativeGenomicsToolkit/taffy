#!/bin/bash
#
# taffy lift concordance + timing benchmark -- LOCAL, whole-genome.
#
# Lifts REF's whole genome to TGT several ways and measures (a) how well
# they agree -- Jaccard of TGT-genome coverage -- and (b) how fast they are:
#
#   taffy lift   on the base .tui and each chained sidecar (--fast mode)
#   UCSC liftOver  with a precomputed pairwise <REF>_vs_<TGT>.chain.gz
#
# Every method gets the SAME whole-genome input (one 0..chromLen interval
# per REF chrom).  Each output is reduced to sorted+merged BED3 in TGT
# coords; the report is per-method wall / peak-RSS / interval count / TGT
# coverage bp, plus a pairwise bedtools-jaccard concordance matrix.
#
# Defaults reproduce the chicken (GCF_016700215.2) -> quail (GCA_039878825.1)
# lift on the vgp-577way universal alignment.
#
# Sidecar note: `taffy tui-chain`'s default --maxGap is 1000, so the
# `<uni>.chained.tui` sidecar IS the g1000 (axtChain-equivalent) index;
# `<uni>.chained_g10000.tui` is the more aggressive g10000.  A sidecar is
# selected by the SUFFIX appended to taffy lift's -i (it reads <i>.tui):
# "" = base .tui, ".chained" = g1000, ".chained_g10000" = g10000.
#
# Usage:
#   taffy_lift_concordance_local.sh [-u UNI.maf.gz] [-r REF] [-g TGT]
#       [-c CHAIN.chain.gz] [-o OUTDIR] [-S label:suffix,label:suffix,...]
#   Binaries via env: TAFFY, LIFTOVER, BEDTOOLS.
#
# Re-run the default (chicken->quail) bench:   ./taffy_lift_concordance_local.sh
#
set -uo pipefail

UNI=/home/hickey/dev/work/unitaf/vgp-577way-v1.uni.maf.gz       # reads <UNI>.tui (+ sidecars)
REF=GCF_016700215.2                                            # chicken (universal row-0)
TGT=GCA_039878825.1                                            # quail
CHAIN=/home/hickey/dev/work/unitaf/GCF_016700215.2_vs_GCA_039878825.1.chain.gz
OUTDIR=/tmp/cq_bench
# Sidecars to sweep, as "label:suffix" (suffix appended to UNI for -i).
# "" = base .tui; ".chained" = g1000; ".chained_g10000" = g10000.
SIDECARS_CSV="base:,g1000:.chained,g10000:.chained_g10000"
TAFFY=${TAFFY:-/home/hickey/dev/taffy/bin/taffy}
LIFTOVER=${LIFTOVER:-$(command -v liftOver || true)}
BEDTOOLS=${BEDTOOLS:-$(command -v bedtools || true)}

usage() {
    sed -n '2,/^set -uo/p' "$0" | sed 's/^# \{0,1\}//; /^set -uo/d'
    exit "${1:-0}"
}
while getopts "u:r:g:c:o:S:h" opt; do case "$opt" in
    u) UNI=$OPTARG;;  r) REF=$OPTARG;;  g) TGT=$OPTARG;;  c) CHAIN=$OPTARG;;
    o) OUTDIR=$OPTARG;;  S) SIDECARS_CSV=$OPTARG;;  h) usage 0;;  *) usage 1;;
esac; done

[[ -f "$UNI.tui" ]] || { echo "ERROR: $UNI.tui not found" >&2; exit 1; }
[[ -f "$CHAIN"   ]] || { echo "ERROR: chain not found: $CHAIN" >&2; exit 1; }
command -v "$TAFFY" >/dev/null 2>&1 || [[ -x "$TAFFY" ]] || { echo "ERROR: taffy not found: $TAFFY" >&2; exit 1; }
[[ -n "$LIFTOVER" ]] || { echo "ERROR: liftOver not on PATH (set \$LIFTOVER)" >&2; exit 1; }
[[ -n "$BEDTOOLS" ]] || { echo "ERROR: bedtools not on PATH (set \$BEDTOOLS)" >&2; exit 1; }
mkdir -p "$OUTDIR"

# Whole-genome REF bed from the chain's tName/tSize (one 0..len line/chrom).
# NATIVE names feed liftOver; REF.-prefixed names feed taffy lift (row-0 coords).
zcat "$CHAIN" | awk '$1=="chain"{print $3"\t0\t"$4}' | sort -u > "$OUTDIR/ref.native.bed"
awk -v p="$REF." '{print p$1"\t"$2"\t"$3}' "$OUTDIR/ref.native.bed" > "$OUTDIR/ref.prefixed.bed"
echo ">> $REF -> $TGT : $(wc -l < "$OUTDIR/ref.native.bed") chroms, $(awk '{s+=$3-$2}END{print s}' "$OUTDIR/ref.native.bed") bp input"

run() {  # $1=name  $2..=command (timed)
    local name=$1; shift
    echo ">> [$name] $*" >&2
    /usr/bin/time -v "$@" > "$OUTDIR/$name.out" 2> "$OUTDIR/$name.time"; local rc=$?
    local wall rss
    wall=$(grep -oE "Elapsed .*: .*" "$OUTDIR/$name.time" | sed 's/.*: //')
    rss=$(grep "Maximum resident" "$OUTDIR/$name.time" | grep -oE "[0-9]+")
    printf "%s\t%s\t%s\t%s\n" "$name" "${wall:-NA}" "${rss:-NA}" "$rc" >> "$OUTDIR/timings.tsv"
}

printf "method\twall\trss_kb\trc\n" > "$OUTDIR/timings.tsv"
METHODS=()
IFS=, read -ra SC <<< "$SIDECARS_CSV"
for entry in "${SC[@]}"; do
    label=${entry%%:*}; sfx=${entry#*:}
    run "$label" "$TAFFY" lift -i "$UNI$sfx" -b "$OUTDIR/ref.prefixed.bed" -g "$TGT" --fast -o "$OUTDIR/$label.bed"
    METHODS+=("$label")
done
run liftover "$LIFTOVER" -minMatch=0 -multiple "$OUTDIR/ref.native.bed" "$CHAIN" "$OUTDIR/liftover.bed" "$OUTDIR/liftover.unmapped"
METHODS+=(liftover)

# Reduce each output to sorted+merged BED3 in TGT coords; report coverage.
printf "method\trows\tmerged\ttgt_cov_bp\n" > "$OUTDIR/coverage.tsv"
for m in "${METHODS[@]}"; do
    [[ -s "$OUTDIR/$m.bed" ]] || { echo "  $m: EMPTY output (lift failed -- see $OUTDIR/$m.time)"; continue; }
    cut -f1-3 "$OUTDIR/$m.bed" | sort -k1,1 -k2,2n -S 2G | "$BEDTOOLS" merge > "$OUTDIR/$m.m.bed"
    printf "%s\t%s\t%s\t%s\n" "$m" "$(wc -l < "$OUTDIR/$m.bed")" "$(wc -l < "$OUTDIR/$m.m.bed")" \
        "$(awk '{s+=$3-$2}END{print s}' "$OUTDIR/$m.m.bed")" >> "$OUTDIR/coverage.tsv"
done

# Pairwise Jaccard of TGT coverage (1.0 = identical coverage).
printf "a\tb\tintersection\tunion\tjaccard\tn\n" > "$OUTDIR/jaccard.tsv"
for a in "${METHODS[@]}"; do for b in "${METHODS[@]}"; do
    [[ "$a" < "$b" ]] || continue
    [[ -s "$OUTDIR/$a.m.bed" && -s "$OUTDIR/$b.m.bed" ]] || continue
    line=$("$BEDTOOLS" jaccard -a "$OUTDIR/$a.m.bed" -b "$OUTDIR/$b.m.bed" 2>/dev/null | tail -1)
    printf "%s\t%s\t%s\n" "$a" "$b" "$line" >> "$OUTDIR/jaccard.tsv"
done; done

echo ">> DONE.  results in $OUTDIR/{timings,coverage,jaccard}.tsv"
echo; column -t "$OUTDIR/timings.tsv"; echo; column -t "$OUTDIR/coverage.tsv"
echo; echo "Jaccard (cols: a b intersection union jaccard n):"; column -t "$OUTDIR/jaccard.tsv"
