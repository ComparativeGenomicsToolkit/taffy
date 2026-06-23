# Merge `taffy depth --coords named --bin` partials into a mean-depth bedGraph.
#
# Input  (must be `LC_ALL=C sort -k1,1 -k2,2n` first): per-(name,bin) records
#          name <TAB> start <TAB> end <TAB> sum <TAB> cnt
#        where a single (name,bin) may appear in several records when its
#        columns straddle SLURM shard boundaries.
# Output: name <TAB> start <TAB> end <TAB> mean   (one row per (name, start/N)),
#        summing (sum,cnt) across the split records and dropping all-zero bins.
# Usage:  ... | awk -v N=<bin width in bp> -f uni_depth_merge.awk
BEGIN { OFS = "\t" }
{
  b = int($2 / N)
  if ($1 == pn && b == pb) {
    if ($2 < ps) ps = $2
    if ($3 > pe) pe = $3
    s += $4; c += $5
  } else {
    if (pn != "" && s > 0) print pn, ps, pe, s / c
    pn = $1; pb = b; ps = $2; pe = $3; s = $4; c = $5
  }
}
END { if (pn != "" && s > 0) print pn, ps, pe, s / c }
