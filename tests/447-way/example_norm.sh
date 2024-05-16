#!/usr/bin/env zsh

### As an example of how to use the taffy commands to make different kinds of
### normalized alignment files suitable for use by the taffy.ml iterators.

# Get the taffy root directory from which this script will run
taffy_root=$1

# exit when any command fails
set -e

# log the commands its running
set -x

# The prefix for the alignment output files
alignment_file_prefix=${taffy_root}/tests/447-way/447-mammalian-2022v1_chr22_22000000_22100000

## Combinations of output:
# (1) No masking bases
no_masking_alignment_file=${alignment_file_prefix}.no_masking_with_ancestors.taf.gz
echo "No masking of bases file " ${no_masking_alignment_file}
# (2) Reference masking of bases
reference_masking_with_ancestors_alignment_file=${alignment_file_prefix}.reference_masking_with_ancestors.taf.gz
echo "Reference masking with ancestors file " ${reference_masking_with_ancestors_alignment_file}
# (3) Lineage masking of bases
lineage_masking_alignment_file=${alignment_file_prefix}.lineage_masking.taf.gz
echo "Lineage masking file ", ${lineage_masking_alignment_file}
# (4) No masking bases, no ancestors
no_masking_no_ancestors_alignment_file=${alignment_file_prefix}.no_masking_no_ancestors.taf.gz
echo "No masking, no ancestors file " ${no_masking_no_ancestors_alignment_file}
# (5) Reference masking of bases, no ancestors
reference_masking_no_ancestors_alignment_file=${alignment_file_prefix}.reference_masking_no_ancestors.taf.gz
echo "Reference masking, no ancestors file", ${reference_masking_no_ancestors_alignment_file}

# The tree file
tree_file=${taffy_root}/tests/447-way/447-mammalian-2022v1.nh

# Rerooted tree file
rerooted_tree_file=${taffy_root}/tests/447-way/447-mammalian-2022v1.rerooted.nh

# Wiggle file
wig_file=${taffy_root}/tests/447-way/447-mammalian-2022v1_hg38_chr22_22000000_22100000.phyloP.wig

# Raw alignment file
alignment_file=${taffy_root}/tests/447-way/447-mammalian-2022v1_hg38_chr22_22000000_22100000.anc.norm.taf.gz

# Sort/dup-filter file
sort_file=${taffy_root}/tests/447-way/447-sort.txt

# Filter file
filter_file=${taffy_root}/tests/447-way/447-filter.txt

# Ref
ref=hg38

# Reroot the tree
time ${taffy_root}/scripts/manipulate_tree.py --reroot $ref  --out_file $rerooted_tree_file $tree_file

# Make the sort file
time ${taffy_root}/scripts/tree_to_sort_file.py --traversal pre --reroot $ref --out_file $sort_file --suffix_to_append_to_labels . $tree_file

# Make the filter file
time ${taffy_root}/scripts/tree_to_sort_file.py --out_file $filter_file --no_leaf_nodes --suffix_to_append_to_labels . $tree_file

# Sort/dedup/filter the alignment to create a temporary alignment file in which each species is present exactly once, then
# apply taffy annotate to add the phyloP tags
# This creates the no masking alignment file
time ${taffy_root}/bin/taffy sort -i ${alignment_file} -n ${sort_file} -p ${sort_file} -d ${sort_file} | ${taffy_root}/bin/taffy annotate -w ${wig_file} --tagName phyloP --refPrefix hg38. -c | ${taffy_root}/bin/taffy view -c --runLengthEncodeBases > ${no_masking_alignment_file}
time ${taffy_root}/bin/taffy index -i ${no_masking_alignment_file}

# Now build the additional final alignment files and indexes

# (2) Reference masking of bases
time ${taffy_root}/bin/taffy view -i ${no_masking_alignment_file} -a -c > ${reference_masking_with_ancestors_alignment_file}
time ${taffy_root}/bin/taffy index -i ${reference_masking_with_ancestors_alignment_file}

# (3) Lineage masking of bases
time ${taffy_root}/bin/taffy view -i ${no_masking_alignment_file} -t $rerooted_tree_file -b -c > ${lineage_masking_alignment_file}
time ${taffy_root}/bin/taffy index -i ${lineage_masking_alignment_file}

# (4) No masking bases, no ancestors
time ${taffy_root}/bin/taffy sort -i ${no_masking_alignment_file} -f ${filter_file} -c > ${no_masking_no_ancestors_alignment_file}
time ${taffy_root}/bin/taffy index -i ${no_masking_no_ancestors_alignment_file}

# (5) Reference masking of bases, no ancestors
time ${taffy_root}/bin/taffy sort -i ${no_masking_alignment_file} -f ${filter_file} -c | ${taffy_root}/bin/taffy view -a -c > ${reference_masking_no_ancestors_alignment_file}
time ${taffy_root}/bin/taffy index -i ${reference_masking_no_ancestors_alignment_file}

# Print stats about the starting and final alignment as a sanity check

echo "Starting alignment"
time ${taffy_root}/bin/taffy stats -i ${alignment_file} -a

echo "Stats of final alignments"
time ${taffy_root}/bin/taffy stats -i ${no_masking_alignment_file} -a

