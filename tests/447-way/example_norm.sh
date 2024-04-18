#!/usr/bin/env zsh

# Run this from the taffy/tests dir for it to function correctly

# exit when any command fails
set -e

# log the commands its running
set -x

# The tree file
tree_file=./447-way/447-mammalian-2022v1.nh

# Raw alignment file
alignment_file=./447-way/447-mammalian-2022v1_chr22_22000000_22100000.taf.gz

# Normalized taf file
norm_alignment_file=./447-way/447-mammalian-2022v1_chr22_22000000_22100000.norm.taf.gz

# Sorted/deduped/filters taf file
final_alignment_file=./447-way/447-mammalian-2022v1_chr22_22000000_22100000.final.taf.gz

# Sort/dup-filter file
sort_file=./447-way/447-sort.txt

# Filter file
filter_file=./447-way/447-filter.txt

# Ref
ref=hg38

# Make the sort file
time ../scripts/tree_to_sort_file.py --traversal pre --reroot $ref --out_file $sort_file --no_internal_nodes --no_check_labels $tree_file

# Make the filter file
time ../scripts/tree_to_sort_file.py --out_file $filter_file --no_leaf_nodes --no_check_labels $tree_file

# Normalize the alignment
time ../bin/taffy norm -i $alignment_file -c > $norm_alignment_file

# Now sort/dedup/filter the alignment
time ../bin/taffy sort -i $norm_alignment_file -n $sort_file -p $sort_file -d $sort_file -f $filter_file  | ../bin/taffy view -c > $final_alignment_file

