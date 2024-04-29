#!/usr/bin/env zsh

# Get the taffy root directory from which this script will run
taffy_root=$1

# exit when any command fails
set -e

# log the commands its running
set -x

# The tree file
tree_file=${taffy_root}/tests/447-way/447-mammalian-2022v1.nh

# Rerooted tree file
rerooted_tree_file=${taffy_root}/tests/447-way/447-mammalian-2022v1.rerooted.nh

# Wiggle file
wig_file=${taffy_root}/tests/447-way/447-mammalian-2022v1_hg38_chr22_22000000_22100000.phyloP.wig

# Raw alignment file
alignment_file=${taffy_root}/tests/447-way/447-mammalian-2022v1_hg38_chr22_22000000_22100000.anc.norm.taf.gz

# Sorted/deduped/filtered taf file
rearranged_alignment_file=${taffy_root}/tests/447-way/447-mammalian-2022v1_chr22_22000000_22100000.rearranged.taf.gz

# Sorted/deduped/filtered/annotated taf file
final_alignment_file=${taffy_root}/tests/447-way/447-mammalian-2022v1_chr22_22000000_22100000.final.taf.gz

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

# Now sort/dedup/filter the alignment, using the -b option to taffy view will also only show lineage differences
time ${taffy_root}/bin/taffy sort -i $alignment_file -n $sort_file -p $sort_file -d $sort_file --logLevel DEBUG | ../bin/taffy view -t $rerooted_tree_file -b -c --runLengthEncodeBases > $rearranged_alignment_file

# Add annotations
time ${taffy_root}/bin/taffy annotate -i $rearranged_alignment_file -w $wig_file --tagName phyloP --refPrefix hg38. -c > $final_alignment_file

# Print stats about the starting and final alignment as a sanity check

echo "Starting alignment"
time ${taffy_root}/bin/taffy stats -i $alignment_file -a

echo "Final alignment"
time ${taffy_root}/bin/taffy stats -i $final_alignment_file -a

