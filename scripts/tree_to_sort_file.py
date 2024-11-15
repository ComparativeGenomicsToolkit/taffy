#!python3

""" Script to create a sort file for taffy sort from a phylogenetic tree
"""

from taffy.newick import PhyloTree
import argparse


def traverse_tree(phylo_tree, traversal_type, include_internal_labels=True, include_leaf_labels=True,
                  suffix_to_append_to_labels=""):
    """ Gets a list of the labels in a tree in order of a traversal of the tree.
    """
    assert traversal_type in ("pre", "mid", "post")
    ordered_labels = []

    def traverse_tree_recursive(t):
        # The root, if it has single child is considered a leaf - this means it won't be dropped
        # if we exclude internal labels and won't be included if we exclude internal nodes
        include_label = (include_internal_labels and
                         ((t.parent is not None and len(t.children) > 0) or len(t.children) > 1)) or \
                        (include_leaf_labels and (len(t.children) == 0 or (len(t.children) == 1 and t.parent is None)))
        if traversal_type == "pre" and include_label:
            ordered_labels.append(t.label + suffix_to_append_to_labels)
        mid_point = len(t.children) // 2
        for i in range(mid_point):
            traverse_tree_recursive(t.children[i])
        if traversal_type == "mid" and include_label:
            ordered_labels.append(t.label + suffix_to_append_to_labels)
        for i in range(mid_point, len(t.children)):
            traverse_tree_recursive(t.children[i])
        if traversal_type == "post" and include_label:
            ordered_labels.append(t.label + suffix_to_append_to_labels)

    traverse_tree_recursive(phylo_tree)
    return ordered_labels


def main():
    # Set script arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "tree",
        type=str,
        help="Path of newick tree file"
    )
    parser.add_argument(
        "--traversal",
        type=str,
        default="mid",
        help="Either 'pre', 'mid' or 'post'. Specifies how the tree is traversed to make the sort file ordering."
    )
    parser.add_argument(
        "--no_internal_nodes",
        dest="internal_labels",
        action='store_false',
        default=True,
        help="Toggle off internal node labels"
    )
    parser.add_argument(
        "--no_leaf_nodes",
        dest="leaf_labels",
        action='store_false',
        default=True,
        help="Toggle off leaf node labels"
    )
    parser.add_argument(
        "--reroot",
        type=str,
        default=None,
        help="Reroot the tree so that the given node is the root."
    )
    parser.add_argument(
        "--no_check_labels",
        dest="check_labels",
        action='store_false',
        default=True,
        help="Don't check that no label is a prefix of another label."
    )
    parser.add_argument(
        "--suffix_to_append_to_labels",
        type=str,
        default="",
        help="Append this string to every label to prevent one label being a prefix of another."
    )
    parser.add_argument(
        "--out_file",
        type=str,
        default=None,
        help="Write to the given file instead of standard out."
    )

    # Parse arguments
    args = parser.parse_args()

    # Parse the newick tree
    with open(args.tree) as fh:
        tree = PhyloTree.newick_tree_parser(fh.read())

    # Re-root if needed
    if args.reroot:
        tree = tree.get_descendant(args.reroot)
        if tree is None:
            raise RuntimeError(f"Couldn't fine {args.reroot} to reroot the tree")
        tree.reroot()

    print("wooooo", args.suffix_to_append_to_labels)

    # Get the labels
    ordered_labels = traverse_tree(tree, args.traversal, include_internal_labels=args.internal_labels,
                                   include_leaf_labels=args.leaf_labels,
                                   suffix_to_append_to_labels=args.suffix_to_append_to_labels)

    # Check the labels in the tree
    if args.check_labels:
        for i in range(len(ordered_labels)):
            label = ordered_labels[i]
            for j in range(i+1, len(ordered_labels)):
                label2 = ordered_labels[j]
                if label2.startswith(label):
                    raise RuntimeError(f"Node label: {label} is a prefix of label {label2}")

    # Finally, print the labels
    if args.out_file:  # To a given file
        with open(args.out_file, "w") as fh:
            for label in ordered_labels:
                fh.write(f"{label}\n")
    else:  # Or to standard out
        for label in ordered_labels:
            print(label)


if __name__ == '__main__':
    main()
