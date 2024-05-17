#!python3

""" Script to re-root a tree so that a chosen node is the root.
"""

from taffy.newick import PhyloTree
import argparse


def main():
    # Set script arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "tree",
        type=str,
        help="Path of newick tree file"
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
        "--no_branch_lengths",
        dest="branch_lengths",
        action='store_false',
        default=True,
        help="Toggle off branch lengths"
    )
    parser.add_argument(
        "--reroot",
        type=str,
        default=None,
        help="Reroot the tree so that the given node is the root."
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

    # Finally, print the tree
    tree_string = tree.tree_string(args.branch_lengths, args.internal_labels, args.leaf_labels)
    if args.out_file:  # To a given file
        with open(args.out_file, "w") as fh:
            fh.write(f"{tree_string}\n")
    else:
        print(f"{tree_string}\n")


if __name__ == '__main__':
    main()
