
class PhyloTree:
    def __init__(self, label='', children=None, branch_length=0.0, parent=None):
        """ Represents a phylogenetic tree.

        :param label: Label of node, should be empty string if not specified
        :param children: List of PhyloTree child nodes
        :param branch_length: Length of branch of node
        :param parent: Reference to parent node
        """
        self.label = label
        self.children = [] if not children else children[:]
        self.branch_length = branch_length
        self.parent = parent

    def tree_string(self, include_branch_lengths=True, include_internal_labels=True, include_leaf_labels=True,
                    branch_length_format=(lambda x: f'{x}')):
        """ Create a string representing the tree

        :param include_branch_lengths: Include branch lengths
        :param include_internal_labels:  Include labels of internal nodes
        :param include_leaf_labels: Include labels of leaf nodes
        :param branch_length_format: Formatting function for branch lengths
        :return:
        """
        def f(child):
            return child.tree_string(include_branch_lengths,include_internal_labels,
                                     include_leaf_labels,branch_length_format)
        subtree_string = f'({",".join([f(i) for i in self.children])})' if len(self.children) > 0 else ''
        if len(self.children) == 0:
            label_string = self.label if include_leaf_labels else ""
        else:
            label_string = self.label if include_internal_labels else ""
        branch_length_string = f':{branch_length_format(self.branch_length)}' if include_branch_lengths else ''
        return f'{subtree_string}{label_string}{branch_length_string}{"" if self.parent else ";"}'

    @staticmethod
    def newick_tree_parser(newick_tree, default_distance=1.0, default_label=''):
        """ Lax newick tree parser.

        :param newick_tree: String representing phylogenetic tree in newick tree format
        :param default_distance: Default branch length to set for any branch not specified
        :param default_label: Default branch label if not specified
        :return:
        """
        # Make it so that we can tokenize the grammer symbols
        newick_tree = newick_tree.replace("(", " ( ")
        newick_tree = newick_tree.replace(")", " ) ")
        newick_tree = newick_tree.replace(":", " : ")
        newick_tree = newick_tree.replace(";", "")
        newick_tree = newick_tree.replace(",", " , ")

        tree_tokens = [i for i in newick_tree.split() if i != '']  # Get all tokens, minus any empty strings

        def get_branch_length(i):
            if i < len(tree_tokens) and tree_tokens[i] == ':':
                if i+1 == len(tree_tokens):
                    raise ValueError(f"Tree string incomplete parsing label: {''.join(tree_tokens)}")
                return float(tree_tokens[i+1]), i+2
            else:
                return default_distance, i

        def get_label(i):
            if i < len(tree_tokens):
                j = tree_tokens[i]
                if j != ':' and j != ')' and j != ',':
                    return j, i+1
            return default_label, i

        def traverse_tree(i):  # Recursive method to parse the tree
            children = []
            if tree_tokens[i] == '(':  # Has children
                i += 1
                if i == len(tree_tokens):
                    raise ValueError(f"Tree string incomplete: {''.join(tree_tokens)}")
                while tree_tokens[i] != ')':
                    subtree, i = traverse_tree(i)
                    children.append(subtree)
                    if i == len(tree_tokens):
                        raise ValueError(f"Tree string incomplete: {''.join(tree_tokens)}")
                    if tree_tokens[i] == ',':
                        i += 1
                    if i == len(tree_tokens):
                        raise ValueError(f"Tree string incomplete: {''.join(tree_tokens)}")
                i += 1
            # Get any label and branch lengths
            label, i = get_label(i)
            branch_length, i = get_branch_length(i)

            # Build the tree
            tree = PhyloTree(children=children, label=label, branch_length=branch_length)

            # Set the parent links
            for child in children:
                child.parent = tree

            # Done
            return tree, i

        return traverse_tree(0)[0]

    def get_descendant(self, node_label):
        """ Get a node that is a descendant of this node and has the given node_label. Otherwise None.
        If multiple nodes with label exist will return an arbitrary one.

        :param node_label: Label of node to search for
        :return: PhyloTree node if found, else None
        """
        if self.label == node_label:
            return self
        for node in self.children:
            n = node.get_descendant(node_label)
            if n is not None:
                return n
        return None

    def reroot(self):
        """ Reroot the tree so that this node is the root of the tree.
        """
        # First, if we are not already at the root
        if self.parent:
            self.parent.reroot()  # Recursively re-root the parent node
            self.children.append(self.parent)  # Make the parent a child of this node
            assert self.parent.parent is None
            self.parent.children.remove(self)  # Remove this node from the parent's list of children
            self.parent.parent = self  # Make sure the parent's parent node is now pointing at this node
            self.parent = None  # Set the parent node to None

    def __str__(self):
        return self.tree_string()

