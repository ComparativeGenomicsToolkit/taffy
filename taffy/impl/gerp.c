/*
 * GERP RS conservation scoring.  See gerp.h for the contract; this file
 * holds the tree state + Fitch parsimony + branch-sum + per-block row
 * resolution.
 */

#include "gerp.h"
#include "sonLibTree.h"
#include <ctype.h>
#include <string.h>

/////////////////////////////////////////////////////////////////////////////
// Tree state
//
// Public type GerpTree is opaque; this struct is the impl.
//
// We pre-walk the tree in post-order once so the per-column scorer can
// iterate by integer index instead of pointer-chasing.  Each node carries
// its branch length, parent index (-1 for root), and either a leaf id
// (>= 0) or -1 (internal).  leaf_idx maps the canonical leaf index back to
// its post-order position, used to populate the scratch from a base array
// indexed by leaf id.
/////////////////////////////////////////////////////////////////////////////

struct GerpTree {
    stTree *root;             // we own; freed in destruct
    int64_t n_nodes;
    int64_t n_leaves;

    int64_t *parent;          // [n_nodes]; -1 for root
    double  *branch_len;      // [n_nodes]; branch above this node (root: 0)
    int64_t *leaf_id_of_node; // [n_nodes]; -1 for internal nodes
    int64_t *node_of_leaf;    // [n_leaves]; post-order index of each leaf

    char **leaf_name;         // [n_leaves]; borrowed from stTree labels
    stSet  *internal_labels;  // ancestor set (labels of all internal nodes)
    stSet  *all_labels;       // every labelled node, for extract_genome_name
    stHash *leaf_by_label;    // leaf label -> (int64_t *) leaf id
};

// Post-order walk that fills parent[], branch_len[], leaf_id_of_node[],
// node_of_leaf[], and leaf_name[] for one node.  Returns this node's
// post-order index.  Children are visited before parents (so a node's
// index is always > its descendants').  parent[] is patched after each
// child returns, since the child's index isn't known until it claims its
// slot in the recursion below.
static int64_t gt_walk(stTree *node, int64_t parent_idx, GerpTree *gt,
                       int64_t *cur) {
    int64_t nc = stTree_getChildNumber(node);
    int64_t  stack_buf[16];
    int64_t *children = (nc <= 16) ? stack_buf : st_malloc((size_t)nc * sizeof(int64_t));
    for (int64_t i = 0; i < nc; i++) {
        // Pass parent_idx=-1 as placeholder; we backfill once we know
        // our own index below.
        children[i] = gt_walk(stTree_getChild(node, i), -1, gt, cur);
    }
    int64_t my_idx = (*cur)++;
    for (int64_t i = 0; i < nc; i++) gt->parent[children[i]] = my_idx;
    if (nc > 16) free(children);

    gt->parent[my_idx]      = parent_idx;   // root caller passes -1
    gt->branch_len[my_idx]  = stTree_getBranchLength(node);
    if (nc == 0) {
        int64_t leaf_id = gt->n_leaves++;
        gt->leaf_id_of_node[my_idx] = leaf_id;
        gt->node_of_leaf[leaf_id]   = my_idx;
        const char *lbl = stTree_getLabel(node);
        gt->leaf_name[leaf_id] = lbl ? (char *) lbl : NULL;
    } else {
        gt->leaf_id_of_node[my_idx] = -1;
    }
    return my_idx;
}

// First pre-pass: count nodes + leaves so we can size all arrays.
static void gt_count(stTree *node, int64_t *n_nodes, int64_t *n_leaves) {
    int64_t nc = stTree_getChildNumber(node);
    (*n_nodes)++;
    if (nc == 0) (*n_leaves)++;
    for (int64_t i = 0; i < nc; i++) gt_count(stTree_getChild(node, i), n_nodes, n_leaves);
}

// Collect ALL labels (leaves + internals) into all_labels, and internal
// labels into internal_labels.  Borrows the label strings from stTree.
static void gt_collect_labels(stTree *node, stSet *all_labels, stSet *internal_labels) {
    int64_t nc = stTree_getChildNumber(node);
    const char *lbl = stTree_getLabel(node);
    if (lbl != NULL && *lbl != '\0') {
        // stSet_insert with stHash_stringKey expects char* keys; copy so
        // we own the strings (tree owns its labels and we'll be queried
        // long after).
        stSet_insert(all_labels, stString_copy(lbl));
        if (nc > 0) stSet_insert(internal_labels, stString_copy(lbl));
    }
    for (int64_t i = 0; i < nc; i++) {
        gt_collect_labels(stTree_getChild(node, i), all_labels, internal_labels);
    }
}

GerpTree *gerp_tree_construct(const char *newick) {
    if (newick == NULL || *newick == '\0') return NULL;
    // stTree_parseNewickString takes a non-const pointer in the sonLib
    // signature but doesn't mutate it -- safe to cast.
    stTree *root = stTree_parseNewickString((char *) newick);
    if (root == NULL) return NULL;

    GerpTree *gt = st_calloc(1, sizeof(GerpTree));
    gt->root = root;
    gt_count(root, &gt->n_nodes, &gt->n_leaves);

    gt->parent          = st_malloc((size_t)gt->n_nodes * sizeof(int64_t));
    gt->branch_len      = st_malloc((size_t)gt->n_nodes * sizeof(double));
    gt->leaf_id_of_node = st_malloc((size_t)gt->n_nodes * sizeof(int64_t));
    gt->node_of_leaf    = st_malloc((size_t)(gt->n_leaves ? gt->n_leaves : 1) * sizeof(int64_t));
    gt->leaf_name       = st_calloc((size_t)(gt->n_leaves ? gt->n_leaves : 1), sizeof(char *));

    int64_t cur = 0;
    int64_t expected_n_leaves = gt->n_leaves;
    gt->n_leaves = 0;  // gt_walk increments as it goes
    gt_walk(root, -1, gt, &cur);
    assert(cur == gt->n_nodes && gt->n_leaves == expected_n_leaves);

    // Label sets.
    gt->all_labels      = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    gt->internal_labels = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    gt_collect_labels(root, gt->all_labels, gt->internal_labels);

    // leaf label -> leaf id, for fast row resolution.  Values are owned
    // int64_t allocations.
    gt->leaf_by_label = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free);
    for (int64_t i = 0; i < gt->n_leaves; i++) {
        if (gt->leaf_name[i] != NULL) {
            int64_t *v = st_malloc(sizeof(int64_t));
            *v = i;
            stHash_insert(gt->leaf_by_label, gt->leaf_name[i], v);
        }
    }
    return gt;
}

void gerp_tree_destruct(GerpTree *gt) {
    if (gt == NULL) return;
    if (gt->root != NULL) stTree_destruct(gt->root);
    free(gt->parent);
    free(gt->branch_len);
    free(gt->leaf_id_of_node);
    free(gt->node_of_leaf);
    free(gt->leaf_name);  // entries borrowed; just free the array
    if (gt->internal_labels) stSet_destruct(gt->internal_labels);
    if (gt->all_labels)      stSet_destruct(gt->all_labels);
    if (gt->leaf_by_label)   stHash_destruct(gt->leaf_by_label);
    free(gt);
}

bool gerp_tree_is_ancestor(const GerpTree *gt, const char *name) {
    return stSet_search(gt->internal_labels, (void *) name) != NULL;
}

int64_t gerp_tree_n_leaves(const GerpTree *gt) {
    return gt->n_leaves;
}

/////////////////////////////////////////////////////////////////////////////
// Per-column scoring
//
// One Fitch pass + branch sum, fused into a single post-order walk.  Char
// sets are 4-bit bitmasks (A=1, C=2, G=4, T=8); intersect = AND, union =
// OR, empty test = (==0).  "Present" means the node's subtree contains at
// least one surviving leaf in this column.
/////////////////////////////////////////////////////////////////////////////

// Per-node Fitch state for one column.  `count[4]` is the votes-per-base
// tally built up by folding children into their parent.  At each internal
// node's own post-order visit we collapse count[] into cset using the
// textbook Hartigan/Fitch batch rule:
//     max_k    = max over bases of count[b]
//     cset     = bitmask of bases with count[b] == max_k
//     cost   += children_present - max_k
// The accumulator-style "intersect-or-union as each child arrives" form
// gives the wrong cost on polytomies (>= 3 children with non-pairwise-
// disjoint sets) -- audit caught this on real cactus trees that have
// non-binary internal nodes.  count[] is uint16_t since cactus polytomies
// at vertebrate scale stay well under 2^16 children.
struct GerpScratch {
    uint8_t  *cset;             // [n_nodes]; 4-bit base set
    uint8_t  *present;          // [n_nodes]; 1 if subtree has a leaf in this column
    uint16_t *count;            // [n_nodes * 4]; per-base vote tally at internals
    uint16_t *children_present; // [n_nodes]; # of present children folded into me
    int64_t   n_nodes;
};

GerpScratch *gerp_scratch_construct(const GerpTree *gt) {
    GerpScratch *sc = st_calloc(1, sizeof(GerpScratch));
    sc->n_nodes          = gt->n_nodes;
    sc->cset             = st_malloc((size_t)sc->n_nodes);
    sc->present          = st_malloc((size_t)sc->n_nodes);
    sc->count            = st_malloc((size_t)sc->n_nodes * 4 * sizeof(uint16_t));
    sc->children_present = st_malloc((size_t)sc->n_nodes * sizeof(uint16_t));
    return sc;
}

void gerp_scratch_destruct(GerpScratch *sc) {
    if (sc == NULL) return;
    free(sc->cset);
    free(sc->present);
    free(sc->count);
    free(sc->children_present);
    free(sc);
}

static inline uint8_t base_to_bit(char b) {
    switch (b) {
        case 'A': return 0x1;
        case 'C': return 0x2;
        case 'G': return 0x4;
        case 'T': return 0x8;
        default:  return 0;
    }
}

bool gerp_score_column(const GerpTree *gt, GerpScratch *sc,
                       const char *leaf_bases, int64_t min_leaves,
                       double branch_scale,
                       double *out_rs, int64_t *out_depth) {
    // Initialise per-node state.  Leaves get their cset/present set from
    // the column; internals are zeroed (count[] votes accumulated during
    // the post-order fold below).
    int64_t depth = 0;
    for (int64_t v = 0; v < gt->n_nodes; v++) {
        if (gt->leaf_id_of_node[v] >= 0) {
            int64_t i = gt->leaf_id_of_node[v];
            char b = leaf_bases[i];
            uint8_t bit = (b == 0) ? 0 : base_to_bit((char)toupper((unsigned char)b));
            sc->cset[v]    = bit;
            sc->present[v] = bit ? 1 : 0;
            if (bit) depth++;
            // Leaves never receive child folds, but zero anyway so an
            // accidentally-routed fold would be detected as garbage.
            sc->count[4*v + 0] = sc->count[4*v + 1] = 0;
            sc->count[4*v + 2] = sc->count[4*v + 3] = 0;
            sc->children_present[v] = 0;
        } else {
            sc->cset[v]    = 0;
            sc->present[v] = 0;
            sc->count[4*v + 0] = sc->count[4*v + 1] = 0;
            sc->count[4*v + 2] = sc->count[4*v + 3] = 0;
            sc->children_present[v] = 0;
        }
    }
    *out_depth = depth;
    if (depth < min_leaves) { *out_rs = 0.0; return false; }

    // Post-order walk.  Two responsibilities per node v:
    //   1. If v is internal: derive v.cset / v.present from v.count[]
    //      via the batch Hartigan/Fitch rule (max-count bases, cost +=
    //      children_present - max_k).
    //   2. If v.present and v has a parent p: fold v into p (count each
    //      base in v.cset onto p.count[b], bump p.children_present, and
    //      add v's branch-to-parent into branch_sum).
    //
    // Order matters: leaves' fold step uses their already-set cset; an
    // internal node's step-1 must complete before step-2 (since step-2
    // reads v.cset).  The single post-order pass below handles both
    // because children precede their parent in the iteration order.
    double cost = 0;
    double branch_sum = 0;

    for (int64_t v = 0; v < gt->n_nodes; v++) {
        // Step 1 (internals only): collapse count[] -> cset + cost.
        if (gt->leaf_id_of_node[v] < 0) {
            uint16_t cp = sc->children_present[v];
            if (cp == 0) {
                // empty subtree: stays present=0, cset=0
                continue;
            }
            uint16_t c0 = sc->count[4*v + 0];
            uint16_t c1 = sc->count[4*v + 1];
            uint16_t c2 = sc->count[4*v + 2];
            uint16_t c3 = sc->count[4*v + 3];
            uint16_t mx = c0;
            if (c1 > mx) mx = c1;
            if (c2 > mx) mx = c2;
            if (c3 > mx) mx = c3;
            uint8_t cs = 0;
            if (c0 == mx) cs |= 0x1;
            if (c1 == mx) cs |= 0x2;
            if (c2 == mx) cs |= 0x4;
            if (c3 == mx) cs |= 0x8;
            sc->cset[v]    = cs;
            sc->present[v] = 1;
            cost += (double)(cp - mx);
        }
        // Step 2: fold v into parent (regardless of leaf/internal).  Only
        // present nodes contribute -- absent subtrees aren't in the
        // pruned tree at all.
        if (!sc->present[v]) continue;
        int64_t p = gt->parent[v];
        if (p < 0) continue;                    // root: no parent, no branch
        branch_sum += gt->branch_len[v];
        uint8_t cs = sc->cset[v];
        if (cs & 0x1) sc->count[4*p + 0]++;
        if (cs & 0x2) sc->count[4*p + 1]++;
        if (cs & 0x4) sc->count[4*p + 2]++;
        if (cs & 0x8) sc->count[4*p + 3]++;
        sc->children_present[p]++;
    }

    *out_rs = branch_sum * branch_scale - cost;
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Per-block row resolution
/////////////////////////////////////////////////////////////////////////////

// Resolve a "genome.seq" style name to its genome label by trying each
// '.' as a split point and accepting the first prefix that's in
// gt->all_labels.  Mirrors extract_genome_name's algorithm but returns
// NULL on failure instead of calling st_errAbort (the version in taf.c
// aborts the process when called with hal_species != NULL).  Returns a
// borrowed pointer into gt->all_labels (do NOT free).
static const char *gerp_resolve_genome(const GerpTree *gt, const char *sequence_name) {
    size_t n = strlen(sequence_name);
    if (n == 0) return NULL;
    // Walk through every '.' position in the name and ask the label set.
    // We can't strdup-and-test because all_labels has its own copies, so
    // we lookup with a temporary buffer (stack for short names).
    char  stack_buf[512];
    char *buf = (n + 1 <= sizeof(stack_buf)) ? stack_buf : st_malloc(n + 1);
    size_t off = 0;
    const char *match = NULL;
    while (off < n) {
        const char *dot = strchr(sequence_name + off, '.');
        if (dot == NULL || dot == sequence_name + n - 1) break;  // no usable dot
        if (dot == sequence_name) { off++; continue; }           // leading dot, skip
        size_t pre_n = (size_t)(dot - sequence_name);
        memcpy(buf, sequence_name, pre_n);
        buf[pre_n] = '\0';
        if (stSet_search(gt->all_labels, buf) != NULL) {
            // stSet stores a copy of the label; look up the canonical
            // pointer via the leaf or internal-label set.  We just want
            // *some* stable pointer into the tree's label space.
            match = stSet_search(gt->all_labels, buf);
            break;
        }
        off = pre_n + 1;
    }
    if (buf != stack_buf) free(buf);
    return match;
}

int gerp_block_resolve_rows(const GerpTree *gt, const Alignment *aln,
                            bool skip_paralogs,
                            Alignment_Row **row_by_leaf,
                            const char **unknown_seq) {
    for (int64_t i = 0; i < gt->n_leaves; i++) row_by_leaf[i] = NULL;
    if (unknown_seq) *unknown_seq = NULL;

    Alignment_Row *row = aln->row;
    bool paralog_seen = false;
    while (row != NULL) {
        const char *genome = gerp_resolve_genome(gt, row->sequence_name);
        if (genome == NULL) {
            if (unknown_seq) *unknown_seq = row->sequence_name;
            return GERP_BLOCK_UNKNOWN_SPECIES;
        }
        // Ancestor row: drop silently.
        if (stSet_search(gt->internal_labels, (void *)genome) != NULL) {
            row = row->n_row;
            continue;
        }
        // Leaf row: look up its tree-leaf index.
        int64_t *p = stHash_search(gt->leaf_by_label, (void *)genome);
        if (p == NULL) {
            // genome resolved to a label in all_labels but isn't a leaf or
            // an internal -- shouldn't be possible since we built both
            // sets from the same tree, but treat as hard error to surface
            // any future drift.
            if (unknown_seq) *unknown_seq = row->sequence_name;
            return GERP_BLOCK_UNKNOWN_SPECIES;
        }
        int64_t leaf_id = *p;
        if (row_by_leaf[leaf_id] != NULL) {
            paralog_seen = true;
            // Keep first-seen row; v1 paralog policy is "skip the block",
            // not "use one of the paralogs".
        } else {
            row_by_leaf[leaf_id] = row;
        }
        row = row->n_row;
    }
    if (paralog_seen && skip_paralogs) return GERP_BLOCK_SKIP;
    return GERP_BLOCK_OK;
}
