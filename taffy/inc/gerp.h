#ifndef TAFFY_GERP_H_
#define TAFFY_GERP_H_

#include <stdint.h>
#include <stdbool.h>
#include "sonLib.h"
#include "taf.h"

/*
 * GERP RS (Rejected Substitutions) conservation scoring.
 *
 * For each MAF/TAF block:
 *   - row 0 is the "reference" (whose sequence_name + start position keys the
 *     wig records).  In a normal hg38-anchored MAF that's hg38; in a
 *     universal MAF it's the anchor ancestor of that block (and so the wig
 *     ends up with many chroms, one per ancestor used as anchor).
 *   - rows whose extracted genome name is an internal tree node label are
 *     dropped (they're ML reconstructions, not observations).
 *   - the remaining "leaf" rows feed the per-column scorer.
 *
 * For each column:
 *   - collect (leaf, base) where base is A/C/G/T (case-folded); gap or N
 *     drops that leaf for this column.
 *   - if surviving leaf count < min_leaves, skip the column (no wig record).
 *   - otherwise:
 *       branch_sum = sum of branch lengths in the tree pruned to surviving
 *                    leaves, scaled by branch_scale.
 *       cost      = Fitch parsimony cost of the column on the pruned tree.
 *       RS         = branch_sum - cost.
 *
 * The wig record's position is row-0's start + (count of non-gap row-0 bases
 * to the left of this column).  Columns where row-0 itself is a gap get no
 * wig record (they have no anchor coordinate).
 *
 * Unknown-species handling: if a (non-ancestor) leaf row's extracted genome
 * isn't found in the tree, gerp_score_block returns a non-zero error code
 * (hard error -- per the v1 decision).
 */

/////////////////////////////////////////////////////////////////////////////
// Tree state
/////////////////////////////////////////////////////////////////////////////

typedef struct GerpTree GerpTree;

/*
 * Parse a Newick string and build the gerp tree state.  All internal-node
 * labels become the "ancestor set"; all leaf labels become the "species
 * set" (looked up via extract_genome_name on a row's sequence_name).
 *
 * Returns NULL on parse failure.
 */
GerpTree *gerp_tree_construct(const char *newick);

void gerp_tree_destruct(GerpTree *gt);

/* Was `name` an internal-node label in the tree?  Used to drop ancestral
 * rows from a MAF block. */
bool gerp_tree_is_ancestor(const GerpTree *gt, const char *name);

/* Number of leaves in the parsed tree. */
int64_t gerp_tree_n_leaves(const GerpTree *gt);

/*
 * Walk up the tree from the node labelled `name` until we hit a node
 * whose label is in `targets` (a stSet of char* labels, e.g. the
 * lineage-roots set for gerp-stats clade attribution).  Returns the
 * matching label (borrowed pointer into the tree's label storage) or
 * NULL if no ancestor of `name` (including itself) is in `targets`.
 *
 * Used by taffy gerp-stats to assign each column's anchor ancestor to
 * a user-defined lineage bucket.  O(tree depth) per call; the caller
 * is expected to cache results per unique chrom name.
 */
const char *gerp_tree_walk_to_set(const GerpTree *gt,
                                   const char *name,
                                   stSet *targets);

/*
 * Depth-from-root (= edge count from root to the node labelled `name`).
 * Returns -1 if `name` isn't in the tree.  Used by gerp-stats to sort
 * lineage roots by depth for the TSV output.
 */
int64_t gerp_tree_depth_from_root(const GerpTree *gt, const char *name);

/*
 * Resolve a "<genome>.<seq>" style name to its genome label by trying
 * each '.' as a split point and accepting the first prefix that matches
 * a tree label.  Same algorithm as the internal MAF-row resolver, but
 * exposed for gerp-stats's chrom -> clade attribution.  Returns a
 * borrowed pointer into the tree's label storage (or NULL).
 */
const char *gerp_tree_resolve_genome(const GerpTree *gt, const char *seq_name);

/////////////////////////////////////////////////////////////////////////////
// Per-block scorer state
/////////////////////////////////////////////////////////////////////////////

/*
 * Per-thread scratch -- reused across blocks.  Holds the per-node character
 * sets + presence flags for one Fitch pass over the tree.  Allocate once
 * (sized to gt->n_nodes), feed to gerp_score_column on each column.
 */
typedef struct GerpScratch GerpScratch;

GerpScratch *gerp_scratch_construct(const GerpTree *gt);
void gerp_scratch_destruct(GerpScratch *sc);

/*
 * Map A/C/G/T (any case) -> 4-bit set bit (A=1, C=2, G=4, T=8).  Anything
 * else (gap, N, NUL) -> 0.  Inlined in the header so the per-column score
 * loop in taf_gerp.c can build leaf csets without a cross-TU call.
 */
static inline uint8_t gerp_base_to_bit(char b) {
    switch (b) {
        case 'A': case 'a': return 0x1;
        case 'C': case 'c': return 0x2;
        case 'G': case 'g': return 0x4;
        case 'T': case 't': return 0x8;
        default:  return 0;
    }
}

/*
 * Score one column.  `leaf_bases` is an array of length gt->n_leaves
 * indexed by the tree's leaf order; entries are A/C/G/T (case-folded) for
 * present leaves and 0 (NUL) for absent (gap, N, or species not in this
 * block).  See gerp_block_state below for how to populate it efficiently.
 *
 * Returns the RS value via *out_rs and the count of A/C/G/T leaves via
 * *out_depth.  If depth < min_leaves the column is unscored: *out_rs is 0
 * and the function returns false.
 *
 * Thin wrapper over gerp_score_column_csets; suitable when each leaf
 * contributes at most one base (no paralog union).
 */
bool gerp_score_column(const GerpTree *gt, GerpScratch *sc,
                       const char *leaf_bases, int64_t min_leaves,
                       double branch_scale,
                       double *out_rs, int64_t *out_depth);

/*
 * Score one column from per-leaf character SETS (4-bit bitmasks, one per
 * leaf, indexed by tree leaf order).  Used by taffy gerp's UNION paralog
 * policy: multiple rows of the same species in one block contribute their
 * bases OR'd into the leaf's cset, so a paralog leaf with rows {A, T}
 * becomes `0x9` (Hartigan/Fitch treats {A,T} as a polymorphic leaf -- the
 * column doesn't count a substitution against any species that agrees
 * with EITHER paralog).
 *
 * depth is the count of leaves with non-empty cset (= count of species
 * contributing any base, regardless of how many paralog rows folded in).
 * Returns false (out_rs = 0) if depth < min_leaves.
 */
bool gerp_score_column_csets(const GerpTree *gt, GerpScratch *sc,
                             const uint8_t *leaf_csets, int64_t min_leaves,
                             double branch_scale,
                             double *out_rs, int64_t *out_depth);

/////////////////////////////////////////////////////////////////////////////
// Per-block driver
/////////////////////////////////////////////////////////////////////////////

/*
 * Paralog policy controls how a block with > 1 row from the same species
 * is handled.
 *
 *   UNION: keep every leaf row.  The per-column scorer ORs paralog bases
 *          into a multi-bit cset for the leaf -- biologically "any of these
 *          paralogs' bases count as that species' character set".  This is
 *          the default: it preserves the column instead of dropping it,
 *          and stays consistent with Hartigan/Fitch's multi-state-leaf
 *          semantics.
 *   SKIP : drop the entire block on any paralog.  Matches strict GERP++
 *          semantics and the v1 default; keep for reproducibility.
 *   FIRST: keep the first-seen row per leaf, drop subsequent paralogs.
 *          A pragmatic middle ground.
 */
typedef enum {
    GERP_PARALOG_UNION = 0,
    GERP_PARALOG_SKIP  = 1,
    GERP_PARALOG_FIRST = 2,
} GerpParalogPolicy;

/*
 * One (leaf, row) pair from the block.  In UNION mode there may be > 1
 * entry with the same leaf_id (the paralogs).  Caller iterates the entries
 * array per column and OR's each row's base bit into leaf_csets[leaf_id].
 *
 * `row` is borrowed from the Alignment (do NOT free).
 */
typedef struct {
    int64_t        leaf_id;
    Alignment_Row *row;
} GerpRowEntry;

/*
 * Map a block's rows to tree leaves.  Built once per block in the main
 * scorer:
 *   - ancestor rows (genome in the tree's internal-label set) are dropped
 *   - the surviving rows are appended to entries[] as (leaf_id, row) pairs
 *   - paralog handling depends on `policy` (see GerpParalogPolicy above)
 *   - if a surviving row's genome isn't in the tree at all, this is a
 *     hard error (the v1 decision -- caller checks the return code).
 *
 * Caller owns `entries` and sizes `entries_cap` to at least aln->row_number
 * (enough to hold every row in the block even in UNION mode).  On return:
 *   *n_active        = number of entries filled
 *   *n_paralog_dups  = number of paralog rows beyond first-per-leaf
 *                      (always 0 in SKIP mode -- SKIP returns early on
 *                       first paralog)
 *
 * On unknown-species error (rc == GERP_BLOCK_UNKNOWN_SPECIES), the offending
 * row's full sequence_name is returned via *unknown_seq for the caller to
 * include in its abort message.  Caller does NOT free it (it's borrowed
 * from the Alignment).
 *
 * Returns:
 *   GERP_BLOCK_OK        (0) -- block ok, score it
 *   GERP_BLOCK_SKIP      (1) -- skip block (paralog with policy=SKIP)
 *   GERP_BLOCK_UNKNOWN_SPECIES (2) -- hard error; *unknown_seq is set
 */
#define GERP_BLOCK_OK                0
#define GERP_BLOCK_SKIP              1
#define GERP_BLOCK_UNKNOWN_SPECIES   2
int gerp_block_resolve_rows(const GerpTree *gt, const Alignment *aln,
                            GerpParalogPolicy policy,
                            GerpRowEntry *entries, int64_t entries_cap,
                            int64_t *n_active, int64_t *n_paralog_dups,
                            const char **unknown_seq);

#endif /* TAFFY_GERP_H_ */
