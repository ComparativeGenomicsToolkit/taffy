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
 * Score one column.  `leaf_bases` is an array of length gt->n_leaves
 * indexed by the tree's leaf order; entries are A/C/G/T (case-folded) for
 * present leaves and 0 (NUL) for absent (gap, N, or species not in this
 * block).  See gerp_block_state below for how to populate it efficiently.
 *
 * Returns the RS value via *out_rs and the count of A/C/G/T leaves via
 * *out_depth.  If depth < min_leaves the column is unscored: *out_rs is 0
 * and the function returns false.
 */
bool gerp_score_column(const GerpTree *gt, GerpScratch *sc,
                       const char *leaf_bases, int64_t min_leaves,
                       double branch_scale,
                       double *out_rs, int64_t *out_depth);

/////////////////////////////////////////////////////////////////////////////
// Per-block driver
/////////////////////////////////////////////////////////////////////////////

/*
 * Map a block's rows to tree leaves.  Built once per block in the main
 * scorer:
 *   - ancestor rows (genome in the tree's internal-label set) are dropped
 *   - the surviving rows are bucketed by tree-leaf index
 *   - if any bucket has > 1 row AND skip_paralogs is true, the whole block
 *     is skipped
 *   - if a surviving row's genome isn't in the tree at all, this is a
 *     hard error (the v1 decision -- caller checks the return code).
 *
 * `row_by_leaf` is owned by the caller and sized to gerp_tree_n_leaves(gt);
 * on return, row_by_leaf[i] points to the Alignment_Row for leaf i (or NULL
 * if leaf i isn't in this block).
 *
 * On unknown-species error (rc == GERP_BLOCK_UNKNOWN_SPECIES), the offending
 * row's full sequence_name is returned via *unknown_seq for the caller to
 * include in its abort message.  Caller does NOT free it (it's borrowed
 * from the Alignment).
 *
 * Returns:
 *   GERP_BLOCK_OK        (0) -- block ok, score it
 *   GERP_BLOCK_SKIP      (1) -- skip block (paralog with skip_paralogs)
 *   GERP_BLOCK_UNKNOWN_SPECIES (2) -- hard error; *unknown_seq is set
 */
#define GERP_BLOCK_OK                0
#define GERP_BLOCK_SKIP              1
#define GERP_BLOCK_UNKNOWN_SPECIES   2
int gerp_block_resolve_rows(const GerpTree *gt, const Alignment *aln,
                            bool skip_paralogs,
                            Alignment_Row **row_by_leaf,
                            const char **unknown_seq);

#endif /* TAFFY_GERP_H_ */
