/*
 * view_chain_dup_filter: per-region chain-based duplicate-row filter for
 * taffy view.  For each non-row-0 target genome that appears multiple
 * times across the buffered blocks, partition its rows into chains via
 * taffy_chain() (using block-index as the query axis and the target row
 * coords as the target axis) and apply the overlap-frac filter (chain.h)
 * to keep only chains whose q-coverage doesn't overlap higher-scoring
 * survivors by more than overlap_frac.
 *
 * Per-block-only chaining is degenerate (every row in a single block has
 * the same column range, so chain-length collapses to row length).  The
 * chaining signal lives ACROSS blocks: same-target-genome rows in
 * adjacent blocks should join into one chain if they are collinear w.r.t.
 * row-0.  Hence the buffered-region API.
 *
 * Row-0 is always preserved.  Target genomes appearing in only one block
 * are always preserved (no chain pass).  This mutates each Alignment's
 * row linked list in place (drops non-survivor rows via
 * alignment_row_destruct + row_number--).
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#ifndef TAF_VIEW_CHAIN_DUP_FILTER_H_
#define TAF_VIEW_CHAIN_DUP_FILTER_H_

#include "taf.h"
#include "sonLib.h"

/* In-place chain dup filter.
 *
 *   blocks:        stList of Alignment* (caller owns).  Each Alignment's
 *                  row list may be mutated; the Alignment pointer itself
 *                  is unchanged.
 *   overlap_frac:  threshold in [0, 1].  0 = drop on ANY q-overlap with
 *                  higher-scoring kept chains (recommended for paralogy);
 *                  1 = essentially keep-all.
 *   cap:           optional hard cap on survivor count per target genome
 *                  (0 = no cap).
 */
void view_chain_dup_filter(stList *blocks, double overlap_frac, int64_t cap);

#endif /* TAF_VIEW_CHAIN_DUP_FILTER_H_ */
