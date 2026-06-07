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
 *   blocks:        stList of Alignment* (caller owns).  ALL blocks are
 *                  used as chain-pass context, but row mutation is only
 *                  applied to blocks at indices in [apply_lo, apply_hi).
 *                  This split lets a windowed caller defer mutation of
 *                  the K-block carryover until later windows give it
 *                  more forward context.  Pass apply_lo=0,
 *                  apply_hi=stList_length(blocks) for the simple
 *                  "filter the whole buffer" case.
 *   overlap_frac:  threshold in [0, 1].  0 = drop on ANY q-overlap with
 *                  higher-scoring kept chains (recommended for paralogy);
 *                  1 = essentially keep-all.
 *   cap:           optional hard cap on survivor count per target genome
 *                  (0 = no cap).
 *   apply_lo,      half-open mutation range [apply_lo, apply_hi).  Pass
 *   apply_hi       (0, stList_length(blocks)) to mutate all blocks.
 */
void view_chain_dup_filter(stList *blocks, double overlap_frac, int64_t cap,
                           int64_t apply_lo, int64_t apply_hi);

#endif /* TAF_VIEW_CHAIN_DUP_FILTER_H_ */
