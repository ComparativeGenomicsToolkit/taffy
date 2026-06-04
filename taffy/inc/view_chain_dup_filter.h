/*
 * view_chain_dup_filter: per-region chain-based duplicate-row filter for
 * taffy view.  For each non-row-0 target genome that appears multiple times
 * across the buffered blocks, partition its rows into chains via
 * taffy_chain() (using row-0 as the query axis and the target row as the
 * target axis) and keep only rows belonging to the top N chains.
 *
 * Per-block-only chaining is degenerate (every row in a single block has
 * the same column range, so chain-length collapses to row length).  The
 * chaining signal lives ACROSS blocks: same-target-genome rows in adjacent
 * blocks should join into one chain if they are collinear w.r.t. row-0.
 * Hence the buffered-region API.
 *
 * Row-0 is always preserved.  Target genomes appearing in only one block
 * are always preserved (no chain pass; trivially top-N).  This mutates
 * each Alignment's row linked list in place (drops non-survivor rows via
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
 *   blocks: stList of Alignment* (caller owns).  Each Alignment's row list
 *           may be mutated; the Alignment pointer itself is unchanged.
 *   top_n:  number of top-scoring chains per target genome to keep
 *           (clamped to >= 1).
 */
void view_chain_dup_filter(stList *blocks, int64_t top_n);

#endif /* TAF_VIEW_CHAIN_DUP_FILTER_H_ */
