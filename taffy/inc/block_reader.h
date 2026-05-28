#ifndef TAF_BLOCK_READER_H_
#define TAF_BLOCK_READER_H_

/*
 * Format-agnostic block reader. Sniffs MAF vs TAF on construction, reads the
 * appropriate header, and dispatches subsequent block reads to the right
 * underlying reader. For MAF input, also handles cross-block linkage via
 * alignment_link_adjacent so that downstream tools see the same coordinate
 * continuity they get from a TAF input.
 *
 * Intended to let CLI tools that previously hardcoded TAF (norm, sort,
 * coverage, annotate, add-gap-bases, stats, ...) accept either format with
 * minimal per-tool churn.
 */

#include "taf.h"

typedef struct _BlockReader BlockReader;

/*
 * Open a reader on an existing LI handle. Sniffs the input format from the
 * first header line, reads the header, and stashes any flags (e.g. RLE) that
 * tools may need to consult. Caller retains ownership of li.
 *
 * Returns NULL on unrecognized format (with an error printed to stderr).
 */
BlockReader *block_reader_open(LI *li);

/*
 * Take ownership of the parsed header tag. Subsequent calls return NULL.
 * The caller is responsible for tag_destruct'ing it (or passing it to a
 * header writer which will). After this call the BlockReader still owns
 * its other state and must still be destructed.
 */
Tag *block_reader_take_header(BlockReader *r);

bool block_reader_is_maf(const BlockReader *r);
bool block_reader_run_length_encoded(const BlockReader *r);

/*
 * Read the next block. For MAF inputs, prev_block is used to link adjacent
 * rows (matching the behavior of taffy view's MAF read path). Pass NULL for
 * the first block. Returns NULL at EOF.
 *
 * Lifetime note: the returned alignment is owned by the caller and must be
 * destructed (alignment_destruct) when no longer needed. The reader does not
 * retain a reference to it.
 */
Alignment *block_reader_next(BlockReader *r, Alignment *prev_block);

void block_reader_destruct(BlockReader *r);

#endif
