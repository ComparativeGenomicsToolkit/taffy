#ifndef TAF_INDEX_H_
#define TAF_INDEX_H_

/*
 * Index a TAF file.  Modeled on .fai from samtools faidx
 * Index format is 
 * Sequence Length Offset Mean-Row-Length TAF-header
 * where TAF header is everything we need to get going parsing the TAF file at that exact position
 *
 */


#include "line_iterator.h"

/*
 * Make an index of a TAF in "idx_fp" 
 * The index is made on the "reference" first contig of each block
 * For each such contig, the index will have one line for each index_block_size
 * region of it found in the TAF. 
 */
int index_taf(LI *li, FILE* idx_fp, int64_t index_block_size);

#endif
