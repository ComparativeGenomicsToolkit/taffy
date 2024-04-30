#ifndef TAF_TAI_H_
#define TAF_TAI_H_

/*
 * Index a TAF file.  Modeled on .fai from samtools faidx
 * 
 * Each line represents a position on a reference contig as mapped to a position in the TAF/MAF file
 * 
 * The lines must be increasing in [Contig-Name, Start-position] order. 
 *
 * Lines can reflect absolute coordinates or relative coordinates. 
 *
 * Absoluate Coordinates:
 *
 * Column 1: Contig Name
 * Column 2: Position in Contig (0-based)
 * Column 3: Offset in bytes in MAF/TAF file being indexed (can be bgzipped)
 *
 * Relative coordinates (to the previous line):
 *
 * Column 1: *
 * Column 2: Positiion in Contig *relative to previous line*
 * Column 3: Offset in file *relative to previous line*
 *
 * For example, ths index:

Anc0.Anc0refChr11	0	1421634189773
*	10007	25769835729
*	10001	30064728450
*	10077	30064782171

* is equivalent to:

Anc0.Anc0refChr11	0	1421634189773
Anc0.Anc0refChr11	10007	1447404025502
Anc0.Anc0refChr11	20008	30064728450
Anc0.Anc0refChr11	30085	30064782171
 *
 */

#include "line_iterator.h"

typedef struct _Tai {
    stSortedSet *idx;
    stList *names; // just to keep track of memory -- we only keep one instance of each sequence name
    bool maf;
} Tai;

typedef struct _TaiIt {
    char *name;
    // these are bed-like 0-based half-open
    int64_t start;
    int64_t end;
    Alignment *alignment;
    Alignment *p_alignment;
    bool run_length_encode_bases;
    bool maf;
} TaiIt;


/* Return taf_path + .tai 
 */
char *tai_path(const char *taf_path);

/*
 * Parse a region into contig / start / length, where subrange is optional
 * chr1:10-13 -> chr1 / 10 / 3
 * chr1:10 -> chr1 / 10 / 1
 * chr1 -> chr1 / 0 / LONG_MAX
 */
char *tai_parse_region(const char* region, int64_t *start, int64_t *length);

/*
 * Make an index of a TAF in "idx_fh" 
 * The index is made on the "reference" first contig of each block
 * For each such contig, the index will have one line for each index_block_size
 * region of it found in the TAF. 
 */
int tai_create(LI *li, FILE* idx_fh, int64_t index_block_size);

/*
 * Load the index from disk
 */
Tai *tai_load(FILE* idx_fh, bool maf);

/*
 * Free the index
 */
void tai_destruct(Tai* idx);

/*
 * Query the taf index
 * start is 0-based
 * length can be -1 for everything
 * If the given start point is not in the index, returns an "empty" iterator (tai_next returns NULL).
 */
TaiIt *tai_iterator(Tai* idx, LI *li, bool run_length_encode_bases, const char *contig, int64_t start, int64_t length);

/*
 * Iterate through a region as obtained via the iterator
 */
Alignment *tai_next(TaiIt *tai_it, LI *li);

/*
 * Check if tai_next() is not NULL
 */
bool tai_has_next(TaiIt *tai_it);

/*
 * Free a tai iterator
 */
void tai_iterator_destruct(TaiIt *tai_it);

/**
 * Return a map of Sequence name to Length. Only reference (ie indexed) sequences are returned
 */
stHash *tai_sequence_lengths(Tai *idx, LI *li);

#endif
