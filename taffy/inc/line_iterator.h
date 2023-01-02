#ifndef STLINE_ITERATOR_H_
#define STLINE_ITERATOR_H_

#include <stdio.h>
#include <stdint.h>
#include "sonLib.h"

typedef struct BGZF BGZF;

/*
 * Simple wrapper around a FILE handle that allows you to iterator over lines from a text file
 * and "peek" at lines before choosing to get them.
 */

struct BGZF;

typedef struct _LI {
#ifdef USE_HTSLIB
    BGZF *bgzf;
#else
    FILE *fh;
#endif
    char *line;
    int64_t prev_pos; // position before reading the current buffer
    int64_t pos;      // position after reading the curent buffer    
} LI;


LI *LI_construct(FILE *fh);

void LI_destruct(LI *li);

/*
 * Check if the underlying file is indexable
 * Will be true for bgzipped and uncompressed, false for gzipped
 */
bool LI_indexable(LI *li);

/*
 * Get the next line from the file or NULL if at EOF.
 */
char *LI_get_next_line(LI *li);

/*
 * Peek at the next line, a call to peek will not get the next line from the file. In this way it can
 * be used to look ahead at the next line in the iteration sequence.
 */
char *LI_peek_at_next_line(LI *li);


/*
 * Go to position in file
 */
void LI_seek(LI *li, int64_t position);

/*
 * Tell the position in the file (ie where current line buf was read from)
 */
int64_t LI_tell(LI *li);


/*
 * Writer for maf and taf block and header writing
 */

typedef struct _LW {
    FILE *fh;
#ifdef USE_HTSLIB
    BGZF *bgzf;
#endif
} LW;

/*
 * Make a LW object. If use_compression is true and compiled with htslib will use bgzf compression on the stream.
 */
LW *LW_construct(FILE *fh, bool use_compression);

void LW_destruct(LW *lw, bool clean_up_file_handle);

int LW_write(LW *lw, const char *string, ...);

#endif /* STLINE_ITERATOR_H_ */

