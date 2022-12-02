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

typedef struct _LI LI;

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

#endif /* STLINE_ITERATOR_H_ */

