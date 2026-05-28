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
    // Coalescing byte buffer for the LW_put* family.  Pending bytes here
    // MUST be flushed before any LW_write / bgzf_write / fprintf so output
    // order is preserved.  Buffer grows as needed up to LW_BUF_TARGET; we
    // flush to the underlying stream when at least that much is pending.
    char  *buf;
    size_t buf_pos;
    size_t buf_cap;
} LW;

/*
 * Make a LW object. If use_compression is true and compiled with htslib will use bgzf compression on the stream.
 */
LW *LW_construct(FILE *fh, bool use_compression);

void LW_destruct(LW *lw, bool clean_up_file_handle);

int LW_write(LW *lw, const char *string, ...);

/*
 * Byte-level emitters that bypass vsnprintf.  Bytes accumulate in an
 * internal buffer and flush when full or at LW_flush / LW_write /
 * LW_destruct.  Use these on hot emit paths (TAF column writer, etc.);
 * LW_write remains correct for headers and other cold paths.
 *
 * `LW_puti64` writes a decimal int64 (no padding) using a hand-rolled
 * two-digits-per-iteration table -- much cheaper than vsnprintf("%"PRIi64).
 *
 * `LW_putc` is inlined: the fast path is one branch + one store, and the
 * full-buffer slow path tail-calls out to LW_putc_slow which flushes and
 * appends.  Each emitted byte going through a function call shows up in
 * callgrind on TAF writes, so the inline matters.
 */
void LW_putc_slow(LW *lw, char c);
static inline void LW_putc(LW *lw, char c) {
    if (lw->buf_pos < lw->buf_cap) {
        lw->buf[lw->buf_pos++] = c;
    } else {
        LW_putc_slow(lw, c);
    }
}
void LW_puts(LW *lw, const char *s);
void LW_putn(LW *lw, const char *s, size_t n);
/* Append `c` repeated `n` times.  Replaces a hot `for(n) LW_putc(c)` loop
 * in the unencoded-bases path of the TAF column writer; one memset is
 * dramatically cheaper than n function calls + branches. */
void LW_putrep(LW *lw, char c, size_t n);
void LW_puti64(LW *lw, int64_t v);
void LW_flush(LW *lw);

#endif /* STLINE_ITERATOR_H_ */

