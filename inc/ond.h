#ifndef STOND_H_
#define STOND_H_

#include <inttypes.h>
#include <stdio.h>
#include <stdbool.h>

typedef struct _WFA WFA;

WFA *WFA_construct(void *string1, void *string2, int64_t string1_length, int64_t string2_length,
                   size_t element_size, bool (*elements_equal)(void *, void *),
                   int64_t gap_score, int64_t mismatch_score);

void WFA_destruct(WFA *wfa);

int64_t WFA_get_alignment_score(WFA *wfa);

void WFA_get_alignment(WFA *wfa, int64_t *elements_aligned_to_string1);

#endif /* STOND_H_ */

