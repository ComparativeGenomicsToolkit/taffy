/*
 * Released under the MIT license, see LICENSE.txt
 * 
 * Wrapper to create a thread-safe sort. Deals with the fact 
 * that qsort_r has different signature on OS/X and Linux
 *
 * WARNING: this should be include first, or at least before the
 * include of stdlib.h, otherwise the qsort_r function will not be
 * defined on Linux.
 */
#ifndef SAFESORT_H
#define SAFESORT_H

#ifdef __linux__
/* linux version */
#define _GNU_SOURCE
#include <stdlib.h>
static inline
void safesort(void* base, size_t nmemb, size_t size,
              int (*compar)(const void *, const void *, void *arg),
              void *arg) {
    qsort_r(base, nmemb, size, compar, arg);
      
}
#else
/* OS/X, probably BSDs */
#include <stdlib.h>
/* this swaps order of args */
struct _safesort_arg {
    int (*compar)(const void *, const void *, void *arg);
    void* realArg;
};

static int _safesort_compar(void *arg, const void *el1, const void *el2) {
    struct _safesort_arg *sarg = arg;
    return sarg->compar(el1, el2, sarg->realArg);
}

static inline
void safesort(void* base, size_t nmemb, size_t size,
              int (*compar)(const void *, const void *, void *arg),
              void *arg) {
    struct _safesort_arg sarg = {compar, arg};
    qsort_r(base, nmemb, size, &sarg, _safesort_compar);
}
#endif

#endif
