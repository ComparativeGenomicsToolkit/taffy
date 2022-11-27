/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/**
 * Wrappers for C library functions that exit on errors.
 * @defgroup stSafeC Robust C library functions
 */
#ifndef stSafeC_h
#define stSafeC_h
#include <stdlib.h>
#include <stdarg.h>
#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

//@{

/**
 * Exception id for numeric conversion errors.
 * @ingroup
 */
//const char *ST_SAFEC_NUM_CONVERT_EXCEPTION_ID;

/**
 * Abort function that doesn't allocate any memory
 * @ingroup safec
 */
void stSafeCErr(const char *msg, ...)
#if defined(__GNUC__)
__attribute__((format(printf, 1, 2)))
#endif
;

/**
 * Allocate uninitialized memory, exiting with using minimal resources if it
 * can't be allocated. Use stSafeCCalloc for cleared memory.
 * @ingroup stSafeC
 */
void *stSafeCMalloc(size_t size);

/**
 * Allocate zeroed memory, exiting with using minimal resources if it can't be
 * allocated.  Should be used for small objects
 * @ingroup stSafeC
 */
void *stSafeCCalloc(size_t size);

/**
 * Reallocated memory, exiting with using minimal resources if it can't be
 * allocated.
 * @ingroup stSafeC
 */
void *stSafeCRealloc(void *mem, size_t size);

/**
 * wrapper around free, to be consistent with other alloc functions.
 * @ingroup stSafeC
 */
static inline void stSafeCFree(void *mem) {
    if (mem != NULL) {
        free(mem);
    }
}

/* copy a block of memory */
void *stSafeCCopyMem(void *mem, size_t size);

/**
 * sprintf format with buffer overflow checking.  The resulting string is
 * always terminated with zero byte.
 * @ingroup stSafeC
 */
int stSafeCFmtv(char *buffer, int bufSize, const char *format, va_list args);

/**
 * sprintf format with buffer overflow checking.  The resulting string is
 * always terminated with zero byte.
 * @ingroup stSafeC
 */
int stSafeCFmt(char* buffer, int bufSize, const char *format, ...)
#if defined(__GNUC__)
__attribute__((format(printf, 3, 4)))
#endif
;

/**
 * sprintf formatting, returning a dynamically allocated string.
 * @ingroup stSafeC
 */
char *stSafeCDynFmtv(const char *format, va_list args);

/**
 * sprintf formatting, returning a dynamically allocated string.
 * @ingroup stSafeC
 */
char *stSafeCDynFmt(const char *format, ...)
#if defined(__GNUC__)
__attribute__((format(printf, 1, 2)))
#endif
;

/* convert a string to a 32 unsigned int, exception if invalid */
uint64_t stSafeStrToUInt32(const char *str);

/* convert a string to a 64 int, exception if invalid */
int64_t stSafeStrToInt64(const char *str);

#ifdef __cplusplus
}
#endif
#endif
