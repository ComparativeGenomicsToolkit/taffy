/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibLogging.c
 *
 *  Created on: 24 May 2010
 *      Author: benedictpaten
 */

#define _POSIX_C_SOURCE 1 // needed for fileno()

#include "sonLibGlobalsInternal.h"
#include <errno.h>
// We need features.h to check if __GNU_LIBRARY__ is
// defined. Unfortunately, features.h is only available in glibc. The
// common workaround is to include limits.h, which includes features.h
// in glibc.
#include <limits.h>
#include <stdio.h>
#ifdef __GNU_LIBRARY__
#define USE_BACKTRACE 1
#define MAX_BACKTRACE_DEPTH 256
#include <execinfo.h>
#endif // __GLIBC__

static enum stLogLevel LOG_LEVEL = critical;

void *st_malloc(size_t i) {
    void *j;
    j = malloc(i);
    if (i > 0 && j == NULL) {
        st_errnoAbort("Malloc failed with a request for: %zu bytes", i);
    }
    return j;
}

void *st_realloc(void *buffer, size_t desiredSize) {
    void *newBuffer = realloc(buffer, desiredSize);
    if(desiredSize > 0 && newBuffer == NULL) {
        if(buffer != NULL) {
            free(buffer); // Free the old buffer
        }
        st_errnoAbort("Realloc failed with a request for: %zu bytes", desiredSize);
    }
    return newBuffer;
}

void *st_calloc(int64_t elementNumber, size_t elementSize) {
    void *k;
    k = calloc(elementNumber, elementSize);
    if (elementNumber > 0 && k == NULL) {
        st_errnoAbort("Calloc failed with request for %lld lots of %zu bytes",
                      (long long int) elementNumber, elementSize);
    }
    return k;
}

void st_setLogLevelFromString(const char *string) {
    if (string != NULL) {
        char *string2 = stString_copy(string);
        for (int64_t i = 0; i < strlen(string); i++) {
            string2[i] = tolower(string2[i]);
        }
        if (strcmp(string2, "off") == 0) {
            LOG_LEVEL = off;
        } else if (strcmp(string2, "critical") == 0) {
            LOG_LEVEL = critical;
        } else if (strcmp(string2, "error") == 0) {
            LOG_LEVEL = error;
        } else if (strcmp(string2, "warning") == 0) {
            LOG_LEVEL = warning;
        } else if (strcmp(string2, "info") == 0) {
            LOG_LEVEL = info;
        } else if (strcmp(string2, "debug") == 0) {
            LOG_LEVEL = debug;
        } else {
            st_errAbort("Unrecognised logging string %s", string);
        }
        free(string2);
    }
}

void st_setLogLevel(enum stLogLevel level) {
    LOG_LEVEL = level;
}

enum stLogLevel st_getLogLevel(void) {
    return LOG_LEVEL;
}

void st_logCritical(const char *string, ...) {
    if (st_getLogLevel() >= critical) {
        va_list ap;
        va_start(ap, string);
        vfprintf(stderr, string, ap);
        va_end(ap);
    }
}

void st_logInfo(const char *string, ...) {
    if (st_getLogLevel() >= info) {
        va_list ap;
        va_start(ap, string);
        vfprintf(stderr, string, ap);
        va_end(ap);
    }
}

void st_logDebug(const char *string, ...) {
    if (st_getLogLevel() >= debug) {
        va_list ap;
        va_start(ap, string);
        vfprintf(stderr, string, ap);
        va_end(ap);
    }
}

void st_uglyf(const char *string, ...) {
    va_list ap;
    va_start(ap, string);
    vfprintf(stderr, string, ap);
    va_end(ap);
}

int64_t st_system(const char *string, ...) {
    char cA[100000];
    int64_t i;
    va_list ap;
    va_start(ap, string);
    vsprintf(cA, string, ap);
    va_end(ap);
    assert(strlen(cA) < 100000);
    st_logDebug("Running command %s\n", cA);
    i = system(cA);
    //vfprintf(stdout, string, ap);
    return i;
}

static void st_backtrace(void) {
#ifdef USE_BACKTRACE
    void **backtraceArray = st_malloc(MAX_BACKTRACE_DEPTH * sizeof(void *));
    int depth = backtrace(backtraceArray, MAX_BACKTRACE_DEPTH);
    backtrace_symbols_fd(backtraceArray, depth, fileno(stderr));
    free(backtraceArray);
#endif // USE_BACKTRACE
}

void st_errAbort(char *format, ...) {
    st_backtrace();
    fprintf(stderr, "ERROR: ");
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

void st_errnoAbort(char *format, ...) {
    st_backtrace();
    fprintf(stderr, "ERROR: ");
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fprintf(stderr, ": %s\n", strerror(errno));
    exit(1);
}

static int64_t reverse8Bytes(int64_t in) {
    int64_t out;
    char *inByte = (char *) &in;
    char *outByte = (char *) &out;
    // looped for conciseness, hopefully it will be unrolled in
    // practice.
    for (int64_t i = 0; i < sizeof(int64_t); i++) {
        outByte[i] = inByte[sizeof(int64_t) - i];
    }
    return out;
}

int64_t st_nativeInt64FromLittleEndian(int64_t in) {
    static const int64_t one = 1;
    if (*(char *)&one == 0) {
        // big-endian
        return reverse8Bytes(in);
    } else {
        // little-endian
        return in;
    }
}

int64_t st_nativeInt64ToLittleEndian(int64_t in) {
    static const int64_t one = 1;
    if (*(char *)&one == 0) {
        // big-endian
        return reverse8Bytes(in);
    } else {
        // little-endian
        return in;
    }
}

int64_t st_nativeInt64ToBigEndian(int64_t in) {
    static const int64_t one = 1;
    if (*(char *)&one == 0) {
        // big-endian
        return in;
    } else {
        // little-endian
        return reverse8Bytes(in);
    }
}

int64_t st_nativeInt64FromBigEndian(int64_t in) {
    static const int64_t one = 1;
    if (*(char *)&one == 0) {
        // big-endian
        return in;
    } else {
        // little-endian
        return reverse8Bytes(in);
    }
}
