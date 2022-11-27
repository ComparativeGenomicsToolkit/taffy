/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/**
 * Wrappers for C library functions that exit on errors.
 */
#include "stSafeC.h"
#include "sonLibExcept.h"
#include <stdbool.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

const char *ST_SAFEC_NUM_CONVERT_EXCEPTION_ID = "ST_SAFEC_NUM_CONVERT_EXCEPTION_ID";

/* Abort function that doesn't allocate any memory */
void stSafeCErr(const char *msg, ...)  {
    char buf[256];
    va_list args;
    va_start(args, msg);
    vsnprintf(buf, sizeof(buf), msg, args);
    va_end(args);
    // n.b.: check number of bytes here is a hack around not being able
    // to turn off gcc -Wunused-result for this special case
    int nbytes = write(2, buf, strlen(buf));
    nbytes += write(2, "\n", 1);
    _exit(128);
}

/* Allocate uninitialized memory. */
void *stSafeCMalloc(size_t size) {
    void *mem = malloc(size);
    if (mem == NULL) {
        stSafeCErr("can't allocate %zd bytes of memory", size);
    }
    return mem;
}

/* Allocate zeroed memory. */
void *stSafeCCalloc(size_t size) {
    void *mem = stSafeCMalloc(size);
    memset(mem, 0, size);
    return mem;
}

/* Reallocated memory. */
void *stSafeCRealloc(void *mem, size_t size) {
    if ((mem = realloc(mem, size)) == NULL) {
        stSafeCErr("can't reallocate %zd bytes of memory", size);
    }
    return mem;
}

/* copy a block of memory */
void *stSafeCCopyMem(void *mem, size_t size) {
    void *mem2 = stSafeCMalloc(size);
    memcpy(mem2, mem, size);
    return mem2;
}

/* sprintf format with buffer overflow checking. */
int stSafeCFmtv(char *buffer, int bufSize, const char *format, va_list args) {
    int sz = vsnprintf(buffer, bufSize, format, args);
    /* note that some version return -1 if too small */
    if ((sz < 0) || (sz >= bufSize))
        stSafeCErr("buffer overflow, size %d, format: %s", bufSize, format);
    return sz;
}

/* sprintf format with buffer overflow checking. */
int stSafeCFmt(char* buffer, int bufSize, const char *format, ...) {
    int sz;
    va_list args;
    va_start(args, format);
    sz = stSafeCFmtv(buffer, bufSize, format, args);
    va_end(args);
    return sz;
}

/* sprintf formatting, returning a dynamically allocated string. */
char *stSafeCDynFmtv(const char *format, va_list args) {
    int bufSize = 2*strlen(format);
    char *buf = stSafeCMalloc(bufSize); 
    while (true) {
        va_list argscp;
        va_copy(argscp, args);
        int sz = vsnprintf(buf, bufSize, format, argscp);
        va_end(argscp);

        /* note that some version return -1 if too small */
        if ((sz < 0) || (sz >= bufSize)) {
            bufSize *= 2;
            buf = stSafeCRealloc(buf, bufSize);
        } else {
            break;
        }
    }
    return buf;
}

/* sprintf formatting, returning a dynamically allocated string. */
char *stSafeCDynFmt(const char *format, ...) {
    va_list args;
    va_start(args, format);
    char *buf = stSafeCDynFmtv(format, args);
    va_end(args);
    return buf;
}


/* convert a string to a 32 unsigned int, exception if invalid */
uint64_t stSafeStrToUInt32(const char *str) {
    char *end;
    long num = strtol(str, &end, 10);
    if ((end == str) || (*end != '\0') || (num < 0)) {
        stThrowNew(ST_SAFEC_NUM_CONVERT_EXCEPTION_ID, "invalid 32-bit unsigned integer: \"%s\"", str);
    }
    return num;
}

/* convert a string to a 64 int, exception if invalid */
int64_t stSafeStrToInt64(const char *str) {
    char *end;
    long long num = strtoll(str, &end, 10);
    if ((end == str) || (*end != '\0')) {
        stThrowNew(ST_SAFEC_NUM_CONVERT_EXCEPTION_ID, "invalid 64-bit signed integer: \"%s\"", str);
    }
    return num;
}
