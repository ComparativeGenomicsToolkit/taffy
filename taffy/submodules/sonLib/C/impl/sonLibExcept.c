/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLibExcept.h"
#include <stdlib.h>
#include "stSafeC.h"

/**
 * Exception object that is thrown.
 * @ingroup stExcept
 */
struct stExcept {
    const char *id;          /**< symbolic exception id */
    const char *msg;         /**< natural-language error message */
    struct stExcept *cause;  /**< error stack, NULL if no causing errors */
};

/*
 * Exception content Top Of Stack.
 * FIXME: should be in thread-local storage.
 */
struct _stExceptContext *_cexceptTOS = NULL;

stExcept *stExcept_newv(const char *id, const char *msg, va_list args) {
    stExcept *except = stSafeCCalloc(sizeof(stExcept));
    except->id = id;
    except->msg = stSafeCDynFmtv(msg, args);
    return except;
}

stExcept *stExcept_new(const char *id, const char *msg, ...) {
    va_list args;
    va_start(args, msg);
    stExcept *except = stExcept_newv(id, msg, args);
    va_end(args);
    return except;
}

stExcept *stExcept_newCausev(stExcept *cause, const char *id, const char *msg, va_list args) {
    stExcept *except = stExcept_newv(id, msg, args);
    except->cause = cause;
    return except;
}

stExcept *stExcept_newCause(stExcept *cause, const char *id, const char *msg, ...) {
    va_list args;
    va_start(args, msg);
    stExcept *except = stExcept_newCausev(cause, id, msg, args);
    va_end(args);
    return except;
}

void stExcept_free(stExcept *except) {
    if (except != NULL) {
        if (except->cause != NULL) {
            stExcept_free(except->cause);
        }
        stSafeCFree((char*)except->msg);
        stSafeCFree(except);
    }
}

const char* stExcept_getId(const stExcept *except) {
    return except->id;
}

const char* stExcept_getMsg(const stExcept *except) {
    return except->msg;
}

stExcept *stExcept_getCause(const stExcept *except) {
    return except->cause;
}

static void handleUncaught(stExcept *except) {
    fflush(stdout);
    fprintf(stderr, "Exception: %s: %s\n", except->id, except->msg);
    for (stExcept *cause = except->cause; cause != NULL; cause = cause->cause) {
        fprintf(stderr, "\tcaused by: %s: %s\n", cause->id, cause->msg);
    }
    fflush(stderr);
    if (getenv("ST_ABORT_UNCAUGHT") != NULL) {
        abort();
    }
    stSafeCErr("Uncaught exception");
}

void stThrow(stExcept *except) {
    assert(except != NULL);
    if (getenv("ST_ABORT") != NULL) {
        // print message and call abort to enter debugger
        fprintf(stderr, "Exception: %s: %s\n", except->id, except->msg);
        abort();
    }

    if (_cexceptTOS == NULL) {
        handleUncaught(except);
    }
    _cexceptTOS->except = except;
    longjmp(_cexceptTOS->env, 1);
}

void stThrowNew(const char *id, const char *msg, ...) {
    va_list args;
    va_start(args, msg);
    stExcept *except = stExcept_newv(id, msg, args);
    va_end(args);
    stThrow(except);
}

void stThrowNewCause(stExcept *cause, const char *id, const char *msg, ...) {
    va_list args;
    va_start(args, msg);
    stExcept *except = stExcept_newCausev(cause, id, msg, args);
    va_end(args);
    stThrow(except);
}

