/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/**
 * stExcept - a C exception mechanism used by SonTrace.
 *
 * This provides a simple setjmp based exception mechanism for C.  All errors
 * in the API are reported via exceptions.  An exception object, containing
 * both a symbolic and natural language error description is used to represent
 * errors.  Exception chaining is supported.
 *
 * The usage idiom for catching errors is:
 * \code
 *  stTry {
 *      // user code here
 *  } stCatch(except) {
 *      // user exception handling here
 *      stExcept_free(except);
 *  } stTryEnd;
 * \endcode
 *
 * To return within the try block, \code stTryReturn(value) \endcode must
 * be used.  There is no macro to return within the try block without
 * returning a value, as this is considered tacky. A standard
 * \code return \endcode can be used in the catch block.
 *
 * Any variable who's scope is outside of the try block and is modified
 * within the try block must be declared \code volatile \endcode.  For
 * example:
 * \code
 *  volatile int cnt = 0;
 *  ...
 *  stTry {
 *      cnt++;
 *      if (cnt > 10) {
 *         stTryReturn(cnt);
 *      }
 *  } stCatch(except) {
 *      stExcept_free(except);
 *      return -cnt;
 *  } stTryEnd;
 *  return cnt;
 * \endcode
 *
 * If environment variable ST_ABORT is set, throwing an error causes the
 * message to be printed to stderr and an abort(), which is useful for
 * stopping under a debugger. Setting ST_ABORT_UNCAUGHT environment
 * causes abort only on uncaught exceptions.
 *
 * @defgroup CExceptions Exceptions for C
 */

#ifndef sonLibExcept_h
#define sonLibExcept_h
#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <setjmp.h>
#include <stdarg.h>
#include <string.h>
//@{
/// @defgroup CExcept class CExcept
/// @ingroup CExceptions
//@{

/**
 * Construct a new stExcept object.
 * 
 * @param id symbolic exception id.  This should be a constant, static string.
 * @param msg print style format used to generate errMsg
 * @param args arguments to format into errMsg
 * @ingroup stExcept
 */
stExcept *stExcept_newv(const char *id, const char *msg, va_list args);

/**
 * Construct a new stExcept object.
 * 
 * @param id symbolic exception id.  This should be a constant, static string.
 * @param msg print style format used to generate errMsg
 * @param ... arguments to format into errMsg
 * @ingroup stExcept
 */
stExcept *stExcept_new(const char *id, const char *msg, ...);

/**
 * Construct a new stExcept object, setting cause.
 * 
 * @param cause causing error cause, ownership passed to the new object.
 * @param id symbolic exception id.  This should be a constant, static string.
 * @param msg print style format used to generate errMsg
 * @param args arguments to format into errMsg
 * @ingroup stExcept
 */
stExcept *stExcept_newCausev(stExcept *cause, const char *id, const char *msg, va_list args);

/**
 * Construct a new stExcept object, setting cause.
 * 
 * @param cause causing error cause, ownership passed to the new object.
 * @param id symbolic exception id.  This should be a constant, static string.
 * @param msg print style format used to generate errMsg
 * @param ... arguments to format into errMsg
 * @ingroup stExcept
 */
stExcept *stExcept_newCause(stExcept *cause, const char *id, const char *msg, ...);

/**
 * Free an stExcept object.  Frees assocated causes as well.
 * @ingroup stExcept
 */
void stExcept_free(stExcept *except);

/**
 * Get the id for a stExcept.
 * @ingroup stExcept
 */
const char* stExcept_getId(const stExcept *except);

/**
 * Determine if the id in an exception is equal to a specified id.
 * @ingroup stExcept
 */
static inline bool stExcept_idEq(const stExcept *except, const char *id) {
    return strcmp(stExcept_getId(except), id) == 0;
}

/**
 * Get the message for a stExcept.
 * @ingroup stExcept
 */
const char* stExcept_getMsg(const stExcept *except);

/**
 * Get the cause for a stExcept.
 * @ingroup stExcept
 */
stExcept *stExcept_getCause(const stExcept *except);
//@}

/*
 * Context allocated on the stack for an exception.
 * (Internal structure, don't use directly)
 */
struct _stExceptContext {
    jmp_buf env;                     // setjmp environment 
    struct _stExceptContext *prev;   // previous context on the stack.
    stExcept *except;                // exception thrown into this context
};

/* 
 * Exception content Top Of Stack.
 * (Internal structure, don't use directly)
 */
extern struct _stExceptContext *_cexceptTOS;

/// @defgroup CMacros C try/catch macros
/// @ingroup stExceptions
//@{

/**
 * Begin a try block.
 */
#define stTry {\
    struct _stExceptContext _cexceptContext;\
    _cexceptContext.prev = _cexceptTOS;\
    _cexceptTOS = &_cexceptContext;\
    _cexceptContext.except = NULL;\
    if (setjmp(_cexceptContext.env) == 0)
    
/**
 * Catch following a try block.
 */
#define stCatch(exceptVar) \
    if (_cexceptTOS->except != NULL) {\
        stExcept *exceptVar = _cexceptTOS->except;\
        _cexceptTOS = _cexceptTOS->prev;

/**
 * End of try/cache block
 */
#define stTryEnd \
        stExcept_free(_cexceptContext.except);\
    } else {\
        _cexceptTOS = _cexceptTOS->prev;\
    }}

/**
 * Return a value inside of a ceTry.
 */
#define stTryReturn(val) {\
    _cexceptTOS = _cexceptTOS->prev;\
    return val;\
    }

/**
 * Raise an exception.  This function is not inlined to allow setting
 * breakpoints.
 * @param except object describing the exception.  Ownership of object
 *  is passed to exception mechanism.
 */
void stThrow(stExcept *except);

/**
 * Construct and raise an exception.
 * @param id symbolic exception id.  This should be a constant, static string.
 * @param msg print style format used to generate errMsg
 * @param ... arguments to format into errMsg
 */
void stThrowNew(const char *id, const char *msg, ...);

/**
 * Construct and raise an exception, setting cause.
 * @param cause causing error cause, ownership passed to the new object.
 * @param id symbolic exception id.  This should be a constant, static string.
 * @param msg print style format used to generate errMsg
 * @param ... arguments to format into errMsg
 */
void stThrowNewCause(stExcept *cause, const char *id, const char *msg, ...);
//@}
//@}
#ifdef __cplusplus
}
#endif
#endif
