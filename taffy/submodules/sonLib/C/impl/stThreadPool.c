// A very simple pthreads thread pool implementation.
//
// A work stack feeds a number of threads. Every work unit is worked
// on by a single thread, which executes the work function, then
// fetches a new work unit when finished.
//
// Often threads only need to write once to some single "result"
// variable, after the heavy computational work is done. To avoid
// needing to manage your own locks for this simple case, an optional
// "finish" function is available, which is guaranteed to run serially
// (although the work function and main thread will still run
// alongside it). If you don't want to use a finishing function, then
// the value returned by the work function is completely ignored.
//
// Note that the work "queue" is actually a stack, not a queue. This
// shouldn't matter, as you have absolutely no ordering guarantees on
// what order the work is completed in regardless.
#include <stdlib.h>
#include <pthread.h>
#include "sonLib.h"

struct _stThreadPool {
    pthread_mutex_t stackLock;  // Locks the stack so work can be
                                // pushed/popped serially.

    pthread_cond_t stackCond;   // Signals threads waiting on the
                                // stack that there's new data. Locked
                                // by stackLock.

    int64_t numFinishedThreads; // Number of threads currently waiting
                                // on the stack. Locked by stackLock
                                // as well.

    pthread_cond_t finishedCond;// Signals the parent thread that work
                                // just finished. Locked by stackLock.

    stList *stack;              // Work units stored here. Locked by
                                // stackLock.

    pthread_mutex_t finishLock; // Lock to ensure that all "finish" functions
                                // are done serially.

    pthread_t *threads;         // The threads in the pool.

    int64_t numThreads;         // Number of threads in the pool.

    void *(*workFunc)(void *);  // Function that takes in a work unit
                                // and returns a result.

    void (*finishFunc)(void *); // Function that takes a result and
                                // does something with it.

    bool killFlag;              // Set when the thread pool wants the
                                // threads to stop. pthread_cancel
                                // requires that they hit a
                                // cancellation point, which isn't
                                // guaranteed to happen.
};

// Worker function for each thread spawned. This function should never
// return until it is cancelled in the destructor.
static void worker(stThreadPool *threadPool) {
    int oldtype;
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &oldtype);
    // We don't need the old cancellation type but we have to specify
    // a non-NULL argument to be POSIX-compatible
    (void) oldtype;
    while (1) {
        pthread_mutex_lock(&threadPool->stackLock);
        // Signal that we're waiting for work.
        threadPool->numFinishedThreads++;
        while (stList_length(threadPool->stack) == 0 && !threadPool->killFlag) {
            // Wake up the main thread in case it's waiting for us to be
            // done.
            pthread_cond_signal(&threadPool->finishedCond);
            // This unlocks the lock and waits to be signaled on stackCond.
            pthread_cond_wait(&threadPool->stackCond, &threadPool->stackLock);
            // We have the lock again now.
        }

        if (threadPool->killFlag) {
            // Time to go away now
            pthread_mutex_unlock(&threadPool->stackLock);
            pthread_exit(NULL);
        }

        void *workUnit = stList_pop(threadPool->stack);
        // Before unlocking the lock, un-signal that we're waiting for work.
        threadPool->numFinishedThreads--;
        pthread_mutex_unlock(&threadPool->stackLock);

        // Do the actual work.
        void *result = threadPool->workFunc(workUnit);

        // Finish up serially if needed.
        if (threadPool->finishFunc != NULL) {
            pthread_mutex_lock(&threadPool->finishLock);
            threadPool->finishFunc(result);
            pthread_mutex_unlock(&threadPool->finishLock);
        }
    }
}

// Initialize the thread pool. finishFunc can be NULL if you are
// managing your own locks for output.
stThreadPool *stThreadPool_construct(int64_t numThreads,
                                     void *(*workFunc)(void *),
                                     void (*finishFunc)(void *)) {
    assert(numThreads > 0);
    stThreadPool *ret = st_calloc(1, sizeof(stThreadPool));
    ret->threads = st_calloc(numThreads, sizeof(pthread_t));
    ret->stack = stList_construct();
    ret->workFunc = workFunc;
    ret->finishFunc = finishFunc;
    ret->numThreads = numThreads;

    // Set up the locks.
    int pthreadError;
    pthreadError = pthread_mutex_init(&ret->stackLock, NULL);
    if (pthreadError) {
        st_errAbort("stThreadPool: pthread_mutex_init failed: %s",
                    strerror(pthreadError));
    }
    pthreadError = pthread_mutex_init(&ret->finishLock, NULL);
    if (pthreadError) {
        st_errAbort("stThreadPool: pthread_mutex_init failed: %s",
                    strerror(pthreadError));
    }
    pthreadError = pthread_cond_init(&ret->stackCond, NULL);
    if (pthreadError) {
        st_errAbort("stThreadPool: pthread_cond_init failed: %s",
                    strerror(pthreadError));
    }
    pthreadError = pthread_cond_init(&ret->finishedCond, NULL);
    if (pthreadError) {
        st_errAbort("stThreadPool: pthread_cond_init failed: %s",
                    strerror(pthreadError));
    }

    // Start the threads. All initialization of the thread pool struct
    // must happen before this, as the threads will look for work
    // right away.
    for (int64_t i = 0; i < numThreads; i++) {
        pthreadError = pthread_create(&ret->threads[i], NULL,
                                      (void *(*)(void *)) worker, ret);
        if (pthreadError) {
            st_errAbort("stThreadPool: pthread_create failed: %s",
                        strerror(pthreadError));
        }
    }
    return ret;
}

// Push work onto the stack to be consumed by a thread.
void stThreadPool_push(stThreadPool *threadPool, void *workUnit) {
    pthread_mutex_lock(&threadPool->stackLock);
    stList_append(threadPool->stack, workUnit);
    pthread_cond_signal(&threadPool->stackCond);
    pthread_mutex_unlock(&threadPool->stackLock);
}

// Block until all work currently in the stack is complete. Can block
// indefinitely if something goes wrong. Use stThreadPool_waitSafe to
// ensure that you can get execution back after a timeout.
void stThreadPool_wait(stThreadPool *threadPool) {
    pthread_mutex_lock(&threadPool->stackLock);
    while (stList_length(threadPool->stack) != 0 || threadPool->numFinishedThreads != threadPool->numThreads) {
        pthread_cond_wait(&threadPool->finishedCond, &threadPool->stackLock);
    }
    pthread_mutex_unlock(&threadPool->stackLock);
}

// Doesn't wait on all the work to complete, just returns whether the
// work in the stack is all done or not.
bool stThreadPool_done(stThreadPool *threadPool) {
    pthread_mutex_lock(&threadPool->stackLock);
    bool ret = stList_length(threadPool->stack) == 0 && threadPool->numFinishedThreads == threadPool->numThreads;
    pthread_mutex_unlock(&threadPool->stackLock);
    return ret;
}

// Destroys the thread pool and destroys all threads. They must finish
// running their work before they exit, but any remaining work in the
// queue will be unfinished.
void stThreadPool_destruct(stThreadPool *threadPool) {
    // Wake all currently running threads so they know that they need
    // to die.
    threadPool->killFlag = true;
    for (int64_t i = 0; i < threadPool->numThreads; i++) {
        pthread_cond_signal(&threadPool->stackCond);
    }
    // Ensure that all threads are dead before freeing the memory out
    // from under them.
    for (int64_t i = 0; i < threadPool->numThreads; i++) {
        pthread_join(threadPool->threads[i], NULL);
    }
    free(threadPool->threads);

    stList_destruct(threadPool->stack);

    pthread_mutex_destroy(&threadPool->stackLock);
    pthread_mutex_destroy(&threadPool->finishLock);
    pthread_cond_destroy(&threadPool->finishedCond);
    pthread_cond_destroy(&threadPool->stackCond);

    free(threadPool);
}
