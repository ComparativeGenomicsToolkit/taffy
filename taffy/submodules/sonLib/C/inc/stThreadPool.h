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
#ifndef SONLIB_THREADPOOL_H_
#define SONLIB_THREADPOOL_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef struct _stThreadPool stThreadPool;

// Initialize the thread pool. finishFunc can be NULL if you are
// managing your own locks for output.
stThreadPool *stThreadPool_construct(int64_t numThreads,
                                     void *(*workFunc)(void *),
                                     void (*finishFunc)(void *));

// Push work onto the stack to be consumed by a thread.
void stThreadPool_push(stThreadPool *threadPool, void *workUnit);

// Block until all work currently in the stack is complete. Can block
// indefinitely if something goes wrong. Use stThreadPool_waitSafe to
// ensure that you can get execution back after a timeout.
void stThreadPool_wait(stThreadPool *threadPool);

// Doesn't wait on all the work to complete, just returns whether the
// work in the stack is all done or not.
bool stThreadPool_done(stThreadPool *threadPool);

// Destroys the thread pool and destroys all threads *without*
// checking if they are still processing work or if there is work in
// the queue.
void stThreadPool_destruct(stThreadPool *threadPool);

#ifdef __cplusplus
}
#endif
#endif // SONLIB_THREADPOOL_H_
