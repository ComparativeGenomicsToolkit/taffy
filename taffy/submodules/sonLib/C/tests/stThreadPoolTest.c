#include "CuTest.h"
#include "sonLib.h"

// Test sorting a few sublists into a larger list to try to catch out
// any race conditions.
static stList *sorted;

static stList *sortSubList(stList *input) {
    stList_sort(input, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
    return input;
}

// Ideally this would be a true merge sort. But since this is just to
// test for race conditions we do something slower.
static void insertSubList(stList *result) {
    int64_t result_i = 0, sorted_i = 0;
    stList *newSorted = stList_construct();
    while (result_i != stList_length(result) || sorted_i != stList_length(sorted)) {
        if (result_i != stList_length(result) && (sorted_i == stList_length(sorted) || stIntTuple_cmpFn(stList_get(result, result_i), stList_get(sorted, sorted_i)) == -1)) {
            stList_append(newSorted, stList_get(result, result_i));
            result_i++;
        } else if (sorted_i != stList_length(sorted)) {
            stList_append(newSorted, stList_get(sorted, sorted_i));
            sorted_i++;
        }
    }

    stList_destruct(sorted);
    sorted = newSorted;
}

static void testStThreadPoolSort(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 3; testNum++) {
        // Create two lists, one containing a bunch of random numbers
        // and the other containing a bunch of sublists that contain
        // the same numbers, in aggregate.
        stList *truth = stList_construct();
        sorted = stList_construct();
        // This list and its sublists are the only ones that own the
        // stIntTuples.
        stList *lists = stList_construct3(0, (void (*)(void *)) stList_destruct);
        int64_t numSubLists = st_randomInt64(1, 100);
        for (int64_t i = 0; i < numSubLists; i++) {
            stList *subList = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
            int64_t numNums = st_randomInt64(0, 50000);
            for (int64_t j = 0; j < numNums; j++) {
                stList_append(subList, stIntTuple_construct1(st_randomInt64(-100000, 100000)));
            }
            stList_appendAll(truth, subList);
            stList_append(lists, subList);
        }

        // First sort our list.
        stList_sort(truth, (int (*)(const void *, const void *)) stIntTuple_cmpFn);

        // Now create a bunch of threads and ask them to sort theirs.
        stThreadPool *threadPool = stThreadPool_construct(st_randomInt64(1, 6),
                                                          (void *(*)(void *)) sortSubList,
                                                          (void (*)(void *)) insertSubList);
        for (int64_t i = 0; i < stList_length(lists); i++) {
            stThreadPool_push(threadPool, stList_get(lists, i));
        }

        // Wait for the process to complete.
        stThreadPool_wait(threadPool);
        CuAssertTrue(testCase, stThreadPool_done(threadPool));

        CuAssertIntEquals(testCase, stList_length(truth), stList_length(sorted));
        for (int64_t i = 0; i < stList_length(sorted); i++) {
            CuAssertIntEquals(testCase,
                              stIntTuple_get(stList_get(sorted, i), 0),
                              stIntTuple_get(stList_get(truth, i), 0));
        }

        stThreadPool_destruct(threadPool);
        stList_destruct(sorted);
        stList_destruct(lists);
        stList_destruct(truth);
    }
}

CuSuite *sonLib_stThreadPoolTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStThreadPoolSort);
    return suite;
}
