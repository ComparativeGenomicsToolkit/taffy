/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLibGlobalsTest.h"

static char *tempFileDir = "sonLibFileTestTempDir";
static char *tempFileName1 =
        "sonLibFileTestTempDir/sonLibFileTestTempFile1.txt";
static char *tempFileName2 =
        "sonLibFileTestTempDir/sonLibFileTestTempFile2.txt";

static void teardown() {
    if (stFile_exists(tempFileDir)) {
        stFile_rmrf(tempFileDir);
    }
}

static void setup() {
    teardown();

    //Make the temp dir
    stFile_mkdir(tempFileDir);
    //Now the containing files

    //Temp file 1
    FILE *fileHandle = st_fopen(tempFileName1, "w");
    fprintf(fileHandle, "hello world\n");
    fprintf(fileHandle, "foo bar 123456\n");
    fprintf(fileHandle, " \n");
    fprintf(fileHandle, "\n");
    fprintf(fileHandle, "bye bye\n");
    fprintf(fileHandle, "\t");
    fclose(fileHandle);

    //Temp file 2
    fileHandle = st_fopen(tempFileName2, "w");
    fclose(fileHandle);
}

static void test_stFile_getLineFromFile(CuTest *testCase) {
    setup();
    FILE *fileHandle = st_fopen(tempFileName1, "r");
    char *s = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "hello world", s);
    free(s);
    s = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "foo bar 123456", s);
    free(s);
    s = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, " ", s);
    free(s);
    s = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "", s);
    free(s);
    s = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "bye bye", s);
    free(s);
    s = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "\t", s);
    free(s);
    CuAssertTrue(testCase, NULL == stFile_getLineFromFile(fileHandle));
    CuAssertTrue(testCase, NULL == stFile_getLineFromFile(fileHandle));
    CuAssertTrue(testCase, NULL == stFile_getLineFromFile(fileHandle));
    CuAssertTrue(testCase, NULL == stFile_getLineFromFile(fileHandle));
    CuAssertTrue(testCase, NULL == stFile_getLineFromFile(fileHandle));
    fclose(fileHandle);
    teardown();
}

static void test_stFile_exists(CuTest *testCase) {
    teardown();
    CuAssertTrue(testCase, !stFile_exists(tempFileDir));
    CuAssertTrue(testCase, !stFile_exists(tempFileName1));
    CuAssertTrue(testCase, !stFile_exists(tempFileName2));
    setup();
    CuAssertTrue(testCase, stFile_exists(tempFileDir));
    CuAssertTrue(testCase, stFile_exists(tempFileName1));
    CuAssertTrue(testCase, stFile_exists(tempFileName2));
    teardown();
}

static void test_stFile_isDir(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stFile_isDir(tempFileDir));
    CuAssertTrue(testCase, !stFile_isDir(tempFileName1));
    CuAssertTrue(testCase, !stFile_isDir(tempFileName2));
    teardown();
    stTry {
            stFile_isDir(tempFileDir);
            CuAssertTrue(testCase, 0);
        }
        stCatch(except)
            {
                CuAssertStrEquals(testCase, stExcept_getId(except),
                        ST_FILE_EXCEPTION);
            }stTryEnd;
}

static void test_stFile_getFileNamesInDirectory(CuTest *testCase) {
    setup();
    stList *childFiles = stFile_getFileNamesInDirectory(tempFileDir);
    CuAssertIntEquals(testCase, stList_length(childFiles), 2);
    char *file1 = stFile_pathJoin(tempFileDir, stList_get(childFiles, 0));
    char *file2 = stFile_pathJoin(tempFileDir, stList_get(childFiles, 1));
    CuAssertTrue(testCase, stString_eq(file1, tempFileName1) || stString_eq(file1, tempFileName2));
    CuAssertTrue(testCase, stString_eq(file2, tempFileName1) || stString_eq(file2, tempFileName2));
    CuAssertTrue(testCase, !stString_eq(file1, file2));
    stList_destruct(childFiles);
    free(file1);
    free(file2);
    stTry {
            stFile_getFileNamesInDirectory(tempFileName1);
            CuAssertTrue(testCase, 0);
        }
        stCatch(except)
            {
                CuAssertStrEquals(testCase, stExcept_getId(except),
                        ST_FILE_EXCEPTION);
            }stTryEnd;
    teardown();
    stTry {
            stFile_getFileNamesInDirectory(tempFileDir);
            CuAssertTrue(testCase, 0);
        }
        stCatch(except)
            {
                CuAssertStrEquals(testCase, stExcept_getId(except),
                        ST_FILE_EXCEPTION);
            }stTryEnd;
}

static void test_stFile_pathJoin(CuTest *testCase) {
    char *s = stFile_pathJoin("hello", "world");
    CuAssertStrEquals(testCase, "hello/world", s);
    free(s);
    s = stFile_pathJoin("hello/", "world");
    CuAssertStrEquals(testCase, "hello/world", s);
    free(s);
    s = stFile_pathJoin("", "world");
    CuAssertStrEquals(testCase, "world", s);
    free(s);
    s = stFile_pathJoin("hello", "");
    CuAssertStrEquals(testCase, "hello/", s);
    free(s);
    s = stFile_pathJoin("", "");
    CuAssertStrEquals(testCase, "", s);
    free(s);
}

CuSuite* sonLibFileTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stFile_getLineFromFile);
    SUITE_ADD_TEST(suite, test_stFile_pathJoin);
    SUITE_ADD_TEST(suite, test_stFile_exists);
    SUITE_ADD_TEST(suite, test_stFile_isDir);
    SUITE_ADD_TEST(suite, test_stFile_getFileNamesInDirectory);

    return suite;
}
