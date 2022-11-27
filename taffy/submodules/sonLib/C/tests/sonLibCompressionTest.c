/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibRandomTest.c
 *
 *  Created on: 22-Jun-2010
 *      Author: benedictpaten
 */

#include "sonLibGlobalsTest.h"
#include <time.h>

static void test_stCompression_compressAndDecompressP(CuTest *testCase, int64_t rounds, int64_t minSize, int64_t maxSize,
        void *(*compress)(void *, int64_t, int64_t *, int64_t), void *(*decompress)(void *, int64_t, int64_t *)) {
    /*
     * Exercises the two compression functions.
     */
    int64_t bytesCompressed = 0, bytesDecompressed = 0;
    for(int64_t i=0; i<rounds; i++) {
        int64_t size = st_randomInt64(minSize, maxSize);
        char *randomString = st_malloc(size);
        for(int64_t j=0; j<size; j++) {
            randomString[j] = (char)st_randomInt(0, 4); //nearly random string
        }
        int64_t level = -1; //st_randomInt(-1, 10);
        int64_t compressedSizeInBytes;
        //Do the compression
        time_t startTime = time(NULL);
        void *compressedString = compress(randomString, sizeof(char)*size, &compressedSizeInBytes, level);
        CuAssertTrue(testCase, compressedSizeInBytes >= 0);
        st_logDebug("I did a round of compression: I got %" PRIi64  " bytes to compress and %" PRIi64 " bytes compressed at level %" PRIi64 " in %f seconds\n", size, compressedSizeInBytes, level, (float)difftime(time(NULL), startTime));
        //Now decompress
        int64_t size2;
        char *randomString2 = decompress(compressedString, compressedSizeInBytes, &size2);
        st_logDebug("I did a round of compression/decompression: I got %" PRIi64  " bytes to compress and %" PRIi64 " bytes compressed at level %" PRIi64 " in %f seconds\n", size, compressedSizeInBytes, level, (float)difftime(time(NULL), startTime));
        bytesCompressed += compressedSizeInBytes;
        bytesDecompressed += size;
        CuAssertTrue(testCase, size == size2);
        for(int64_t j=0; j<size; j++) {
            CuAssertIntEquals(testCase, randomString[j], randomString2[j]);
        }
        free(randomString);
        free(randomString2);
        free(compressedString);
    }
    st_logDebug("Got a compression ration of %f from %" PRIi64 " bytes compressed and %" PRIi64 " bytes uncompressed \n", (double)bytesCompressed/bytesDecompressed, bytesCompressed, bytesDecompressed);
}

/*
 * Does a large number of rounds of compressions and decompression of small strings.
 */
static void test_stCompression_compressAndDecompress_Lots(CuTest *testCase) {
    test_stCompression_compressAndDecompressP(testCase, 1000, 50, 1000, stCompression_compress, stCompression_decompress);
}

/*
 * Does a small number of rounds of large strings.
 */
static void test_stCompression_compressAndDecompress_Big(CuTest *testCase) {
    test_stCompression_compressAndDecompressP(testCase, 5, 10000000, 50000000, stCompression_compress, stCompression_decompress);
}

/*
 * Does a large number of rounds of compressions and decompression of small strings.
 */
static void test_stCompression_compressAndDecompress_Lots_Zlib(CuTest *testCase) {
    test_stCompression_compressAndDecompressP(testCase, 1000, 50, 100, stCompression_compressZlib, stCompression_decompressZlib);
}

/*
 * Does a small number of rounds of large strings.
 */
static void test_stCompression_compressAndDecompress_Big_Zlib(CuTest *testCase) {
    test_stCompression_compressAndDecompressP(testCase, 5, 10000000, 50000000, stCompression_compressZlib, stCompression_decompressZlib);
}

CuSuite* sonLib_stCompressionTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stCompression_compressAndDecompress_Lots);
    SUITE_ADD_TEST(suite, test_stCompression_compressAndDecompress_Big);
    SUITE_ADD_TEST(suite, test_stCompression_compressAndDecompress_Lots_Zlib);
        SUITE_ADD_TEST(suite, test_stCompression_compressAndDecompress_Big_Zlib);
    return suite;
}
