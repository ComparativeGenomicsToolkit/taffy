/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibCompression.h
 *
 *  Created on: 03-Sep-2010
 *      Author: benedictpaten
 */

#ifndef SONLIBCOMPRESSION_H_
#define SONLIBCOMPRESSION_H_

#include "sonLibTypes.h"
#ifdef __cplusplus
extern "C" {
#endif

//The exception string
extern const char *ST_COMPRESSION_EXCEPTION_ID;

/*
 * Compresses the data and returns it. sizeInBytes in the size of the uncompressed data array, the pointer
 * compressedSizeInBytes is given the size of the compressed string. The level is a value between 0 and 10 giving the
 * degree of required compression. If -1 is given then the default level is used.
 */
void *stCompression_compress(void *data, int64_t sizeInBytes, int64_t *compressedSizeInBytes, int64_t level);

/*
 * Decompresses the compressed data string of size compressedSizeInBytes, initialises the sizeInBytes point to the size
 * of the decompressed string.
 */
void *stCompression_decompress(void *compressedData, int64_t compressedSizeInBytes, int64_t *sizeInBytes);

/*
 * Uses Zlib.
 */
void *stCompression_compressZlib(void *data, int64_t sizeInBytes, int64_t *compressedSizeInBytes, int64_t level);

/*
 * Uses Zlib.
 */
void *stCompression_decompressZlib(void *compressedData, int64_t compressedSizeInBytes, int64_t *sizeInBytes);

#ifdef __cplusplus
}
#endif
#endif /* SONLIBCOMPRESSION_H_ */
