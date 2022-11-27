/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibFile.h
 *
 *  Created on: 7 Sep 2010
 *      Author: benedictpaten
 */

#ifndef SONLIBFILE_H_
#define SONLIBFILE_H_

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

//The exception string
extern const char *ST_FILE_EXCEPTION;

/*
 * Reads a line from a file (which may be terminated by a newline char or EOF),
 * returning the line excluding the newline character.
 * If the file has hit the EOF then it returns NULL.
 */
char *stFile_getLineFromFile(FILE *fileHandle);

/*
 * Reads the line from the give file, returning lines in a list in order. Newlines/EOF characters are removed from the lines.
 */
stList *stFile_getLinesFromFile(char *fileName);

/*
 * Joins together two strings.
 */
char *stFile_pathJoin(const char *pathPrefix, const char *pathSuffix);

/*
 * Returns non-zero iff the file exists.
 */
bool stFile_exists(const char *fileName);

/*
 * Returns non-zero iff the file is a directory. Raises an exception if the file does not exist.
 */
bool stFile_isDir(const char *fileName);

/*
 * Get list of file names (as strings) in a directory. Raises an exception if dir is not a directory.
 */
stList *stFile_getFileNamesInDirectory(const char *dir);

/*
 * Creates a directory with 777 access permissions, throws exceptions if unsuccessful.
 */
void stFile_mkdir(const char *dirName);

/*
 * Creates a directory with 777 access permissions, throws exceptions if unsuccessful.
 */
void stFile_mkdir(const char *dirName);

/*
 * Creates a directory and missing parents with 777 access permissions, throws exceptions if unsuccessful.
 */
void stFile_mkdirp(const char *dirName);

/*
 * Forceably remove a file. If a dir, removes dir and children. Be careful.
 * The file must exist.
 */
void stFile_rmrf(const char *fileName);

/*
 * Forceably remove a file. If a dir, removes dir and children. Be careful.
 * The file does not have to exist.
 */
void stFile_rmtree(const char *fileName);

/*
 * Thin wrapper around fopen that dies with a perror if the operation
 * fails. Meant to be more informative than a random segfault later
 * down the road.
 */
FILE *st_fopen(const char *fileName, const char *mode);

#ifdef __cplusplus
}
#endif
#endif /* SONLIBFILE_H_ */
