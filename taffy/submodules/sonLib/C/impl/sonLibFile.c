/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibFile.c
 *
 *  Created on: 7 Sep 2010
 *      Author: benedictpaten
 */
#include "sonLibGlobalsInternal.h"
#include <dirent.h>
#include <sys/stat.h>
#include <inttypes.h>

const char *ST_FILE_EXCEPTION = "ST_FILE_EXCEPTION";

static int64_t stFile_getLineFromFileP(char **s, int64_t *n, FILE *f) {
    register int64_t nMinus1 = ((*n) - 1), i = 0;

    char *s2 = *s;
    while (1) {
        register int64_t ch = (char) getc(f);

        if (ch == '\r') {
            ch = getc(f);
        }

        if (i == nMinus1) {
            *n = 2 * (*n) + 1;
            *s = realloc(*s, (*n + 1) * sizeof(char));
            assert(*s != NULL);
            s2 = *s + i;
            nMinus1 = ((*n) - 1);
        }

        if ((ch == '\n') || (ch == EOF)) {
            *s2 = '\0';
            return (feof(f) ? -1 : i);
        } else {
            *s2 = ch;
            s2++;
        }
        ++i;
    }
    return 0;
}

char *stFile_getLineFromFile(FILE *fileHandle) {
    int64_t length = 100;
    char *cA = st_malloc(length * sizeof(char));
    int64_t i = stFile_getLineFromFileP(&cA, &length, fileHandle);
    if (i == -1 && strlen(cA) == 0) {
        free(cA);
        return NULL;
    }
    char *cA2 = stString_copy(cA);
    free(cA);
    return cA2;
}

stList *stFile_getLinesFromFile(char *fileName) {
    FILE *fh = st_fopen(fileName, "r");
    stList *lines = stList_construct3(0, free);
    char *line;
    while ((line = stFile_getLineFromFile(fh)) != NULL) {
        stList_append(lines, line);
    }
    fclose(fh);
    return lines;
}

char *stFile_pathJoin(const char *pathPrefix, const char *pathSuffix) {
    return strlen(pathPrefix) == 0 || pathPrefix[strlen(pathPrefix) - 1] == '/' ? stString_print(
            "%s%s", pathPrefix, pathSuffix)
                                                                                : stString_print("%s/%s", pathPrefix,
                                                                                                 pathSuffix);
}

bool stFile_exists(const char *fileName) {
    struct stat info;
    return !stat(fileName, &info);
}

bool stFile_isDir(const char *fileName) {
    if (!stFile_exists(fileName)) {
        stThrowNew(ST_FILE_EXCEPTION, "File does not exist exist: %s\n",
                fileName);
    }
    struct stat info;
    if (stat(fileName, &info)) {
        stThrowNew(ST_FILE_EXCEPTION,
                "Failed to get information about the file: %s\n", fileName);
    }
    return S_ISDIR(info.st_mode);
}

void stFile_mkdir(const char *dirName) {
    if (stFile_exists(dirName)) {
        stThrowNew(ST_FILE_EXCEPTION,
                "Directory already appears to exist: %s\n", dirName);
    }
    int i = mkdir(dirName, S_IRWXU | S_IRWXG | S_IRWXO);
    if (i) {
        stThrowNew(ST_FILE_EXCEPTION, "Error code from mkdir function: %" PRIi64 "\n",
                i);
    }
}

void stFile_mkdirp(const char *dirName) {
    st_system("mkdir -m 0777 -p %s", dirName);
}


stList *stFile_getFileNamesInDirectory(const char *dir) {
    if (!stFile_exists(dir)) {
        stThrowNew(ST_FILE_EXCEPTION, "Directory does not exist exist: %s\n",
                dir);
    }
    if (!stFile_isDir(dir)) {
        stThrowNew(ST_FILE_EXCEPTION,
                "Specified file is not a directory: %s\n", dir);
    }
    stList *files = stList_construct3(0, free);
    DIR *dh = opendir(dir);
    struct dirent *file;//a 'directory entity' AKA file
    while ((file = readdir(dh)) != NULL) {
        if (file->d_name[0] != '.') {
            stList_append(files, stString_copy(file->d_name));
        }
    }
    closedir(dh);
    return files;
}

void stFile_rmrf(const char *fileName) {
    if (!stFile_exists(fileName)) {
        stThrowNew(ST_FILE_EXCEPTION, "File to remove does not exist: %s\n",
                fileName);
    }
    stFile_rmtree(fileName);
}

void stFile_rmtree(const char *fileName) {
    st_system("rm -rf %s", fileName);
}

FILE *st_fopen(const char *fileName, const char *mode) {
    FILE *ret = fopen(fileName, mode);
    if (ret == NULL) {
        st_errnoAbort("Could not open file %s", fileName);
    }
    return ret;
}
