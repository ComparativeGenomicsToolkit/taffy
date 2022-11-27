/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef kvDatabaseTestCommon_h
#define  kvDatabaseTestCommon_h

/* Parse options for specifying database to tests.  Fill in positional argument vector, setting unused ones to NULL.
 * The positional and numPositionalRet maybe NULL. */
stKVDatabaseConf *kvDatabaseTestParseOptions(int argc, char *const *argv, const char *desc, int minNumPositional, int maxNumPositional,
                                             char **positional, int *numPositionalRet);

#endif
