/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * common code for database tests.
 */
#include "sonLibGlobalsTest.h"
#include "kvDatabaseTestCommon.h"
#include "stSafeC.h"
#include <getopt.h>

/* usage message and exit */
static void usage(const char *desc) {
    static const char *help = 
        "Options:\n"
        "\n"
        "-t --type=dbtype - one of 'KyotoTycoon', 'TokyoCabinet' or 'MySql'.\n"
        "    Values area case-insensitive, defaults to TokyoCabinet.\n"
        "-d --db=database - database directory for TokyoCabinet or database name\n"
        "    for SQL databases. Defaults to testTCDatabase for TokyoCabinet,\n"
        "    SQL databases must specify.\n"
        "--host=host - Tycoon or SQL database host, defaults to localhost\n"
        "--port=port - Tycoon or SQL database port.\n"
        "-u, --user=user - SQL database user.\n"
        "-p, --pass=pass - SQL database password.\n"
        "-h, --help - print this message.\n";
    fprintf(stderr, "%s\n%s\n", desc, help);
    exit(1);
}

static stKVDatabaseType parseDbType(const char *dbTypeStr) {
    if (stString_eqcase(dbTypeStr, "TokyoCabinet")) {
        return stKVDatabaseTypeTokyoCabinet;
    } else if (stString_eqcase(dbTypeStr, "KyotoTycoon")) {
        return stKVDatabaseTypeKyotoTycoon;
    } else if (stString_eqcase(dbTypeStr, "MySql")) {
        return stKVDatabaseTypeMySql;
    } else if (stString_eqcase(dbTypeStr, "Redis")) {
        return stKVDatabaseTypeRedis;
    } else {
        fprintf(stderr, "Error: invalid value for --type: %s\n", dbTypeStr);
        exit(1);
    }
}

/* Parse options for specifying database to tests.  Fill in positional argument vector, setting unused ones to NULL.
 * The positional and numPositionalRet maybe NULL. */
stKVDatabaseConf *kvDatabaseTestParseOptions(int argc, char *const *argv, const char *desc, int minNumPositional, int maxNumPositional,
                                             char **positional, int *numPositionalRet) {
    static struct option longOptions[] = {
        {"type", required_argument, NULL, 't'},
        {"db", required_argument,   NULL, 'd'},
        {"host", required_argument, NULL, 'H'},
        {"port", required_argument, NULL, 'P'},
        {"timeout", required_argument, NULL, 'i'},
        {"maxKTRecordSize", required_argument, NULL, 'r'},
        {"maxKTBulkSetSize", required_argument, NULL, 'b'},
	{"maxRedisBulkSetSize", required_argument, NULL, 'R'},
        {"user", required_argument, NULL, 'u'},
        {"pass", required_argument, NULL, 'p'},
        {"name", required_argument, NULL, 'n'},
        {"help", no_argument,       NULL, 'h'},
        {NULL, 0, NULL, '\0'}
    };
    stKVDatabaseType optType = stKVDatabaseTypeTokyoCabinet;
    const char *optDb = "testTCDatabase";
    const char *optHost = "localhost";
    unsigned int optPort = 0;
    int optTimeout = -1;
    int64_t optMaxKTRecordSize = (int64_t) 1U << 27;
    int64_t optMaxKTBulkSetSize = (int64_t) 1U << 27;
    int64_t optMaxKTBulkSetNumRecords = (int64_t) 1U << 27;
    int64_t optMaxRedisBulkSetSize = (int64_t) 1U << 30;
    const char *optUser = NULL;
    const char *optPass = NULL;
    const char *optName = NULL;
    int optKey, optIndex;
    while ((optKey = getopt_long(argc, argv, "t:d:H:P:i:r:b:R:u:p:h", longOptions, &optIndex)) >= 0) {
        switch (optKey) {
        case 't':
            optType = parseDbType(optarg);
            break;
        case 'd':
            optDb = optarg;
            break;
        case 'H':
            optHost = optarg;
            break;
        case 'P':
            optPort = stSafeStrToUInt32(optarg);
            break;
        case 'i':
            optTimeout = stSafeStrToUInt32(optarg);
            break;
        case 'r':
        	optMaxKTRecordSize = stSafeStrToInt64(optarg);
			break;
        case 'b':
        	optMaxKTBulkSetSize = stSafeStrToInt64(optarg);
			break;
	case 'R':
                optMaxRedisBulkSetSize = stSafeStrToInt64(optarg);
                        break;
        case 'u':
            optUser = optarg;
            break;
        case 'p':
            optPass = optarg;
            break;
        case 'n':
        	optName = optarg;
        	break;
        case 'h':
            usage(desc);
            break;
        default:
            fprintf(stderr, "Error: invalid  option: %s\n", argv[optind]); 
            usage(desc);
        }
    }
    int numPositional = argc - optind;
    if (!((minNumPositional <= numPositional) && (numPositional <= maxNumPositional))) {
        fprintf(stderr, "Error: wrong number of arguments\n");
        usage(desc);
    }
    if (positional != NULL) {
        int i;
        for (i = 0; i < numPositional; i++) {
            positional[i] = argv[optind+i];
        }
        for (; i < maxNumPositional; i++) {
            positional[i] = NULL;
        }
    }
    if (numPositionalRet != NULL) {
        *numPositionalRet = numPositional;
    }
    stKVDatabaseConf *conf = NULL;
    if (optType == stKVDatabaseTypeTokyoCabinet) {
        conf = stKVDatabaseConf_constructTokyoCabinet(optDb);
        fprintf(stderr, "running Tokyo Cabinet sonLibKVDatabase tests\n");
    } else if (optType == stKVDatabaseTypeKyotoTycoon) {
        conf = stKVDatabaseConf_constructKyotoTycoon(optHost, optPort, optTimeout,
        		optMaxKTRecordSize, optMaxKTBulkSetSize, optMaxKTBulkSetNumRecords, optDb, optName);
        fprintf(stderr, "running Kyoto Tycoon sonLibKVDatabase tests\n");
    } else if (optType == stKVDatabaseTypeMySql) {
        conf = stKVDatabaseConf_constructMySql(optHost, 0, optUser, optPass, optDb, "cactusDbTest");
        fprintf(stderr, "running MySQL sonLibKVDatabase tests\n");
    } else if (optType == stKVDatabaseTypeRedis) {
        conf = stKVDatabaseConf_constructRedis(optHost, optPort, optMaxRedisBulkSetSize);
        fprintf(stderr, "running Redis sonLibKVDatabase tests\n");
    }
    return conf;
}

