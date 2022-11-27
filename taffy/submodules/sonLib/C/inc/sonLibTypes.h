/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibGlobals.h
 *
 *  Created on: 21 May 2010
 *      Author: benedictpaten
 */

#ifndef SONLIB_GLOBALS_H_
#define SONLIB_GLOBALS_H_

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
    
#include <inttypes.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic data structure declarations (contents hidden)
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
    
typedef struct _stTree stTree;
typedef struct _stHash stHash;
typedef struct _stSet stSet;
typedef struct hashtable_itr stHashIterator;
typedef struct _stSetIterator stSetIterator;
typedef struct _stSortedSet stSortedSet;
typedef struct _stSortedSetIterator stSortedSetIterator;
typedef struct _stList stList;
typedef struct _stListIterator stListIterator;
typedef int64_t stIntTuple;
typedef double stDoubleTuple;
typedef struct stExcept stExcept;
typedef struct stCache stCache;
typedef struct stKVDatabase stKVDatabase;
typedef struct stKVDatabaseConf stKVDatabaseConf;
typedef struct stKVDatabaseBulkRequest stKVDatabaseBulkRequest;
typedef struct stKVDatabaseBulkResult stKVDatabaseBulkResult;
typedef struct _stEdge stEdge;
typedef struct _stGraph stGraph;
typedef struct _stPosetAlignment stPosetAlignment;
typedef struct _stTreap stTreap;
typedef struct _stEulerVertex stEulerVertex;
typedef struct _stEulerHalfEdge stEulerHalfEdge;
typedef struct _stEulerTour stEulerTour;
typedef struct _stEulerTourIterator stEulerTourIterator;
typedef struct _stEulerTourEdgeIterator stEulerTourEdgeIterator;
typedef struct _stEulerTourComponentIterator stEulerTourComponentIterator;
typedef struct _stConnectivity stConnectivity;
typedef struct _stConnectedComponent stConnectedComponent;
typedef struct _stConnectedComponentIterator stConnectedComponentIterator;
typedef struct _stConnectedComponentNodeIterator stConnectedComponentNodeIterator;
typedef struct _stEdgeContainer stEdgeContainer;
typedef struct _stEdgeContainerIterator stEdgeContainerIterator;
typedef struct _stIncidentEdgeList stIncidentEdgeList;
typedef struct _stNaiveConnectivity stNaiveConnectivity;
typedef struct _stNaiveConnectedComponent stNaiveConnectedComponent;
typedef struct _stNaiveConnectedComponentIterator stNaiveConnectedComponentIterator;
typedef struct _stNaiveConnectedComponentNodeIterator stNaiveConnectedComponentNodeIterator;
typedef struct _stMatrix stMatrix;

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // SONLIBGLOBALS_H_
