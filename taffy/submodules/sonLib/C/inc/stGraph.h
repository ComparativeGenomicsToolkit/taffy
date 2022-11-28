/*
 * stGraph.h
 *
 *  Created on: 15 Mar 2013
 *      Author: benedictpaten
 */

#ifndef STGRAPH_H_
#define STGRAPH_H_

#include "sonLibTypes.h"

/*
 * Make graph with vertices 0, 1 ... vertexNo.
 */
stGraph *stGraph_construct(int64_t vertexNo);

void stGraph_destruct(stGraph *g);

/*
 * Number of vertices in graph.
 */
int64_t stGraph_cardinality(stGraph *g);

/*
 * Get edges incident with v, returned as linked list.
 */
stEdge *stGraph_getEdges(stGraph *graph, int64_t v);

/*
 * Gets next edge incident with vertex, or NULL if no more.
 */
stEdge *stEdge_nextEdge(stEdge *e);

double stEdge_weight(stEdge *e);

int64_t stEdge_to(stEdge *e);

/*
 * Adds an edge to graph.
 */
void stGraph_addEdge(stGraph *graph, int64_t v1, int64_t v2, double weight);

/*
 * Computes dijkstras, returning shortest path distances between
 * chosen vertex and other vertices, as an array of doubles.
 */
double *stGraph_shortestPaths(stGraph *g, int64_t sourceVertex);

#endif /* STGRAPH_H_ */
