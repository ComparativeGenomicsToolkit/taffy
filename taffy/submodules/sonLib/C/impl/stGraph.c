#include "sonLibGlobalsInternal.h"

struct _stEdge {
    int64_t to;
    double weight;
    stEdge *nEdge;
};

struct _stGraph {
    int64_t vertexNo;
    stEdge **adjLists;
};

stGraph *stGraph_construct(int64_t vertexNo) {
    stGraph *graph = st_malloc(sizeof(stGraph));
    graph->vertexNo = vertexNo;
    graph->adjLists = st_calloc(vertexNo, sizeof(stEdge *));
    return graph;
}

void stGraph_destruct(stGraph *g) {
    for(int64_t v=0; v<g->vertexNo; v++) {
        stEdge *e = g->adjLists[v];
        while(e != NULL) {
            stEdge *e2 = e;
            e = e->nEdge;
            free(e2);
        }
    }
    free(g->adjLists);
    free(g);
}

int64_t stGraph_cardinality(stGraph *g) {
    return g->vertexNo;
}

stEdge *stGraph_getEdges(stGraph *graph, int64_t v) {
    return graph->adjLists[v];
}

stEdge *stEdge_nextEdge(stEdge *e) {
    return e->nEdge;
}

double stEdge_weight(stEdge *e) {
    return e->weight;
}

int64_t stEdge_to(stEdge *e) {
    return e->to;
}

static void stGraph_addEdgeP(stGraph *graph, int64_t v1, int64_t v2, double weight) {
    stEdge *e = st_malloc(sizeof(stEdge));
    e->to = v2;
    e->weight = weight;
    stEdge *e2 = graph->adjLists[v1];
    graph->adjLists[v1] = e;
    e->nEdge = e2;
}

void stGraph_addEdge(stGraph *graph, int64_t v1, int64_t v2, double weight) {
    stGraph_addEdgeP(graph, v1, v2, weight);
    stGraph_addEdgeP(graph, v2, v1, weight);
}

typedef struct _VDistance {
    int64_t v;
    double distance;
} VDistance;

int vDistance_cmp(VDistance *vD1, VDistance *vD2) {
    return vD1->distance > vD2->distance ? 1 : (vD1->distance < vD2->distance ? -1 : (vD1->v > vD2->v ? 1 : (vD1->v < vD2->v ? -1 : 0)));
}

double *stGraph_shortestPaths(stGraph *g, int64_t sourceVertex) {
    stSortedSet *orderedDistances = stSortedSet_construct3((int (*)(const void *, const void *))vDistance_cmp, NULL);
    VDistance *distances = st_malloc(sizeof(VDistance) * stGraph_cardinality(g));
    for(int64_t v=0; v<stGraph_cardinality(g); v++) {
        VDistance *vD = &distances[v];
        vD->v = v;
        vD->distance = (v == sourceVertex ? 0 : INT64_MAX);
        stSortedSet_insert(orderedDistances, vD);
    }
    while(stSortedSet_size(orderedDistances) > 0) {
        VDistance *vD = stSortedSet_getFirst(orderedDistances);
        stSortedSet_remove(orderedDistances, vD);
        stEdge *e = stGraph_getEdges(g, vD->v);
        while(e != NULL) {
            double d = vD->distance + e->weight;
            VDistance *vD2 = &distances[e->to];
            if(vD2->distance > d) {
                assert(stSortedSet_search(orderedDistances, vD2) != NULL);
                stSortedSet_remove(orderedDistances, vD2);
                vD2->distance = d;
                stSortedSet_insert(orderedDistances, vD2);
            }
            e = stEdge_nextEdge(e);
        }
    }
    double *dA = st_malloc(sizeof(double) * stGraph_cardinality(g));
    for(int64_t v=0; v<stGraph_cardinality(g); v++) {
        dA[v] = distances[v].distance;
    }
    stSortedSet_destruct(orderedDistances);
    free(distances);
    return dA;
}
