#ifndef EULER_H
#define EULER_H
stEulerVertex *stEulerVertex_construct(void *vertexID);
void stEulerVertex_destruct(stEulerVertex *vertex);
stTreap *stEulerVertex_incidentEdgeA(stEulerVertex *vertex);
stTreap *stEulerVertex_incidentEdgeB(stEulerVertex *vertex);
int stEulerVertex_connected(stEulerVertex *from, stEulerVertex *to);
int stEulerVertex_isSingleton(stEulerVertex *vertex);
//-------------------------------------------------------
stEulerHalfEdge *stEulerHalfEdge_construct(void);
void stEulerHalfEdge_destruct(stEulerHalfEdge *edge);
int stEulerHalfEdge_contains(stEulerHalfEdge *edge, stEulerVertex *vertex);
void *stEulerHalfEdge_getFrom(stEulerHalfEdge *edge);
void *stEulerHalfEdge_getTo(stEulerHalfEdge *edge);
//--------------------------------------------------
stEulerTour *stEulerTour_construct(void);
bool stEulerTour_hasEdge(stEulerTour *et, void *u, void *v);
stEdgeContainer *stEulerTour_getEdges(stEulerTour *et);
bool stEulerTour_isSingleton(stEulerTour *et, void *v);
void *stEulerTour_getConnectedComponent(stEulerTour *et, void *v);
void stEulerTour_printTour(stEulerTour *et, void *v);
void stEulerTour_printEdgeTour(stEulerTour *et, void *v);
int stEulerTour_connected(stEulerTour *et, void *u, void *v);
int stEulerTour_size(stEulerTour *et, void *v);
int stEulerTour_getNComponents(stEulerTour *et);
stTreap *stEulerTour_findRoot(stEulerTour *et, void *v);
void *stEulerTour_findRootNode(stEulerTour *et, void *v);
stEulerVertex *stEulerTour_getVertex(stEulerTour *et, void *v);
stEulerVertex *stEulerTour_createVertex(stEulerTour *et, void *vertexID);
void stEulerTour_removeVertex(stEulerTour *et, void *vertexID);
void stEulerTour_destruct(stEulerTour *et);
void stEulerTour_makeRoot(stEulerTour *et, stEulerVertex *vertex);
void stEulerTour_link(stEulerTour *et, void *u, void *v);
void stEulerTour_cut(stEulerTour *et, void *u, void *v);
//-----------------------------------------------------------
stEulerTourIterator *stEulerTour_getIterator(stEulerTour *et, void *v);
void *stEulerTourIterator_getNext(stEulerTourIterator *it);
void stEulerTourIterator_destruct(stEulerTourIterator *it);
stSet *stEulerTour_getNodesInComponent(stEulerTour *et, void *v);
//-------------------------------------------------------------
stEulerTourEdgeIterator *stEulerTour_getEdgeIterator(stEulerTour *et, void *v);
bool stEulerTourEdgeIterator_getNext(stEulerTourEdgeIterator *it, void **node1, void **node2); 
void stEulerTourEdgeIterator_destruct(stEulerTourEdgeIterator *it);
//------------------------------------------------------------
stEulerTourComponentIterator *stEulerTour_getComponentIterator(stEulerTour *et);
void *stEulerTourComponentIterator_getNext(stEulerTourComponentIterator *it);
void stEulerTourComponentIterator_destruct(stEulerTourComponentIterator *it);
#endif
