/*
 * A data structure for very efficiently keeping track of connected
 * components in an undirected multigraph.
 */
#ifndef SONLIB_CONNECTIVITY_H_
#define SONLIB_CONNECTIVITY_H_

/*
 * Create a new graph.
 */
stConnectivity *stConnectivity_construct(void);

/*
 * Free a graph properly.
 */
void stConnectivity_destruct(stConnectivity *connectivity);

/*
 * Check whether two nodes are connected.
 */
int stConnectivity_connected(stConnectivity *connectivity, void *node1, void *node2);

/*
 * Get the edges in the graph.
 */
stEdgeContainer *stConnectivity_getEdges(stConnectivity *connectivity);

/*
 * Add a new isolated node to the graph.
 */
void stConnectivity_addNode(stConnectivity *connectivity, void *node);

/*
 * Add an edge to the graph. The two nodes must have been added
 * beforehand using addNode.
 */
void stConnectivity_addEdge(stConnectivity *connectivity, void *node1, void *node2);

/*
 * Check whether the graph has at least one edge between node1 and node2.
 */
bool stConnectivity_hasEdge(stConnectivity *connectivity, void *node1, void *node2);

/*
 * Remove a single edge between node1 and node2.
 */
void stConnectivity_removeEdge(stConnectivity *connectivity, void *node1, void *node2);

/*
 * Remove a node, and all its edges, from the graph.
 */
void stConnectivity_removeNode(stConnectivity *connectivity, void *node);

/*
 * Callbacks, to be notified when certain events happen to connected components.
 */

/*
 * A new component was created.
 */
void stConnectivity_setCreationCallback(stConnectivity *connectivity, void (*callback)(void *, stConnectedComponent *), void *extraData);

/*
 * The first component was merged into the second component.
 */
void stConnectivity_setMergeCallback(stConnectivity *connectivity, void (*callback)(void *, stConnectedComponent *, stConnectedComponent *), void *extraData);

/*
 * The ends in the set were cleaved from the first component into the second, new component.
 */
void stConnectivity_setCleaveCallback(stConnectivity *connectivity, void (*callback)(void *, stConnectedComponent *, stConnectedComponent *, stSet *), void *extraData);

/*
 * A component was deleted.
 */
void stConnectivity_setDeletionCallback(stConnectivity *connectivity, void (*callback)(void *, stConnectedComponent *), void *extraData);

stConnectedComponent *stConnectedComponent_construct(stConnectivity *connectivity, void *node);

stConnectedComponent *stConnectivity_getConnectedComponent(stConnectivity *connectivity, void *node);

int stConnectivity_getNComponents(stConnectivity *connectivity);

void *stConnectedComponent_getNodeInComponent(stConnectedComponent *component);

stConnectedComponentNodeIterator *stConnectedComponent_getNodeIterator(stConnectedComponent *component);

void *stConnectedComponentNodeIterator_getNext(stConnectedComponentNodeIterator *it);

void stConnectedComponentNodeIterator_destruct(stConnectedComponentNodeIterator *it);

stConnectedComponentIterator *stConnectivity_getConnectedComponentIterator(stConnectivity *connectivity);

stConnectedComponent *stConnectedComponentIterator_getNext(stConnectedComponentIterator *it);

void stConnectedComponentIterator_destruct(stConnectedComponentIterator *it);
void stConnectedComponent_destruct(stConnectedComponent *comp);
stEdgeContainer *stEdgeContainer_construct(void(*destructEdge)(void *));
void stEdgeContainer_destruct(stEdgeContainer *container);
void *stEdgeContainer_getEdge(stEdgeContainer *container, void *u, void *v); 
void stEdgeContainer_deleteEdge(stEdgeContainer *container, void *u, void *v); 
bool stEdgeContainer_hasEdge(stEdgeContainer *container, void *u, void *v); 
void stEdgeContainer_addEdge(stEdgeContainer *container, void *u, void *v, void *edge); 
stList *stEdgeContainer_getIncidentEdgeList(stEdgeContainer *container, void *v);
stEdgeContainerIterator *stEdgeContainer_getIterator(stEdgeContainer *container); 
bool stEdgeContainer_getNext(stEdgeContainerIterator *it, void **node1, void **node2); 
void stEdgeContainer_destructIterator(stEdgeContainerIterator *it); 
#endif
