#ifndef NAIVE_CONNECTIVITY_H
#define NAIVE_CONNECTIVITY_H
stNaiveConnectivity *stNaiveConnectivity_construct(void);
stSet *stNaiveConnectedComponent_getNodes(stNaiveConnectedComponent *comp);
void stNaiveConnectivity_destruct(stNaiveConnectivity *connectivity);
void stNaiveConnectivity_addNode(stNaiveConnectivity *connectivity, void *node);
void stNaiveConnectivity_addEdge(stNaiveConnectivity *connectivity, void *node1, void *node2);
bool stNaiveConnectivity_hasEdge(stNaiveConnectivity *connectivity, void *node1, void *node2);
void stNaiveConnectivity_removeEdge(stNaiveConnectivity *connectivity, void *node1, void *node2);
void stNaiveConnectivity_removeNode(stNaiveConnectivity *connectivity, void *node);
stNaiveConnectedComponent *stNaiveConnectivity_getConnectedComponent(stNaiveConnectivity 
		*connectivity, void *node);
stNaiveConnectedComponentNodeIterator *stNaiveConnectedComponent_getNodeIterator(stNaiveConnectedComponent 
		*component);
void *stNaiveConnectedComponentNodeIterator_getNext(stNaiveConnectedComponentNodeIterator *it);
void stNaiveConnectedComponentNodeIterator_destruct(stNaiveConnectedComponentNodeIterator *it);
stNaiveConnectedComponentIterator *stNaiveConnectivity_getConnectedComponentIterator(stNaiveConnectivity 
		*connectivity);
stNaiveConnectedComponent *stNaiveConnectedComponentIterator_getNext(stNaiveConnectedComponentIterator *it);
void stNaiveConnectedComponentIterator_destruct(stNaiveConnectedComponentIterator *it); 
#endif
