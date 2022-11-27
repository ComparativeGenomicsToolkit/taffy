// Dynamically keeps track of the connected components in an
// undirected graph. Adapted from the implementation at www.github.com/sambayless/monosat. The algorithm was 
// described by Mikkel Thorup in "Near-optimal Fully Dynamic Graph Connectivity", 2001. 
//
// The graph has log(n_nodes) levels, and each edge in the graph is assigned a level, starting at level
// (n_levels - 1) when the edge is added. On each level i, a spanning forest is maintained using edges of 
// level at least i. The spanning forest is represented by an Euler Tour tree of the graph. The Euler Tour
// tree represents a tour of the graph by a balanced binary search tree of the edges. 
#include "sonLibGlobalsInternal.h"
#include <assert.h>


//Public data structures

struct _stConnectivity {
	// Data structure keeping track of the nodes, edges, and connected
	// components.
	stSet *nodes;
	int nNodes;
	int nEdges;
	stList *et; //list of the Euler Tour Trees for each level in the graph
	int nLevels;
	stHash *connectedComponents; // Keyed by Euler tour root
	
	//stores every edge in the graph. Each edge object stores the level of the edge and
	//whether or not it's in the spanning forest. For nodes a and b, the edge object is stored
	//in either (a, b) or (b, a), but not both. Both possibilities are checked when getting the
	//edge.
	stEdgeContainer *edges;

	void (*creationCallback)(void *, stConnectedComponent *);
	void *creationExtraData;
	void (*mergeCallback)(void *, stConnectedComponent *, stConnectedComponent *);
	void *mergeExtraData;
	void (*cleaveCallback)(void *, stConnectedComponent *, stConnectedComponent *, stSet *);
	void *cleaveExtraData;
	void (*deletionCallback)(void *, stConnectedComponent *);
	void *deletionExtraData;

};

struct _stConnectedComponent {
	// Data structure representing a connected component. 
	void *nodeInComponent; //the root node of the tour
	stConnectivity *connectivity;

};

struct _stConnectedComponentIterator {
	// Iterator data structure for components in the graph.
	stList *componentList;
	stEulerTourComponentIterator *compIt;
	stConnectivity *connectivity;
};

struct _stConnectedComponentNodeIterator {
	// Iterator for nodes in a component
	
	//set of the nodes that have already been visited in the Euler Tour. The Euler Tour visits each node
	//at least twice, but the iterator should only return each node once.
	stSet *seen;

	//Euler Tour iterator for this connected component on the top level. 
	stEulerTourIterator *tourIterator;
};
//Private data structures-----------------------------------
struct DynamicEdge {
	void *edgeID;
	void *from;
	void *to;

	//whether or not this edge is part of the spanning forest on the top level. 
	bool in_forest; 	

	//lower level edges are checked first when searching for a replacement edge.
	//Changing the levels of edges appropriately saves time when removing future edges.
	int level;

	// The number of copies of this edge that exist in the graph.
	int multiplicity;
};
//private methods-------------------------------------------
struct DynamicEdge *DynamicEdge_construct() {
	struct DynamicEdge *edge = st_malloc(sizeof(struct DynamicEdge));
	edge->from = NULL;
	edge->to = NULL;
	edge->in_forest = false;
	edge->level = 0;
	edge->multiplicity = 1;
	return(edge);
}
void DynamicEdge_increment(struct DynamicEdge *edge) {
	edge->multiplicity++;
}
void DynamicEdge_decrement(struct DynamicEdge *edge) {
	edge->multiplicity--;
}
int DynamicEdge_multiplicity(struct DynamicEdge *edge) {
	return edge->multiplicity;
}
void DynamicEdge_destruct(struct DynamicEdge *edge) {
	free(edge);
}
stEulerTour *getTopLevel(stConnectivity *connectivity) {
	return stList_get(connectivity->et, 0);
}

//add new levels to compensate for the addition of new nodes,
//and copy nodes from level 0 to the new levels so that all
//levels have the same nodes
void addLevel(stConnectivity *connectivity) {
	
	stEulerTour *newLevel = stEulerTour_construct();
	stSetIterator *it = stSet_getIterator(connectivity->nodes);
	void *node;
	while((node = stSet_getNext(it))) {
		stEulerTour_createVertex(newLevel, node);
	}
	stList_append(connectivity->et, newLevel);
	stSet_destructIterator(it);
}
void removeLevel(stConnectivity *connectivity) {
	stEulerTour *top = stList_pop(connectivity->et);
	stEulerTour_destruct(top);
}
void resizeIncidentEdgeList(stList *incident, int newsize) {
	for(int i = newsize; i < stList_length(incident); i++) {
		stList_remove(incident, newsize);
	}
}
int getNLevels(int nNodes) {
	return (int) ((int)(log(nNodes)/log(2)) + 1);
}

// Exported methods ----------------------------------------------------------------------

stConnectivity *stConnectivity_construct(void) {
	// Initialize the data structure with an empty graph with 0 nodes.
	stConnectivity *connectivity = st_malloc(sizeof(stConnectivity));
	connectivity->nodes = stSet_construct();

	connectivity->et = stList_construct3(0, (void(*)(void*))stEulerTour_destruct);

	//add level zero Euler Tour
	stList_append(connectivity->et, stEulerTour_construct());
	connectivity->connectedComponents = stHash_construct2(NULL,
			(void(*)(void*))stConnectedComponent_destruct);


	connectivity->nLevels = 0;
	connectivity->nNodes = 0;
	connectivity->nEdges = 0;
	connectivity->edges = stEdgeContainer_construct((void(*)(void*))DynamicEdge_destruct);

	connectivity->creationCallback = NULL;
	connectivity->mergeCallback = NULL;
	connectivity->cleaveCallback = NULL;
	connectivity->deletionCallback = NULL;

	return(connectivity);

}

void stConnectivity_destruct(stConnectivity *connectivity) {
	// Free the memory for the data structure.
	stList_destruct(connectivity->et);
	stSet_destruct(connectivity->nodes);
	stEdgeContainer_destruct(connectivity->edges);
	stHash_destruct(connectivity->connectedComponents);

	free(connectivity);
}



void stConnectivity_addNode(stConnectivity *connectivity, void *node) {
	// Add an isolated node to the graph. This should end up in a new
	// connected component with only one member.
	connectivity->nNodes++;

	//maintain log(n) levels in the graph
	int nLevels = getNLevels(connectivity->nNodes);
	if (nLevels > connectivity->nLevels) {
		// We never shrink the number of levels in the graph,
		// to avoid pathological cases where an entire level
		// needs to be added/deleted several times in succession.
		connectivity->nLevels = nLevels;
	}
	while(stList_length(connectivity->et) < connectivity->nLevels) {
		addLevel(connectivity);
	}


	stSet_insert(connectivity->nodes, node);


	for(int i = 0; i < connectivity->nLevels; i++) {
		//Add a new disconnected node to the Euler tour on level i
		stEulerTour *et_i = stList_get(connectivity->et, i);
		stEulerTour_createVertex(et_i, node);
	}

	if (connectivity->creationCallback) {
		connectivity->creationCallback(connectivity->creationExtraData, stConnectivity_getConnectedComponent(connectivity, node));
	}
}

bool stConnectivity_hasEdge(stConnectivity *connectivity, void *node1, void *node2) {
	return stEdgeContainer_hasEdge(connectivity->edges, node1, node2);
}

void stConnectivity_addEdge(stConnectivity *connectivity, void *node1, void *node2) {
	assert(node1 != node2);
	struct DynamicEdge *edge = stEdgeContainer_getEdge(connectivity->edges, node1, node2);
	if (edge != NULL) {
		// This edge is already present in the graph--just increment the multiplicity.
		DynamicEdge_increment(edge);
		return;
	}
	connectivity->nEdges++;
	struct DynamicEdge *newEdge = DynamicEdge_construct();

	newEdge->from = node1;
	newEdge->to = node2;

	newEdge->level = 0;
	
	//add the edge object as an incident edge to an arbitrary one of the nodes
	stEdgeContainer_addEdge(connectivity->edges, node1, node2, newEdge);


	stEulerTour *et_lowest = getTopLevel(connectivity);
	if(!stEulerTour_connected(et_lowest, node1, node2)) {
		//the two nodes are not already connected, so the new node will be pat of the spanning forest.
		//find the two connected components and invalidate them.
		stConnectedComponent *component1 = stHash_remove(connectivity->connectedComponents, stEulerTour_getConnectedComponent(et_lowest, node1));
		stConnectedComponent *component2 = stHash_remove(connectivity->connectedComponents, stEulerTour_getConnectedComponent(et_lowest, node2));
		stConnectedComponent_destruct(component1);
		//link the level N - 1 Euler Tours together, which corresponds to 
		//adding a tree edge on level N - 1.
		stEulerTour_link(et_lowest, node1, node2);
		newEdge->in_forest = true;

		if (component2) {
			component2->nodeInComponent = stEulerTour_getConnectedComponent(et_lowest, node1);
			stHash_insert(connectivity->connectedComponents, component2->nodeInComponent, component2);
			if (connectivity->mergeCallback) {
				connectivity->mergeCallback(connectivity->mergeExtraData, component1, component2);
			}
		}
	}
	else {
		//the nodes were already connected, so the Euler tour doesn't need 
		//to be updated, and the edge is not part of the spanning forest.
		newEdge->in_forest = false;
	}
	return;

}

int stConnectivity_getNComponents(stConnectivity *connectivity) {
    stEulerTour *topLevel = getTopLevel(connectivity);
    return stEulerTour_getNComponents(topLevel);
}

int stConnectivity_connected(stConnectivity *connectivity, void *node1, void *node2) {
	//check whether node1 and node2 have the same root in the spanning forest
	//on level N - 1.
	stEulerTour *et_lowest = getTopLevel(connectivity);
	return(stEulerTour_connected(et_lowest, node1, node2));
}
stEdgeContainer *stConnectivity_getEdges(stConnectivity *connectivity) {
	return connectivity->edges;
}

void print_list(stList *list) {
	stListIterator *it = stList_getIterator(list);
	void *item;
	while((item = stList_getNext(it))) {
		printf("%p ", item);
	}
	stList_destructIterator(it);
	printf("\n");
}


//Attempts to find a replacement edge connecting two nodes, after a tree edge has been deleted.
//Searches the edges incident to node w for such a replacement edge.
struct DynamicEdge *visit(stConnectivity *connectivity, void *w, void *otherTreeVertex,
		struct DynamicEdge *removedEdge, int level, stSet *seen) {

	if(stSet_search(seen, w)) {
		return(NULL);
	}
	stSet_insert(seen, w);

	//get a list of all edges incident to node w, which are possible replacement tree edges.
	stList *w_incident = stEdgeContainer_getIncidentEdgeList(connectivity->edges, w);
	int w_incident_length = stList_length(w_incident);

	int k, j = 0;

	stEulerTour *et_level = stList_get(connectivity->et, level);
	for(k = 0; k < w_incident_length; k++) {
		void *e_wk_node2 = stList_get(w_incident, k);
		struct DynamicEdge *e_wk = stEdgeContainer_getEdge(connectivity->edges, 
				w, e_wk_node2);
		if (e_wk == removedEdge || e_wk->in_forest) {
			//This edge won't work because it's either already in the forest,
			//or the edge that we're trying to replace.
			continue;
		}
		if (e_wk->level == level) {

			//Get e_wk's other node (not w)
			void *otherNode = e_wk->from == w ? e_wk->to: e_wk->from;
			if (stEulerTour_connected(et_level, otherTreeVertex, otherNode)) {
				//e_wk connects the two components, so it is the desired replacement
				//edge. 
				for(k = k+1; k < w_incident_length; k++) {
					struct DynamicEdge *e_toRemove_node2 = stList_get(w_incident, k);
					struct DynamicEdge *e_toRemove = stEdgeContainer_getEdge(connectivity->edges, 
							w, e_toRemove_node2);
					if(e_toRemove == removedEdge || e_toRemove->in_forest) {
						continue;
					}
					stList_set(w_incident, j++, e_toRemove_node2);

				}
				stList_destruct(w_incident);
				return(e_wk);
			}
			else {
				//the level of this edge can be increased since it does
				//not connect to the other component. This leads to later time 
				//savings because this edge will be tried sooner when searching
				//for a lower level replacement edge. Decreasing the level of
				//this edge pays for considering it as a replacement.
				e_wk->level++;
			}
		}
		stList_set(w_incident, j++, e_wk_node2);

	}
	stList_destruct(w_incident);
	return(NULL);

}

//Remove an edge from the graph and update the connected components. If the edge was
//a tree edge, attempt to find a replacement edge to keep node1 and node2 connected. The
//replacement edge should have the highest possible level.
void stConnectivity_removeEdge(stConnectivity *connectivity, void *node1, void *node2) {
	struct DynamicEdge *edge = stEdgeContainer_getEdge(connectivity->edges, 
			node1, node2);
	assert(edge);
	DynamicEdge_decrement(edge);
	if(DynamicEdge_multiplicity(edge) > 0) {
		// There's still a copy of this edge in the multigraph.
		return;
	}
	if(!edge->in_forest) {
		stEdgeContainer_deleteEdge(connectivity->edges, node1, node2);

		return;
	}

	stConnectedComponent *previousComponent = stHash_search(connectivity->connectedComponents,
                                                                stEulerTour_getConnectedComponent(getTopLevel(connectivity), node1));

	assert(edge->level < connectivity->nLevels - 1);
	for (int i = edge->level + 1; i < connectivity->nLevels; i++) {
		assert(!stEulerTour_connected(stList_get(connectivity->et, i), node1, node2));
	}
	//stEulerTour *et_top = stConnectivity_getTopLevel(connectivity);
	//assert(stEulerTour_connected(et_top, node1, node2));

	struct DynamicEdge *replacementEdge = NULL;
	for (int i = edge->level; !replacementEdge && i >= 0; i--) {
		stSet *seen = stSet_construct();
		stEulerTour *et_i = stList_get(connectivity->et, i);
		assert(stEulerTour_connected(et_i, node1, node2));

		assert(stEulerTour_hasEdge(et_i, node1, node2) == stEulerTour_hasEdge(et_i, node2, node1));
		assert(stEulerTour_hasEdge(et_i, node1, node2));
		stEulerTour_cut(et_i, node1, node2);

		//set node1 equal to id of the vertex in the smaller of the two components that have just
		//been created by deleting the edge
		if(stEulerTour_size(et_i, node2) > stEulerTour_size(et_i, node1)) {
			void *temp = node1;
			node1 = node2; 
			node2 = temp;
		}
		//void *smallerComponent = stEulerTour_size(et_i, node2) > stEulerTour_size(et_i, node1) ? node1 : node2;
		edge->in_forest = false;
		replacementEdge = visit(connectivity, node2, node1, edge, i, seen);

		//go through each edge in the tour on level i
		stEulerTourEdgeIterator *edgeIt = stEulerTour_getEdgeIterator(et_i, node2);
		void *from = NULL;
		void *to = NULL;
		while(stEulerTourEdgeIterator_getNext(edgeIt, &from, &to)) {
			struct DynamicEdge *treeEdge = stEdgeContainer_getEdge(connectivity->edges, from, to);
			assert(treeEdge->in_forest);
			if(treeEdge->level == i) {
				treeEdge->level++;
				assert(treeEdge->level <= connectivity->nLevels - 1);
				stEulerTour *et_te = stList_get(connectivity->et, treeEdge->level);
				assert(!stEulerTour_connected(et_te, treeEdge->from, treeEdge->to));
				stEulerTour_link(et_te, treeEdge->from, treeEdge->to);
				assert(stEulerTour_connected(et_te, treeEdge->from, treeEdge->to));
			}
			for (int n = 0; !replacementEdge && n < 2; n++) {
				void *w = n ? treeEdge->to : treeEdge->from;
				replacementEdge = visit(connectivity, w, node1, edge, i, seen);
			}
		}
		stEulerTourEdgeIterator_destruct(edgeIt);

		if(replacementEdge) {
			assert(replacementEdge != edge);
			assert(replacementEdge->level == i);
			replacementEdge->in_forest = true;

			stEulerTour_link(et_i, replacementEdge->from, replacementEdge->to);
			for(int h = replacementEdge->level - 1; h >= 0; h--) {
				stEulerTour *et_h = stList_get(connectivity->et, h);
				assert(stEulerTour_hasEdge(et_h, node1, node2));
				stEulerTour_cut(et_h, node1, node2);
				assert(!stEulerTour_connected(et_h, node1, node2));
				stEulerTour_link(et_h, replacementEdge->from, replacementEdge->to);
			}
		}
		stSet_destruct(seen);
	}
	if (previousComponent) {
       		//update the component to use the new root
		stHash_remove(connectivity->connectedComponents, previousComponent->nodeInComponent);
	       	previousComponent->nodeInComponent = stEulerTour_getConnectedComponent(getTopLevel(connectivity), node1);
       		stHash_insert(connectivity->connectedComponents, previousComponent->nodeInComponent, previousComponent);
		if (!replacementEdge && connectivity->cleaveCallback) {
			stConnectedComponent *newComponent = stConnectivity_getConnectedComponent(connectivity, node2);
			stSet *newVertices = stSet_construct();
			stConnectedComponentNodeIterator *it = stConnectedComponent_getNodeIterator(newComponent);
			void *vertex;
			while ((vertex = stConnectedComponentNodeIterator_getNext(it)) != NULL) {
				stSet_insert(newVertices, vertex);
			}
			connectivity->cleaveCallback(connectivity->cleaveExtraData, previousComponent, newComponent, newVertices);
			stSet_destruct(newVertices);
			stConnectedComponentNodeIterator_destruct(it);
		}
	}
	stEdgeContainer_deleteEdge(connectivity->edges, node1, node2);
	return;

}

// Might be cool to be able to add or remove several edges at once, if there is a way to
// make that more efficient than adding them one at a time. Depends on the details of 
// the algorithm you use.

void stConnectivity_removeNode(stConnectivity *connectivity, void *node) {
	// Remove a node (and all its edges) from the graph.
	stList *nodeIncident = stEdgeContainer_getIncidentEdgeList(connectivity->edges, node);
	stListIterator *it = stList_getIterator(nodeIncident);
	void *node2;
	while((node2 = stList_getNext(it))) {
		// Work around the fact that edges can have multiple copies.
		while (stConnectivity_hasEdge(connectivity, node, node2)) {
			stConnectivity_removeEdge(connectivity, node, node2);
		}
	}
	stList_destructIterator(it);
	stList_destruct(nodeIncident);
	stSet_remove(connectivity->nodes, node);
	
	connectivity->nNodes--;
	//delete a level if necessary to preserve log(n) levels
	// Temporarily disabled since with the removal of enough
	// nodes, there can be cases where some edges still have a
	// higher level than the graph has.

	// connectivity->nLevels = getNLevels(connectivity->nNodes);
	// if(connectivity->nLevels < stList_length(connectivity->et)) {
	//	removeLevel(connectivity);
	// }

	if (connectivity->deletionCallback) {
		connectivity->deletionCallback(connectivity->deletionExtraData, stConnectivity_getConnectedComponent(connectivity, node));
	}

        stConnectedComponent_destruct(stHash_remove(connectivity->connectedComponents, node));

	for(int i = 0; i < connectivity->nLevels; i++) {
		stEulerTour *et_i = stList_get(connectivity->et, i);
		stEulerTour_removeVertex(et_i, node);
	}

}

void stConnectivity_setCreationCallback(stConnectivity *connectivity, void (*callback)(void *, stConnectedComponent *), void *extraData) {
    connectivity->creationCallback = callback;
    connectivity->creationExtraData = extraData;
}

void stConnectivity_setMergeCallback(stConnectivity *connectivity, void (*callback)(void *, stConnectedComponent *, stConnectedComponent *), void *extraData) {
    connectivity->mergeCallback = callback;
    connectivity->mergeExtraData = extraData;
}

void stConnectivity_setCleaveCallback(stConnectivity *connectivity, void (*callback)(void *, stConnectedComponent *, stConnectedComponent *, stSet *), void *extraData) {
    connectivity->cleaveCallback = callback;
    connectivity->cleaveExtraData = extraData;
}

void stConnectivity_setDeletionCallback(stConnectivity *connectivity, void (*callback)(void *, stConnectedComponent *), void *extraData) {
    connectivity->deletionCallback = callback;
    connectivity->deletionExtraData = extraData;
}

stConnectedComponent *stConnectedComponent_construct(stConnectivity *connectivity, void *node) {
	stConnectedComponent *component = st_malloc(sizeof(stConnectedComponent));
	component->nodeInComponent = node;
	component->connectivity = connectivity;
	return(component);
}

stConnectedComponent *stConnectivity_getConnectedComponent(stConnectivity *connectivity, 
		void *node) {
	// Get the connected component that this node is a member of. If
	// the component is modified, this pointer can be invalidated
	stEulerTour *et_0 = getTopLevel(connectivity);
	void *compNode = stEulerTour_getConnectedComponent(et_0, node);
	stConnectedComponent *comp;
	if ((comp = stHash_search(connectivity->connectedComponents, compNode)) == NULL) {
		comp = stConnectedComponent_construct(connectivity, compNode);
		stHash_insert(connectivity->connectedComponents, compNode, comp);
	}
	return comp;
}
void *stConnectedComponent_getNodeInComponent(stConnectedComponent *component) {
	return(component->nodeInComponent);
}

stConnectedComponentNodeIterator *stConnectedComponent_getNodeIterator(stConnectedComponent 
		*component) {
	// Get an iterator over the nodes in a particular connected
	// component. You can safely assume that the graph won't be
	// modified while this iterator is active.
	stConnectedComponentNodeIterator *it = st_malloc(sizeof(stConnectedComponentNodeIterator));
	stEulerTour *et = getTopLevel(component->connectivity);
	it->seen = stSet_construct();
	it->tourIterator = stEulerTour_getIterator(et, component->nodeInComponent);
	return(it);
}

void *stConnectedComponentNodeIterator_getNext(stConnectedComponentNodeIterator *it) {
	// Return the next node of the connected component, or NULL if all have been traversed.
	//void *currentNode = it->currentNode;
	void *currentNode = stEulerTourIterator_getNext(it->tourIterator);

	if(stSet_search(it->seen, currentNode)) {
		return(stConnectedComponentNodeIterator_getNext(it));
	}
	stSet_insert(it->seen, currentNode);
	return(currentNode);

}

void stConnectedComponentNodeIterator_destruct(stConnectedComponentNodeIterator *it) {
	// Free the iterator data structure.
	stSet_destruct(it->seen);
	stEulerTourIterator_destruct(it->tourIterator);
	free(it);
}

stConnectedComponentIterator *stConnectivity_getConnectedComponentIterator(stConnectivity *connectivity) {
	// Get an iterator over the connected components in the
	// graph. Again, if the graph is modified while the iterator is
	// active, it's ok for it to break.
	stConnectedComponentIterator *it = st_malloc(sizeof(stConnectedComponentIterator));
	it->componentList = stList_construct3(0, (void(*)(void*))stConnectedComponent_destruct);

	stEulerTour *et_0 = getTopLevel(connectivity);
	it->compIt = stEulerTour_getComponentIterator(et_0);
	it->connectivity = connectivity;
	return(it);
}

stConnectedComponent *stConnectedComponentIterator_getNext(stConnectedComponentIterator *it) {
	// Return the next connected component in the graph, or NULL if all have been traversed.
	void *nextNode = stEulerTourComponentIterator_getNext(it->compIt);
	if(!nextNode) return NULL;
	stConnectedComponent *next = stConnectedComponent_construct(it->connectivity, nextNode);
	stList_append(it->componentList, next);
	return(next);
}

void stConnectedComponentIterator_destruct(stConnectedComponentIterator *it) {
	// Free the iterator data structure.
	stEulerTourComponentIterator_destruct(it->compIt);
	stList_destruct(it->componentList);
	free(it);
}
void stConnectedComponent_destruct(stConnectedComponent *comp) {
	free(comp);
}
//---------------------------------------------------------------------------------------------------------
struct _stEdgeContainer {
	stHash *edges;
	stSet *edgeObjects;
	void(*destructEdge)(void *);

};
struct linkedListNode {
	struct linkedListNode *next;
	struct linkedListNode *prev;
	void *key;
	void *value;
};
struct linkedList {
	struct linkedListNode *head;
};
struct linkedList *linkedList_construct() {
	struct linkedList *list = st_malloc(sizeof(struct linkedList));
	list->head = NULL;
	return list;
}
void linkedList_destruct(struct linkedList *list) {
	struct linkedListNode *node = list->head;
	if(!node) {
		free(list);
		return;
	}
	while(node->next) {
		node = node->next;
		free(node->prev);
	}
	free(node);
	free(list);
}
struct linkedListNode *linkedListNode_construct(void *key, void *value) {
	struct linkedListNode *node = st_malloc(sizeof(struct linkedListNode));
	node->next = NULL;
	node->prev = NULL;
	node->key = key;
	node->value = value;
	return node;
}
void *linkedList_search(struct linkedList *list, void *key) {
	struct linkedListNode *node = list->head;
	while(node) {
		if(node->key == key) {
			return node->value;
		}
		node = node->next;
	}
	return NULL;
}
void linkedList_delete(struct linkedList *list, void *key) {
	struct linkedListNode *node = list->head;
	while(node) {
		if(node->key == key) break;
		node = node->next;
	}
	if(!node) return;
	if(node->prev) node->prev->next = node->next;
	if(node->next) node->next->prev = node->prev;
	if(!node->prev) {
		list->head = node->next;
	}
	free(node);
}
void linkedList_insert(struct linkedList *list, void *key, void *value) {
	struct linkedListNode *new = linkedListNode_construct(key, value);
	if(!list->head) {
		list->head = new;
		return;
	}
	struct linkedListNode *node = list->head;
	while(node->next) {
		node = node->next;
	}
	node->next = new;
	new->prev = node;
	new->next = NULL;
}
stList *linkedList_toList(struct linkedList *oldList) {
	stList *newList = stList_construct();
	struct linkedListNode *node = oldList->head;
	while(node) {
		stList_append(newList, node->key);
		node = node->next;
	}
	return newList;
}

stEdgeContainer *stEdgeContainer_construct(void(*destructEdge)(void *)) {
	stEdgeContainer *container = st_malloc(sizeof(stEdgeContainer));
	container->edges = stHash_construct2(NULL, (void(*)(void*))linkedList_destruct);

	container->edgeObjects = stSet_construct2(destructEdge);
	container->destructEdge = destructEdge;
	return container;
}
void stEdgeContainer_destruct(stEdgeContainer *container) {
	stHash_destruct(container->edges);
	stSet_destruct(container->edgeObjects);
	free(container);
}
void *stEdgeContainer_getEdge(stEdgeContainer *container, void *u, void *v) {
	struct linkedList *u_incident = stHash_search(container->edges, u);
	if(!u_incident) return NULL;
	return linkedList_search(u_incident, v);
}
void stEdgeContainer_deleteEdge(stEdgeContainer *container, void *u, void *v) {
	void *edge = stEdgeContainer_getEdge(container, u, v);
	struct linkedList *u_incident = stHash_search(container->edges, u);
	linkedList_delete(u_incident, v);
	struct linkedList *v_incident = stHash_search(container->edges, v);
	linkedList_delete(v_incident, u);
	void *edgeObject = stSet_remove(container->edgeObjects, edge);
	container->destructEdge(edgeObject);
}
bool stEdgeContainer_hasEdge(stEdgeContainer *container, void *u, void *v) {
	struct linkedList *u_incident = stHash_search(container->edges, u);
	if(!u_incident) return false;
	if(linkedList_search(u_incident, v)) return true;
	return false;
}

void stEdgeContainer_addEdge(stEdgeContainer *container, void *u, void *v, void *edge) {
	struct linkedList *u_incident = stHash_search(container->edges, u);
	if(!u_incident) {
		u_incident = linkedList_construct();
		stHash_insert(container->edges, u, u_incident);
	}
	linkedList_insert(u_incident, v, edge);
	struct linkedList *v_incident = stHash_search(container->edges, v);
	if(!v_incident) {
		v_incident = linkedList_construct();
		stHash_insert(container->edges, v, v_incident);
	}
	linkedList_insert(v_incident, u, edge);

	stSet_insert(container->edgeObjects, edge);
}
stList *stEdgeContainer_getIncidentEdgeList(stEdgeContainer *container, void *v) {
	struct linkedList *vincident = stHash_search(container->edges, v);
	if(!vincident) return stList_construct();
	return linkedList_toList(vincident);
}

struct _stEdgeContainerIterator {
	stHashIterator *nodeIterator;
	struct linkedListNode *node;
	stEdgeContainer *container;
};
stEdgeContainerIterator *stEdgeContainer_getIterator(stEdgeContainer *container) {
	stEdgeContainerIterator *it = st_malloc(sizeof(stEdgeContainerIterator));
	it->nodeIterator = stHash_getIterator(container->edges);
	it->container = container;
	it->node = NULL;
	return it;
}
bool stEdgeContainer_getNext(stEdgeContainerIterator *it, void **node1, void **node2) {
	if(!it->node) {
		*node1 = stHash_getNext(it->nodeIterator);
		struct linkedList *nextList = stHash_search(it->container->edges, *node1);
		if(!nextList) return false;
		it->node = nextList->head;
	}
	if(!it->node) return false;
	*node2 = it->node->key;
	it->node = it->node->next;
	return true;
}
void stEdgeContainer_destructIterator(stEdgeContainerIterator *it) {
	stHash_destructIterator(it->nodeIterator);
	free(it);
}
