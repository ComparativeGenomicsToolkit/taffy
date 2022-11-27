// Naively keeps track of the connected components in an
// undirected graph. Very slow and wasteful.
#include "sonLibGlobalsInternal.h"

// Data structures

struct _stNaiveConnectivity {
    // Hash to an adjacency list. Edges are only kept in the adjacency
    // list belonging the node whose pointer value is smallest.
    stHash *nodesToAdjList;
    // What it says on the tin.
    stNaiveConnectedComponent *connectedComponentCache;
};

struct _stNaiveConnectedComponent {
    stSet *nodes;
    stNaiveConnectedComponent *next;
};

struct _stNaiveConnectedComponentIterator {
    stNaiveConnectedComponent *cur;
};

struct _stNaiveConnectedComponentNodeIterator {
    stSetIterator *it;
};

// Linked list for storing edges.
struct adjacency {
    void *toNode;
    // for faster removal of nodes/edges. maybe not worth it, but who
    // cares?
    struct adjacency *inverse;
    struct adjacency *prev;
    struct adjacency *next;
};

// Destroy an adjacency linked list.
static void adjList_destruct(struct adjacency *adjList) {
    assert(adjList == NULL || adjList->prev == NULL);
    while (adjList != NULL) {
        struct adjacency *next = adjList->next;
        free(adjList);
        adjList = next;
    }
}
stSet *stNaiveConnectedComponent_getNodes(stNaiveConnectedComponent *comp) {
	return(comp->nodes);
}

stNaiveConnectivity *stNaiveConnectivity_construct(void) {
    stNaiveConnectivity *connectivity = calloc(1, sizeof(struct _stNaiveConnectivity));
    connectivity->nodesToAdjList = stHash_construct2(NULL, (void (*)(void *)) adjList_destruct);

    return connectivity;
}

static void invalidateCache(stNaiveConnectivity *connectivity) {
    if (connectivity->connectedComponentCache != NULL) {
        stNaiveConnectedComponent *curComponent = connectivity->connectedComponentCache;
        while (curComponent != NULL) {
            stSet_destruct(curComponent->nodes);
            stNaiveConnectedComponent *next = curComponent->next;
            free(curComponent);
            curComponent = next;
        }
        connectivity->connectedComponentCache = NULL;
    }
}

void stNaiveConnectivity_destruct(stNaiveConnectivity *connectivity) {
    invalidateCache(connectivity);
    stHash_destruct(connectivity->nodesToAdjList);
    free(connectivity);
}

void stNaiveConnectivity_addNode(stNaiveConnectivity *connectivity, void *node) {
    invalidateCache(connectivity);

    assert(stHash_search(connectivity->nodesToAdjList, node) == NULL);
    stHash_insert(connectivity->nodesToAdjList, node, NULL);
}

void stNaiveConnectivity_addEdge(stNaiveConnectivity *connectivity, void *node1, void *node2) {
    invalidateCache(connectivity);

    struct adjacency *newEdge1 = malloc(sizeof(struct adjacency));
    struct adjacency *newEdge2 = malloc(sizeof(struct adjacency));
    newEdge1->toNode = node2;
    newEdge2->toNode = node1;
    newEdge1->inverse = newEdge2;
    newEdge2->inverse = newEdge1;
    newEdge1->prev = NULL;
    newEdge2->prev = NULL;

    struct adjacency *adjList1 = stHash_search(connectivity->nodesToAdjList, node1);
    if (adjList1 == NULL) {
        newEdge1->next = NULL;
    } else {
        newEdge1->next = adjList1;
        adjList1->prev = newEdge1;
    }
    stHash_remove(connectivity->nodesToAdjList, node1);
    stHash_insert(connectivity->nodesToAdjList, node1, newEdge1);

    struct adjacency *adjList2 = stHash_search(connectivity->nodesToAdjList, node2);
    if (adjList2 == NULL) {
        newEdge2->next = NULL;
    } else {
        newEdge2->next = adjList2;
        adjList2->prev = newEdge2;
    }
    stHash_remove(connectivity->nodesToAdjList, node2);
    stHash_insert(connectivity->nodesToAdjList, node2, newEdge2);
}

// Remove and free an edge properly.
static void removeEdgeFromAdjList(stNaiveConnectivity *connectivity, void *node, struct adjacency *adj) {
    invalidateCache(connectivity);

    if (adj->next != NULL) {
        adj->next->prev = adj->prev;
    }

    if (adj->prev != NULL) {
        adj->prev->next = adj->next;
    } else {
        stHash_remove(connectivity->nodesToAdjList, node);
        stHash_insert(connectivity->nodesToAdjList, node, adj->next);
    }
    free(adj);
}
bool stNaiveConnectivity_hasEdge(stNaiveConnectivity *connectivity, void *node1, void *node2) {
	struct adjacency *adjList1 = stHash_search(connectivity->nodesToAdjList, node1);
	if(!adjList1) return false;
	while(adjList1 != NULL) {
		if (adjList1->toNode == node2) {
			return true;
		}
		adjList1 = adjList1->next;
	}
	return false;
}

void stNaiveConnectivity_removeEdge(stNaiveConnectivity *connectivity, void *node1, void *node2) {
    invalidateCache(connectivity);

    struct adjacency *adjList1 = stHash_search(connectivity->nodesToAdjList, node1);
    assert(adjList1 != NULL);
    assert(stHash_search(connectivity->nodesToAdjList, node2) != NULL);

    while (adjList1 != NULL) {
        if (adjList1->toNode == node2) {
            break;
        }
        adjList1 = adjList1->next;
    }

    // We can find the link in the other node's adjacency list easily
    struct adjacency *adjList2 = adjList1->inverse;
    assert(adjList2->inverse == adjList1);

    // Now remove the links from the lists, and free them.
    removeEdgeFromAdjList(connectivity, node1, adjList1);
    removeEdgeFromAdjList(connectivity, node2, adjList2);
}

void stNaiveConnectivity_removeNode(stNaiveConnectivity *connectivity, void *node) {
    invalidateCache(connectivity);

    struct adjacency *adjList = stHash_search(connectivity->nodesToAdjList, node);

    while (adjList != NULL) {
        struct adjacency *next = adjList->next; // Have to do this beforehand -- adjList will be freed by next line!
        stNaiveConnectivity_removeEdge(connectivity, node, adjList->toNode);
        adjList = next;
    }

    stHash_remove(connectivity->nodesToAdjList, node);
}

// Remove a connected component from the list of connected components
// properly, then free it.
static void removeComponent(stNaiveConnectedComponent **head, stNaiveConnectedComponent *component) {
    stNaiveConnectedComponent *curComponent = *head;
    assert(curComponent != NULL);
    stNaiveConnectedComponent *prevComponent = NULL;
    while (curComponent != NULL) {
        if (curComponent == component) {
            if (prevComponent != NULL) {
                prevComponent->next = curComponent->next;
            } else {
                *head = curComponent->next;
            }
            stSet_destruct(component->nodes);
            free(component);
            break;
        }
        prevComponent = curComponent;
        curComponent = curComponent->next;
    }
}

// Compute the connected components, if they haven't been computed
// already since the last modification.
static void computeConnectedComponents(stNaiveConnectivity *connectivity) {
    if (connectivity->connectedComponentCache != NULL) {
        // Already computed the connected components.
        return;
    }

    stHashIterator *nodeIt = stHash_getIterator(connectivity->nodesToAdjList);
    void *node;
    stNaiveConnectedComponent *componentsHead = NULL;
    while ((node = stHash_getNext(nodeIt)) != NULL) {
        stSet *myNodeSet = stSet_construct();
        stSet_insert(myNodeSet, node);
        struct adjacency *adjList = stHash_search(connectivity->nodesToAdjList, node);
        if (adjList != NULL) {
            while (adjList != NULL) {
                stSet_insert(myNodeSet, adjList->toNode);
                adjList = adjList->next;
            }
        }

        // Now go through the existing connected components and see if
        // this overlaps any of them. If it's not a full overlap, then
        // this set becomes the union, and we continue looking for
        // additional overlaps, then this becomes a new connected
        // component. If we find that this is a subset of an existing
        // component, we can quit early, since we can't possibly add
        // to it or any others.
        stNaiveConnectedComponent *curComponent = componentsHead;
        while (curComponent != NULL) {
            stNaiveConnectedComponent *next = curComponent->next;

            // Find out whether our node set is a subset of this
            // connected component, or if it shares any overlap.
            bool isSubset = true;
            bool overlap = false;
            stSetIterator *myNodeIt = stSet_getIterator(myNodeSet);
            void *node;
            while ((node = stSet_getNext(myNodeIt)) != NULL) {
                if (stSet_search(curComponent->nodes, node)) {
                    overlap = true;
                } else {
                    isSubset = false;
                }
            }
            stSet_destructIterator(myNodeIt);

            if (isSubset) {
                assert(overlap == true);
                // Quit early.
                stSet_destruct(myNodeSet);
                myNodeSet = NULL;
                break;
            } else if (overlap) {
                stSet *newNodeSet = stSet_getUnion(myNodeSet, curComponent->nodes);
                stSet_destruct(myNodeSet);
                removeComponent(&componentsHead, curComponent);
                myNodeSet = newNodeSet;
            }

            curComponent = next;
        }
        if (myNodeSet != NULL) {
            // We have a new (or possibly merged) connected component to
            // add to the list.
            stNaiveConnectedComponent *newComponent = malloc(sizeof(stNaiveConnectedComponent));
            newComponent->nodes = myNodeSet;
            newComponent->next = componentsHead;
            componentsHead = newComponent;
        }
    }

    stHash_destructIterator(nodeIt);

    connectivity->connectedComponentCache = componentsHead;
}

stNaiveConnectedComponent *stNaiveConnectivity_getConnectedComponent(stNaiveConnectivity 
		*connectivity, void *node) {
    computeConnectedComponents(connectivity);

    stNaiveConnectedComponent *curComponent = connectivity->connectedComponentCache;
    while (curComponent != NULL) {
        if (stSet_search(curComponent->nodes, node)) {
            break;
        }
        curComponent = curComponent->next;
    }
    assert(curComponent != NULL);
    return curComponent;
}

stNaiveConnectedComponentNodeIterator *stNaiveConnectedComponent_getNodeIterator(stNaiveConnectedComponent 
		*component) {
    stNaiveConnectedComponentNodeIterator *ret = malloc(sizeof(stNaiveConnectedComponentNodeIterator));
    ret->it = stSet_getIterator(component->nodes);
    return ret;
}

void *stNaiveConnectedComponentNodeIterator_getNext(stNaiveConnectedComponentNodeIterator *it) {
    return stSet_getNext(it->it);
}

void stNaiveConnectedComponentNodeIterator_destruct(stNaiveConnectedComponentNodeIterator *it) {
    stSet_destructIterator(it->it);
    free(it);
}

stNaiveConnectedComponentIterator *stNaiveConnectivity_getConnectedComponentIterator(stNaiveConnectivity 
		*connectivity) {
    computeConnectedComponents(connectivity);

    stNaiveConnectedComponentIterator *ret = malloc(sizeof(stNaiveConnectedComponentIterator));
    ret->cur = connectivity->connectedComponentCache;
    return ret;
}

stNaiveConnectedComponent *stNaiveConnectedComponentIterator_getNext(stNaiveConnectedComponentIterator *it) {
    stNaiveConnectedComponent *ret = it->cur;
    if (ret != NULL) {
        it->cur = it->cur->next;
    }
    return ret;
}

void stNaiveConnectedComponentIterator_destruct(stNaiveConnectedComponentIterator *it) {
    free(it);
}
