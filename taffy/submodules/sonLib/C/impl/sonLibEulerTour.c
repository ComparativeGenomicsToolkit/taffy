#include "sonLibGlobalsInternal.h"


/*represents the Euler Tour Tree of an undirected graph. The Euler Tour
 * is stored in a balanced binary tree where each node in the BST represents one 
 * traversal of an edge
 * in the graph, and the nodes are ordered by their traversal time on the Euler tour.
 * Each edge in the graph is traversed exactly twice.
 *
 *    Graph:
 *             y       z
 *        C--------D-------E
 *       / \              /
 *     w/   \x           /u
 *     /     \    v     /
 *    A       B--------F
 *    
 *    Euler tour: wxvvxyzzyw
 *
 *    Balanced Euler tour tree of the edges, keyed by traversal order in the Euler tour:
 *                       (y,6)
 *                     /       \
 *                    /         \
 *                   /           \
 *                 (v,4)         (z,8)
 *                 /   \         /   \
 *                /     \       /     \
 *              (x,2)   (x,5) (z,7)  (w,10)
 *              /   \                  /
 *             /     \                /
 *          (w,1)   (v,3)           (y,9)
 *
 *
 * Graph representation:
 *               V_right_in ____  W_right_in  ___
 *tour end  ...<---------- |    |<-----------|   |<---------
 *                         | V  |            | W |
 *tour start ...---------> |____|----------->|___|--------->
 *              V_left_out         W_left_out   
 *
 *  Undirected edges are represented as two directed half-edges.
 */

struct _stEulerVertex {
	stEulerHalfEdge *leftOut; //first time this vertex is visited in tour
	stEulerHalfEdge *rightIn; //second time vertex is visited
	stEulerTour *owner;

	void *vertexID; //index into the Euler Tree's hashtable of vertices
};

struct _stEulerHalfEdge {
	int isForwardEdge;

	stEulerVertex *from;
	stEulerVertex *to;
	stEulerHalfEdge *inverse;

	//node in treap for this half-edge.
	stTreap *node;

};

struct _stEulerTour {
	stHash *vertices;
	stEdgeContainer *edges;
	stSet *connectedComponents;

	int nComponents;
};
struct _stEulerTourIterator {
	void *currentVertex;
	stTreap *currentEdgeNode;
};
struct _stEulerTourEdgeIterator {
	stTreap *currentEdgeNode;
};
struct _stEulerTourComponentIterator {
	stSetIterator *it;
};


/*Euler Vertex methods ----------------------------------------------------- */


stEulerVertex *stEulerVertex_construct(void *vertexID) {
	stEulerVertex *v = st_malloc(sizeof(stEulerVertex));
	v->vertexID = vertexID;
	v->leftOut = NULL;
	v->rightIn = NULL;
	return(v);
}
void stEulerVertex_destruct(stEulerVertex *vertex) {
	free(vertex);
}
/*Returns the half-edge coming into this vertex on the first traversal.*/
stTreap *stEulerVertex_incidentEdgeA(stEulerVertex *vertex) {
	if(vertex->leftOut) {
		return(vertex->leftOut->node);
	}
	else {
		return(NULL);
	}
}
/*Returns the half edge leaving this vertex on the second traversal. */
stTreap *stEulerVertex_incidentEdgeB(stEulerVertex *vertex) {
	if(vertex->rightIn) {
		return(vertex->rightIn->node);
	}
	else {
		return(NULL);
	}
}
int stEulerVertex_connected(stEulerVertex *from, stEulerVertex *to) {
	if(from == to) {
		return(true);
	}
	stTreap *fromNode = stEulerVertex_incidentEdgeA(from);
	stTreap *toNode = stEulerVertex_incidentEdgeA(to);
	if(fromNode == NULL || toNode == NULL) {
		return(false);
	}
	stTreap *a_root = stTreap_findRoot(fromNode);
	stTreap *b_root = stTreap_findRoot(toNode);
	return(a_root == b_root);
}
int stEulerVertex_isSingleton(stEulerVertex *vertex) {
	assert(!vertex->rightIn == !vertex->leftOut);

	return(!(vertex->leftOut));
}

/*Euler Half Edge methods ------------------------------------------------- */
stEulerHalfEdge *stEulerHalfEdge_construct() {
	stEulerHalfEdge *newEdge = st_malloc(sizeof(stEulerHalfEdge));
	newEdge->node = stTreap_construct(newEdge);
	return newEdge;
}


int stEulerHalfEdge_contains(stEulerHalfEdge *edge, stEulerVertex *vertex) {
	return (edge->from == vertex || edge->to == vertex);
}
void stEulerHalfEdge_destruct(stEulerHalfEdge *edge) {
	if(edge != NULL) {
		stTreap_nodeDestruct(edge->node);
	}
	if(edge->inverse != NULL) {
		stTreap_nodeDestruct(edge->inverse->node);
	}
	free(edge->inverse);
	free(edge);
}
void *stEulerHalfEdge_getFrom(stEulerHalfEdge *edge) {
	return(edge->from->vertexID);
}
void *stEulerHalfEdge_getTo(stEulerHalfEdge *edge) {
	return(edge->to->vertexID);
}


/*Euler Tour methods ------------------------------------------------------ */
stEulerTour *stEulerTour_construct() {
	stEulerTour *et = st_malloc(sizeof(stEulerTour));

	et->vertices = stHash_construct2(NULL, (void(*)(void*))stEulerVertex_destruct);
	et->edges = stEdgeContainer_construct((void(*)(void*))stEulerHalfEdge_destruct);
	et->connectedComponents = stSet_construct();

	et->nComponents = 0;
	return(et);
}
bool stEulerTour_hasEdge(stEulerTour *et, void *u, void *v) {
	return stEdgeContainer_hasEdge(et->edges, u, v);
}
int stEulerTour_getNComponents(stEulerTour *et) {
    return stSet_size(et->connectedComponents);
}
stEdgeContainer *stEulerTour_getEdges(stEulerTour *et) {
	return et->edges;
}
bool stEulerTour_isSingleton(stEulerTour *et, void *v) {
	stEulerVertex *vertex = stEulerTour_getVertex(et, v);
	return stEulerVertex_isSingleton(vertex);
}
void stEulerTour_printTour(stEulerTour *et, void *v) {
	stEulerTourIterator *it = stEulerTour_getIterator(et, v);
	printf("tour: ");
	void *node;
	while((node = stEulerTourIterator_getNext(it))) {
		printf("%p ", node);
	}
	stEulerTourIterator_destruct(it);
	printf("\n");
}
void stEulerTour_printEdgeTour(stEulerTour *et, void *v) {
	stEulerTourEdgeIterator *it = stEulerTour_getEdgeIterator(et, v);
	printf("edge tour: ");
	void *node1, *node2;
	while(stEulerTourEdgeIterator_getNext(it, &node1, &node2)) {
		printf("%p->%p  ", node1, node2);
	}
	stEulerTourEdgeIterator_destruct(it);
}

void stEulerTour_destruct(stEulerTour *et) {
	if(et != NULL) {
		stHash_destruct(et->vertices);
		stEdgeContainer_destruct(et->edges);
		stSet_destruct(et->connectedComponents);
		free(et);
	}
}
stTreap *stEulerTour_findRoot(stEulerTour *et, void *v) {
	stEulerVertex *vertex = stEulerTour_getVertex(et, v);
	if(!vertex) {
		fprintf(stderr, "Warning: vertex does not exist\n");
		return(NULL);
	}
	stTreap *treapNode = stEulerVertex_incidentEdgeA(vertex);
	if(!treapNode) {
		return(NULL);
	}
	return(stTreap_findMin(stTreap_findRoot(treapNode)));
}
void *stEulerTour_findRootNode(stEulerTour *et, void *v) {
	stTreap *root = stEulerTour_findRoot(et, v);
	stEulerHalfEdge *edge = (stEulerHalfEdge *) stTreap_getValue(root);
	return(edge->from->vertexID);
}
int stEulerTour_connected(stEulerTour *et, void *u, void *v) {
	stEulerVertex *node1 = stHash_search(et->vertices, u);
	stEulerVertex *node2 = stHash_search(et->vertices, v);
	return(stEulerVertex_connected(node1, node2));
}

/*returns the size of the component that vertex v is in. */
int stEulerTour_size(stEulerTour *et, void *v) {
	stEulerVertex *vertex = stHash_search(et->vertices, v);
	if(!vertex->leftOut) {
		return(1);
	}

	int tour_length = stTreap_size(vertex->leftOut->node);
	return (tour_length/2 + 1);
}
stEulerVertex *stEulerTour_getVertex(stEulerTour *et, void *v) {
	return(stHash_search(et->vertices, v));
}
stEulerVertex *stEulerTour_createVertex(stEulerTour *et, void *vertexID) {
	stEulerVertex *newVertex = stEulerVertex_construct(vertexID);
	newVertex->owner = et;
	et->nComponents++;
	stSet_insert(et->connectedComponents, vertexID);
	stHash_insert(et->vertices, vertexID, newVertex);
	return(newVertex);
}
void stEulerTour_removeVertex(stEulerTour *et, void *vertexID) {
	assert(stEulerTour_isSingleton(et, vertexID)); //remove all edges before removing the vertex
	stEulerVertex *removedVertex = stHash_remove(et->vertices, vertexID);
	stEulerVertex_destruct(removedVertex);
	stSet_remove(et->connectedComponents, vertexID);
	et->nComponents--;
}

void *stEulerTour_getConnectedComponent(stEulerTour *et, void *nodeInComponent) {
	if(stEulerTour_isSingleton(et, nodeInComponent)) {
		return nodeInComponent;
	}
	return stEulerTour_findRootNode(et, nodeInComponent);
}

/*Move a vertex to the beginning of the tour. */
void stEulerTour_makeRoot(stEulerTour *et, stEulerVertex *vertex) {
	if(stEulerVertex_isSingleton(vertex)) {
		return;
	}
	//update this component's root node in the list of connected components
	void *currentRoot = stEulerTour_findRootNode(et, vertex->vertexID);
	stSet_remove(et->connectedComponents, currentRoot);
	stSet_insert(et->connectedComponents, vertex->vertexID);

	
	if(stTreap_size(vertex->leftOut->node) == 2) {
		assert(stTreap_findRoot(stEulerVertex_incidentEdgeA(vertex)) ==
				stTreap_findRoot(stEulerVertex_incidentEdgeB(vertex)));
		stTreap *first = stTreap_findMin(stTreap_findRoot(stEulerVertex_incidentEdgeA(vertex)));
		stEulerHalfEdge *firstEdge = stTreap_getValue(first);
		if(firstEdge->from == vertex) return;
		stTreap *second = stTreap_splitAfter(first);
		stTreap_concat(second, first);
		return;
	}


	stTreap *ftreapNode = stEulerVertex_incidentEdgeA(vertex);
	stTreap *btreapNode = stEulerVertex_incidentEdgeB(vertex);
	stEulerHalfEdge *f = stTreap_getValue(ftreapNode);
	stEulerHalfEdge *b = stTreap_getValue(btreapNode);
	assert(f != b);
	assert(f);
	assert(b);
	if(stTreap_compare(f->node, b->node) > 0) {
		//swap the pointers to the half edges;
		stEulerHalfEdge *temp = f;
		f = b;
		b = temp;
	}
	assert(stTreap_compare(f->node, b->node) < 0);

	stEulerVertex *other = (f->to == vertex) ? f->from: f->to;
	assert(other != vertex);

	//the next edge traversed in the Euler tour
	stTreap *nextNode = stTreap_next(f->node);
	stEulerHalfEdge *next = (stEulerHalfEdge *)stTreap_getValue(nextNode);

	if(!stEulerHalfEdge_contains(next, vertex)) {
		stTreap *fnodeprev = stTreap_prev(f->node);

		//if the previous node is null, vertex comes first in the tour (is the root)
		//so do nothing.
		if (!fnodeprev) {
			assert(stTreap_findRoot(stEulerVertex_incidentEdgeA(vertex)) == 
					stTreap_findRoot(stEulerVertex_incidentEdgeB(vertex)));
			return;
		}
		else {
			f = stTreap_getValue(fnodeprev);
		}
	}
	else if (stEulerHalfEdge_contains(next, other)) {
		stTreap *next_next = stTreap_next(next->node);	
		if(!next_next) {
			next_next = stTreap_prev(f->node);
		}

		if (next_next) {
			stEulerHalfEdge *next_next_edge = 
				(stEulerHalfEdge*)stTreap_getValue(next_next);
			if(stEulerHalfEdge_contains(next_next_edge, vertex)) {
				f = next;
			} else {
				//node is leaf
			}
		}
	} else {
		//next contains node
	}

	//f is now guaranteed to be before b in the tour
	//assert(stTreap_compare(f->node, b->node) < 0);
	stTreap *rightSubtree = stTreap_splitAfter(f->node);

	if(rightSubtree) {
		assert(stTreap_findMax(stTreap_findRoot(f->node)) == f->node);
		assert(stTreap_findRoot(rightSubtree) != stTreap_findRoot(f->node));
		stTreap_concat(rightSubtree, f->node);
	}
	assert(stTreap_findMax(stTreap_findRoot(f->node)) == f->node);
	assert(stTreap_findRoot(stEulerVertex_incidentEdgeA(vertex)) == 
			stTreap_findRoot(stEulerVertex_incidentEdgeB(vertex)));
}
void stEulerTour_link(stEulerTour *et, void *u, void *v) {
	assert(u != v);
	assert(!stEulerTour_connected(et, u, v));
	assert(!stEdgeContainer_hasEdge(et->edges, u, v));
	stEulerVertex *vertex = stHash_search(et->vertices, u);
	stEulerVertex *other = stHash_search(et->vertices, v);
	et->nComponents--;
	

	stEulerHalfEdge *newForwardEdge = stEulerHalfEdge_construct();
	stEulerHalfEdge *newBackwardEdge = stEulerHalfEdge_construct();
	stEdgeContainer_addEdge(et->edges, u, v, newForwardEdge);

	newForwardEdge->isForwardEdge = true;
	newBackwardEdge->isForwardEdge = false;
	newBackwardEdge->inverse = newForwardEdge;
	newForwardEdge->inverse = newBackwardEdge;

	newForwardEdge->from = vertex;
	newForwardEdge->to = other;

	newBackwardEdge->from = other;
	newBackwardEdge->to = vertex;

	
	stEulerTour_makeRoot(et, vertex);
	stEulerTour_makeRoot(et, other);

	stSet_remove(et->connectedComponents, stEulerTour_getConnectedComponent(et, vertex->vertexID));
	stSet_remove(et->connectedComponents, stEulerTour_getConnectedComponent(et, other->vertexID));

	stTreap *f = NULL;
	if (stEulerVertex_incidentEdgeA(vertex)) {
		f = stTreap_findMin(stTreap_findRoot(stEulerVertex_incidentEdgeA(vertex)));
	}


	stTreap *tleft = NULL;
	if (f) {
		stTreap_concat(f, newForwardEdge->node);
	}
	else {
		vertex->leftOut = newForwardEdge;
	}
	if (stEulerVertex_incidentEdgeA(other)) {
		assert(stTreap_findRoot(stEulerVertex_incidentEdgeA(other)) == 
				stTreap_findRoot(stEulerVertex_incidentEdgeB(other)));
		stTreap_concat(newForwardEdge->node, other->leftOut->node);
		assert(stTreap_findRoot(stEulerVertex_incidentEdgeA(other)) == 
				stTreap_findRoot(stEulerVertex_incidentEdgeB(other)));
		assert(stTreap_findRoot(newForwardEdge->node) == stTreap_findRoot(stEulerVertex_incidentEdgeB(other)));
	}
	else {
		other->leftOut = newForwardEdge;
	}
	if (stEulerVertex_incidentEdgeB(other)) {
		assert(stTreap_findRoot(stEulerVertex_incidentEdgeA(other)) == 
				stTreap_findRoot(stEulerVertex_incidentEdgeB(other)));
		stTreap_concat(other->rightIn->node, newBackwardEdge->node);
		assert(stTreap_findRoot(stEulerVertex_incidentEdgeA(other)) == 
				stTreap_findRoot(newBackwardEdge->node));
		assert(stTreap_findRoot(stEulerVertex_incidentEdgeA(vertex)) == 
				stTreap_findRoot(newBackwardEdge->node));

	}
	else {
		other->rightIn = newBackwardEdge;
		stTreap_concat(stEulerVertex_incidentEdgeA(vertex), newBackwardEdge->node);
	}
	if (tleft) {
		stTreap_concat(newBackwardEdge->node, tleft);
		vertex->rightIn = stTreap_getValue(tleft);
	}
	else {
		vertex->rightIn = newBackwardEdge;
	}
	assert(stTreap_findRoot(newForwardEdge->node) == stTreap_findRoot(newBackwardEdge->node));
	stSet_insert(et->connectedComponents, stEulerTour_getConnectedComponent(et, vertex->vertexID));
	assert(stEulerVertex_connected(vertex, other));
	assert(stEulerTour_hasEdge(et, u, v));
	assert(stEulerTour_hasEdge(et, v, u));
}


/*remove an edge from the graph, splitting into two Euler tours.
 *
 * Graph:
 *           y1        z1
 *	   C--------D----------E
 *        / \ y2       z2     /
 *       /   \               /               
 *    w1/w2 x1\x2         u1/u2
 *     /       \     v1    /
 *    A         B---------F
 *                   v2
 *
 * Euler tour: w1->x1->v1->v2->x2->y1->z1->z2->y2->w2
 *                       y1
 *                    /      \
 *                   /        \
 *                  v2        z2
 *                 /   \     /   \
 *                /     \   /     \
 *               x1    x2  z1     w2
 *              /   \             /
 *             /     \           /
 *          w1       v1         y2
 *
 * Remove y:
 * 
 * tree1: w1->x1->v1->v2->x2  tree2: w2
 *
 * concat(tree1, tree2): w1 x1 v1 v2 x2 w2
 *          
 *
 */

void stEulerTour_cut(stEulerTour *et, void *u, void *v) {
	et->nComponents++;
	stSet_remove(et->connectedComponents, stEulerTour_getConnectedComponent(et, u));

	//get the two halves of this edge
	assert(stEulerTour_connected(et, u, v));
	assert(stEulerTour_hasEdge(et, u, v));
	assert(stEulerTour_hasEdge(et, v, u));
	stEulerHalfEdge *f = stEdgeContainer_getEdge(et->edges, u, v);
	assert(f);
	stEulerHalfEdge *b = f->inverse;
	assert(b->inverse == f);
	assert(b->inverse->inverse == b);
	stEulerVertex *from = f->from;
	stEulerVertex *to = f->to;

	assert(stEulerTour_connected(et, f->from->vertexID, f->to->vertexID));
	assert(stTreap_findRoot(f->node) == stTreap_findRoot(b->node));

	//f should point to the first half of this edge to be traversed,
	//and b should point to the half that is traversed second.
	if (stTreap_compare(f->node, b->node) > 0) {
		//swap f and b
		stEulerHalfEdge *temp = f;
		f = b;
		b = temp;
	}


	stTreap *p = stTreap_prev(f->node); //p
	stTreap *n = stTreap_next(b->node); //n
	stTreap *pn = stTreap_next(f->node); //pn
	stTreap *nn = stTreap_prev(b->node); //nn

	assert(stTreap_next(nn) == b->node);

	/*              p        ____      f       ____   pn
	 *tour start ---------->|from|----------->|to  |-------->
	 *tour end   <----------|____|<-----------|____|<--------
	 *              n                  b              nn  
	 *                          
	 *                          
	 *                  
	 *
	 */

	//get the part of the tour that happened before reaching the first half of this
	//edge and the part that happened after reaching the second half.
	stTreap *tree1 = stTreap_splitBefore(f->node);
	assert(!tree1 || stTreap_findMax(tree1) == p);
	if(tree1) {
		assert(stTreap_findRoot(tree1) != stTreap_findRoot(f->node));
	}
	stTreap *tree2 = stTreap_splitAfter(b->node);
	assert(!tree2 || stTreap_findMin(tree2) == n);

	if (tree1 && tree2) {
		stTreap_concat(tree1, tree2);
	}
	assert(!tree1 || stTreap_findRoot(tree1) != stTreap_findRoot(f->node));
	assert(!tree2 || stTreap_findRoot(tree2) != stTreap_findRoot(b->node));
	assert(!tree1 || stTreap_findRoot(tree1) != stTreap_findRoot(b->node));
	assert(!tree2 || stTreap_findRoot(tree2) != stTreap_findRoot(f->node));

	assert(pn);
	assert(nn);

	//choose new incident edges
	if (stEulerHalfEdge_contains(stTreap_getValue(pn), from) && 
			stEulerHalfEdge_contains(stTreap_getValue(pn), to)) {

		if((n || p) && !(n && p)) {
			if (!n) {
				n = stTreap_findMin(stTreap_findRoot(p));
			}
			else {
				p = stTreap_findMax(stTreap_findRoot(n));
			}
			assert(n != p);
		}
		if (n) {
			assert(p);
			assert(stEulerHalfEdge_contains(stTreap_getValue(n), from) || 
					stEulerHalfEdge_contains(stTreap_getValue(n), to));
			assert(stEulerHalfEdge_contains(stTreap_getValue(p), from) || 
					stEulerHalfEdge_contains(stTreap_getValue(p), to));
			if(stEulerHalfEdge_contains(stTreap_getValue(n), from)) {
				assert(stEulerHalfEdge_contains(stTreap_getValue(p), from));
				from->leftOut = stTreap_getValue(n);
				from->rightIn = stTreap_getValue(p);
				to->leftOut = NULL;
				to->rightIn = NULL;
			}
			else {
				assert(stEulerHalfEdge_contains(stTreap_getValue(n), to));
				assert(stEulerHalfEdge_contains(stTreap_getValue(p), to));
				to->leftOut = stTreap_getValue(n);
				to->rightIn = stTreap_getValue(p);
				from->leftOut = NULL;
				from->rightIn = NULL;
			}
		}
		else {
			//to and from are singletons
			from->leftOut = NULL;
			from->rightIn = NULL;
			to->leftOut = NULL;
			to->rightIn = NULL;
		}	
	}
	else if (stEulerHalfEdge_contains(stTreap_getValue(pn), from)) {
		assert(stEulerHalfEdge_contains(stTreap_getValue(nn), from));
		from->leftOut = stTreap_getValue(pn);
		from->rightIn = stTreap_getValue(nn);
		if ((n || p) && !(n && p)) {
			if(!n) {
				n = stTreap_findMin(stTreap_findRoot(p));
			}
			else {
				p = stTreap_findMax(stTreap_findRoot(n));
			}
			assert(n != p);
		}
		if (n) {
			assert(p);
			assert(stEulerHalfEdge_contains(stTreap_getValue(n), to));
			assert(stEulerHalfEdge_contains(stTreap_getValue(p), to));
			to->leftOut = stTreap_getValue(n);
			to->rightIn = stTreap_getValue(p);
		}
		else {
			//to is a singleton
			to->leftOut = NULL;
			to->rightIn = NULL;
		}
	}
	else if (stEulerHalfEdge_contains(stTreap_getValue(pn), to)) {
		assert(stEulerHalfEdge_contains(stTreap_getValue(nn), to));
		to->leftOut = stTreap_getValue(pn);
		to->rightIn = stTreap_getValue(nn);
		if( (n || p) && !(n && p)) {
			if(!n) {
				n = stTreap_findMin(stTreap_findRoot(p));
			}
			else {
				p = stTreap_findMax(stTreap_findRoot(n));
			}
			assert(n != p);
		}
		if (n) {
			assert(p);
			assert(stEulerHalfEdge_contains(stTreap_getValue(n), from));
			assert(stEulerHalfEdge_contains(stTreap_getValue(p), from));
			from->leftOut = stTreap_getValue(n);
			from->rightIn = stTreap_getValue(p);
		}
		else {
			//to is a singleton
			from->leftOut = NULL;
			from->rightIn = NULL;
		}
	}
	stTreap *after = stTreap_splitAfter(f->node);
	stTreap *before = stTreap_splitBefore(b->node);
	if(before) assert(stTreap_next(before) != b->node);
	if(after) assert(stTreap_prev(after) != f->node);

	if(stEulerVertex_incidentEdgeA(from) && 
			(stTreap_size(stEulerVertex_incidentEdgeA(from)) == 1)) {
		//from is a singleton
		from->leftOut = NULL;
		from->rightIn = NULL;
	}
	if(stEulerVertex_incidentEdgeA(to) && 
			(stTreap_size(stEulerVertex_incidentEdgeA(to)) == 1)) {
		to->leftOut = NULL;
		to->rightIn = NULL;
	}
	assert(stTreap_size(f->node) == 1);
	assert(stTreap_size(b->node) == 1);
	assert(stTreap_next(f->node) == NULL);
	assert(stTreap_prev(f->node) == NULL);
	assert(stTreap_next(b->node) == NULL);
	assert(stTreap_prev(b->node) == NULL);
	stEdgeContainer_deleteEdge(et->edges, u, v);

	stSet_insert(et->connectedComponents, stEulerTour_getConnectedComponent(et, u));
	stSet_insert(et->connectedComponents, stEulerTour_getConnectedComponent(et, v));
}
//------------------------------------------------------------------
// Tour iterators
stEulerTourIterator *stEulerTour_getIterator(stEulerTour *et, void *v) {
	stEulerTourIterator *it = st_malloc(sizeof(stEulerTourIterator));
	it->currentEdgeNode = stEulerTour_findRoot(et, v);
	it->currentVertex = v;
	return(it);
}
void *stEulerTourIterator_getNext(stEulerTourIterator *it) {
	if(!it->currentEdgeNode) {
		void *currentVertex = it->currentVertex;
		it->currentVertex = NULL;
		return(currentVertex);
	}
	stEulerHalfEdge *currentEdge = (stEulerHalfEdge *) stTreap_getValue(it->currentEdgeNode);
	void *vertexToReturn = currentEdge->from->vertexID;
	it->currentVertex = currentEdge->to->vertexID;

	it->currentEdgeNode = stTreap_next(it->currentEdgeNode);
	return(vertexToReturn);
}
void stEulerTourIterator_destruct(stEulerTourIterator *it) {
	free(it);
}
stSet *stEulerTour_getNodesInComponent(stEulerTour *et, void *v) {
	stEulerTourIterator *it = stEulerTour_getIterator(et, v);
	stSet *nodes = stSet_construct();
	void *node;
	while((node = stEulerTourIterator_getNext(it))) {
		stSet_insert(nodes, node);
	}
	stEulerTourIterator_destruct(it);
	return(nodes);
}
stEulerTourEdgeIterator *stEulerTour_getEdgeIterator(stEulerTour *et, void *v) {
	stEulerTourEdgeIterator *it = st_malloc(sizeof(stEulerTourEdgeIterator));
	it->currentEdgeNode = stEulerTour_findRoot(et, v);
	return(it);
}
bool stEulerTourEdgeIterator_getNext(stEulerTourEdgeIterator *it, void **node1, void **node2) {
	if(!it->currentEdgeNode) return false;
	stEulerHalfEdge *edgeToReturn = (stEulerHalfEdge *)stTreap_getValue(it->currentEdgeNode);
	it->currentEdgeNode = stTreap_next(it->currentEdgeNode);
	*node1 = edgeToReturn->from->vertexID;
	*node2 = edgeToReturn->to->vertexID;
	return true;
}
void stEulerTourEdgeIterator_destruct(stEulerTourEdgeIterator *it) {
	free(it);
}
stEulerTourComponentIterator *stEulerTour_getComponentIterator(stEulerTour *et) {
	stEulerTourComponentIterator *it = st_malloc(sizeof(stEulerTourComponentIterator));
	it->it = stSet_getIterator(et->connectedComponents);
	return it;
}
void *stEulerTourComponentIterator_getNext(stEulerTourComponentIterator *it) {
	return(stSet_getNext(it->it));
}
void stEulerTourComponentIterator_destruct(stEulerTourComponentIterator *it) {
	stSet_destructIterator(it->it);
	free(it);
}
