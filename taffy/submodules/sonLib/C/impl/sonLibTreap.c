#include "sonLibGlobalsInternal.h"
/*Implementation of a balanced binary tree. The binary tree property is maintained with 
 * respect to the keys, while the heap property is maintained for the priority values, which
 * are generated randomly for each element as it is inserted into the tree. The randomized
 * priority values make the treap highly likely to be balanced. */


stTreap *stTreap_construct(void *value) {
	//strand(time(0));
	stTreap *node = st_malloc(sizeof(stTreap));
	node->key = 0;
	node->priority = rand();
	node->left = node->right = node->parent = NULL;
	node->count = 1;
	node->value = value;
	return(node);
}
void *stTreap_getValue(stTreap *node) {
	return(node->value);
}
void stTreap_destructRecurse(stTreap *root) {
	if(root->left) {
		stTreap_destructRecurse(root->left);
	}
	if(root->right) {
		stTreap_destructRecurse(root->right);
	}
	//free(root->value);
	free(root);
}
void stTreap_destruct(stTreap *node) {
	node = stTreap_findRoot(node);
	stTreap_destructRecurse(node);
}
void stTreap_nodeDestruct(stTreap *node) {
	free(node);
}
char *stTreap_print(stTreap *node) {
	node = stTreap_findRoot(node);
	node = stTreap_findMin(node);
	char *path = st_calloc(stTreap_size(node) + 1, 1);
	while(node) {
		strcat(path, node->value);
		node = stTreap_next(node);
	}
	return(path);
}
char *stTreap_printBackwards(stTreap *node) {
	node = stTreap_findRoot(node);
	node = stTreap_findMax(node);
	char *path = st_calloc(stTreap_size(node) + 1, 1);
	while(node) {
		strcat(path, node->value);
		node = stTreap_prev(node);
	}
	return(path);
}


/*Left tree rotation at node u, which preserves the in-order traversal:
 *      
 *        p                  p
 *       /                  /
 *      u                  w
 *     / \     ----->     / \
 *    a   w              u   c
 *       / \            / \
 *      b   c          a   b
 * 
 */
void stTreap_rotateLeft(stTreap *u) {
	stTreap *w = u->right;
	if(!w) return;
	w->parent = u->parent;
	if(w->parent) {
		if(w->parent->left == u) {
			w->parent->left = w;
		}
		else {
			w->parent->right = w;
		}
	}
	u->right = w->left;
	if(u->right) {
		u->right->parent = u;
	}
	u->parent = w;
	w->left = u;


	u->count -= w->count;
	if(u->right) {
		u->count += u->right->count;
		w->count -= u->right->count;
	}
	w->count += u->count;
}

/*Right tree rotation at node u:
 *                  p                         p      
 *                 /                         /        
 *                u                         w         
 *               / \                       / \
 *              w   c        ---->        a   u       
 *             / \                           / \
 *            a   b                         b   c     
 *                                                     
 */
void stTreap_rotateRight(stTreap *u) {
	stTreap *w = u->left;
	assert(w);
	w->parent = u->parent;
	if(w->parent) {
		if(w->parent->left == u) {
			w->parent->left = w;
		}
		else {
			w->parent->right = w;
		}
	}
	u->left = w->right;
	if(u->left) {
		u->left->parent = u;
	}
	u->parent = w;
	w->right = u;

	u->count -= w->count;
	if(u->left) {
		u->count += u->left->count;
		w->count -= u->left->count;
	}
	w->count += u->count;

}
stTreap *stTreap_findRoot(stTreap *node) {
	while(node->parent) {
		node = node->parent;
	}
	return(node);
}

int stTreap_size(stTreap *node) {
	node = stTreap_findRoot(node);
	return(node->count);
}

/*Finds the closest node to the given key using a binary search. */
stTreap *stTreap_binarySearch(long key, stTreap *node) {
	node = stTreap_findRoot(node);
	while(1) {
		if(key == node->key) {
			break;
		}
		else if(key < node->key) {
			if(node->left) {
				node = node->left;
			}
			else {
				break;
			}
		}
		else {
			if(node->right) {
				node = node->right;
			}
			else {
				break;
			}
		}
	}
	return(node);
}

/* Returns true of the given tree contains the given node. */
int stTreap_contains(stTreap *tree, stTreap *node) {
	stTreap *closestNode = stTreap_binarySearch(node->key, tree);
	return(closestNode->key == node->key);
}
int stTreap_depth(stTreap *node) {
	int depth = 0;
	while(node->parent) {
		node = node->parent;
		depth++;
	}
	return(depth);
}

/*compares two nodes without using their keys. O(lgn) */
int stTreap_compare(stTreap* a, stTreap *b) {
	if (a == b) {
		return 0;
	}
	assert(stTreap_findRoot(a) == stTreap_findRoot(b));
	int depthA = stTreap_depth(a);
	int depthB = stTreap_depth(b);

	while (depthA > depthB) {
		if (a->parent == b) {
			if (a == b->left) {
				return -1;
			} else {
				assert(a == b->right);
				return 1;
			}
		}
		a = a->parent;
		depthA--;

	}
	while (depthB > depthA) {
		if (b->parent == a) {
			if (b == a->left) {
				return 1;
			} else {
				assert(b == a->right);
				return -1;
			}
		}
		b = b->parent;
		depthB--;
	}
	assert(depthA == depthB);
	int depth = depthA;
	while (depth) {
		assert(a != b);
		if (a->parent == b->parent) {
			if (a == a->parent->left) {
				assert(b == a->parent->right);
				return -1;
			} else {
				assert(b == a->parent->left);
				assert(a == a->parent->right);
				return 1;
			}
		}
		a = a->parent;
		b = b->parent;
		depth--;
	}
	assert(a && b);
	assert(!a->parent);
	assert(!b->parent);
	return 0;
}


/*Restores the heap property by using tree 
 * rotations to move a new node up until
 * its priority is less than that of its parent.*/
void stTreap_moveUp(stTreap *node) {
	while(node->parent && (node->priority >= node->parent->priority)) {
		if(node->parent->left == node) {
			stTreap_rotateRight(node->parent);

		}
		else {
			stTreap_rotateLeft(node->parent);
		}
	}
}


/*Insert new key-value pair into the treap at the position determined by
 * a binary search for the key. Then restore the heap property. */
stTreap *stTreap_insert(long key, void *value, stTreap *node) {
	node = stTreap_findRoot(node);	
	stTreap *newNode = stTreap_construct(value);
	newNode->key = key;

	stTreap *parent = stTreap_binarySearch(key, node);
	if(parent->key == key) {
		//duplicate
		free(newNode);
		return(parent);
	}
	newNode->parent = parent;
	if(key < parent->key) {
		parent->left = newNode;
	}
	else {
		parent->right = newNode;
	}

	//update counts for all parents of this node
	while(parent) {
		parent->count += 1;
		parent = parent->parent;
	}

	stTreap_moveUp(newNode);
	return(newNode);
}

/*Perform rotations to make a node into a leaf, which
 * allows it to be deleted without violating the heap property. */
void stTreap_moveDown(stTreap *node) {
	while(node->left || node->right) {
		if(!node->left) {
			stTreap_rotateLeft(node);
		}
		else if(!node->right) {
			stTreap_rotateRight(node);
		}
		else if(node->left->priority < node->right->priority) {
			stTreap_rotateRight(node);
		}
		else {
			stTreap_rotateLeft(node);
		}
	}
}
int stTreap_remove(long key, stTreap *tree) {
	tree = stTreap_findRoot(tree);
	stTreap *rmnode = stTreap_binarySearch(key, tree);
	if(rmnode->key != key) {
		//key not in the treap
		return(0);
	}
	stTreap_moveDown(rmnode);
	if(rmnode->parent) {
		if(rmnode->parent->left == rmnode) {
			rmnode->parent->left = NULL;
		}
		else {
			rmnode->parent->right = NULL;
		}
	}
	free(rmnode);
	//update parent counts
	stTreap *p = rmnode->parent;
	while(p) {
		p->count -= 1;
		p = p->parent;
	}
	//successfully deleted the node
	return(1);
}

/*The node with the minimum key is stored as the leftmost leaf */
stTreap *stTreap_findMin(stTreap *node) {
	assert(node);
	while(node->left) {
		node = node->left;
	}
	return(node);
}
stTreap *stTreap_findMax(stTreap *node) {
	while(node->right) {
		node = node->right;
	}
	return(node);
}
/*splits the treap containing a node into two treaps, one whose
 * keys are all less than that of the node, and one whose keys are all
 * greater. Leaves the left subtree connected to the split point node,
 * but disconnects the right subtree. Returns the right subtree.*/
stTreap *stTreap_splitAfter(stTreap *node) {
	node->priority = INT_MAX;
	stTreap_moveUp(node);
    assert(stTreap_findRoot(node) == node);
	assert(node->parent == NULL); //node should now be the new root of the treap
	stTreap *rightSubtree = node->right;
	if(rightSubtree) {
		node->count -= rightSubtree->count;
		rightSubtree->parent = NULL;
		node->right = NULL;
	}
	stTreap_chooseNewPriority(node);
	return(rightSubtree);
}

/*same as splitAfter, but returns the left subtree.*/
stTreap *stTreap_splitBefore(stTreap *node) {
	node->priority = INT_MAX;
	stTreap_moveUp(node);
	assert(node);
	assert(stTreap_findRoot(node) == node);
	assert(node->parent == NULL);
	stTreap *leftSubtree = node->left;
	if(leftSubtree) {
		node->count -= leftSubtree->count;
		leftSubtree->parent = NULL;
		node->left = NULL;
	}
	stTreap_chooseNewPriority(node);
	return(leftSubtree);
}
//choose a priority value for a node after moving it to the root
void stTreap_chooseNewPriority(stTreap *node) {
	int pleft = 0;
	int pright = 0;
	if (node->left) {
		pleft = node->left->priority;
	}
	if (node->right) {
		pright = node->right->priority;
	}
	if (pleft > pright) {
		node->priority = pleft + 1;
	}
	else {
		node->priority = pright + 1;
	}
}

/*returns the next node in an in-order traversal of the treap
 * starting at node. */
stTreap *stTreap_next(stTreap *node) {
	assert(node);
	if(node->right != NULL) {
		return(stTreap_findMin(node->right));
	}
	stTreap *p = node->parent;
	while(p != NULL && (node == p->right)) {
		node = p;
		p = p->parent;
	}
	return(p);
}
stTreap *stTreap_prev(stTreap *node) {
	assert(node);
	if(node->left) {
		return(stTreap_findMax(node->left));
	}
	stTreap *p = node->parent;
	while(p && (node == p->left)) {
		node = p;
		p = p->parent;
	}
	return(p);
}


stTreap *stTreap_concat(stTreap *a, stTreap *b) {
	if(!a || !b) {
		return(NULL);
	}
	stTreap *ra = stTreap_findRoot(a);
	//stTreap *ta = ra;

	//move to largest element in a
	//while(ta->right) {
	//	ta = ta->right;
	//}
	stTreap *rb = stTreap_findRoot(b);

	//a and b shouldn't already be concatenated
	assert(rb != ra);
	//stTreap *sb = rb;

	//move to smallest element in b
	//while(sb->left) {
	//	sb = sb->left;
	//}
	stTreap *r = stTreap_concatRecurse(ra, rb);
	r->parent = NULL;
	return(r);
}

/* Joins two treaps, preserving the heap property for priorities and the key orderings, 
 * assuming that all elements in the left tree are smaller than all elements in the right tree.
 * This will be true if the trees were produced by splitting a treap.
 *
 *         d10            x8    
 *         /  \          /  \
 *        b7   e5   +   w6   z4   (Here, keys are letters and priorities
 *       /  \                /        are numbers)
 *      a3  c2            y1
 * 
 * concatRecurse(d10, x8)
 *
 *                 d10
 *                /    \
 *               /      \
 *              b7       x8        In-order traversal: abcdewxyz
 *             /  \     /  \
 *            a3   c2  w6  z4
 *                    /    /
 *                   e5   y1
 *
 *	concatRecurse(e5, x8)
 *                    x8
 *                   /  \
 *                  w6   z4
 *                 /    /
 *                e5   y1
 *
 *		concatRecurse(e5, w6)
 *			w6
 *                     /
 *                    e5
 *			concatRecurse(e5, NULL)
 *				
 *				e5
 *
 */
stTreap *stTreap_concatRecurse(stTreap *a, stTreap *b) {
	if (!a) {
		return(b);
	}
	else if (!b) {
		return(a);
	}
	else if (a->priority > b->priority) {
		if(a->right) {
			a->count -= a->right->count;
		}
		a->right = stTreap_concatRecurse(a->right, b);
		a->count += a->right->count;
		a->right->parent = a;
		return(a);
	}
	else {
		if(b->left) {
			b->count -= b->left->count;
		}

		b->left = stTreap_concatRecurse(a, b->left);
		b->count += b->left->count;
		b->left->parent = b;
		return(b);
	}
}
