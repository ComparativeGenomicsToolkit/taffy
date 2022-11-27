#ifndef TREAP_H
#define TREAP_H

struct _stTreap {
	//heap property is maintained for the priority, so highest priority is at the top
	//of the tree.
	int priority;

	//number of nodes below and including this one
	int count;

	//ordering of the keys is maintained from left to right. Keys don't have to be used
	//if all insertions are done through merges, and ordering will still be maintained
	//by the other operations
	long key;
	void *value;

	stTreap *parent;
	stTreap *left;
	stTreap *right;
};

stTreap *stTreap_construct(void *value);
void *stTreap_getValue(stTreap *node);
void stTreap_destruct(stTreap *node);
void stTreap_destructRecurse(stTreap *root);
void stTreap_nodeDestruct(stTreap *node);
void stTreap_rotateLeft(stTreap *u);
void stTreap_rotateRight(stTreap *u);
stTreap *stTreap_findRoot(stTreap *node);
char *stTreap_print(stTreap *node);
char *stTreap_printBackwards(stTreap *node);
int stTreap_size(stTreap *node);
stTreap *stTreap_binarySearch(long key, stTreap *root);
int stTreap_contains(stTreap *tree, stTreap *node);
int stTreap_depth(stTreap *node);
int stTreap_compare(stTreap* a, stTreap *b);
void stTreap_moveUp(stTreap *node);
stTreap *stTreap_insert(long key, void *value, stTreap *root);
void stTreap_moveDown(stTreap *node);
int stTreap_remove(long key, stTreap *root);
stTreap *stTreap_findMin(stTreap *node);
stTreap *stTreap_findMax(stTreap *node);
stTreap *stTreap_splitAfter(stTreap *node);
stTreap *stTreap_splitBefore(stTreap *node);
void stTreap_chooseNewPriority(stTreap *node);
stTreap *stTreap_next(stTreap *node);
stTreap *stTreap_prev(stTreap *node);
stTreap *stTreap_concat(stTreap *a, stTreap *b);
stTreap *stTreap_concatRecurse(stTreap *a, stTreap *b);

#endif
