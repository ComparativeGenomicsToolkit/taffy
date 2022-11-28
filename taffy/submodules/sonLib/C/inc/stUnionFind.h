// An implementation of a fast union-find data structure with
// link-by-rank and path compression. This is equivalent to a
// partially-dynamic version of stConnectivity.

typedef struct _stUnionFind stUnionFind;

typedef struct _stUnionFindIt stUnionFindIt;

// Create an empty union-find structure.
stUnionFind *stUnionFind_construct(void);

// Free the union-find structure.
void stUnionFind_destruct(stUnionFind *unionFind);

// Add a new isolated component containing only "object".
void stUnionFind_add(stUnionFind *unionFind, void *object);

// Merge two components.
void stUnionFind_union(stUnionFind *unionFind, void *object1, void *object2);

// Find the component for this object.
void *stUnionFind_find(stUnionFind *unionFind, void *object);

// Get an iterator over the components of the union-find.
stUnionFindIt *stUnionFind_getIterator(stUnionFind *unionFind);

// Get a component, as a set of objects, from the iterator.
// Returns NULL at the end of iteration.
stSet *stUnionFindIt_getNext(stUnionFindIt *it);

// Free the iterator.
void stUnionFind_destructIterator(stUnionFindIt *it);
