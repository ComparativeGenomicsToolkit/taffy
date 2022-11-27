#include "sonLibGlobalsInternal.h"

struct _stUnionFind {
    stHash *objectToEntry;
};

struct _stUnionFindIt {
    stList *sets;
    int64_t curIndex;
};

typedef struct _stUnionFindEntry {
    struct _stUnionFindEntry *parent;
    void *object;
    int64_t rank;
} stUnionFindEntry;

static void stUnionFindEntry_destruct(stUnionFindEntry *entry) {
    free(entry);
}

stUnionFind *stUnionFind_construct(void) {
    stUnionFind *ret = st_malloc(sizeof(stUnionFind));
    ret->objectToEntry = stHash_construct2(NULL, (void (*)(void *)) stUnionFindEntry_destruct);
    return ret;
}

void stUnionFind_destruct(stUnionFind *unionFind) {
    stHash_destruct(unionFind->objectToEntry);
    free(unionFind);
}

void stUnionFind_add(stUnionFind *unionFind, void *object) {
    stUnionFindEntry *entry = st_malloc(sizeof(stUnionFindEntry));
    entry->object = object;
    entry->parent = NULL;
    entry->rank = 0;
    stHash_insert(unionFind->objectToEntry, object, entry);
}

void stUnionFind_union(stUnionFind *unionFind, void *object1, void *object2) {
    stUnionFindEntry *root1 = stUnionFind_find(unionFind, object1);
    stUnionFindEntry *root2 = stUnionFind_find(unionFind, object2);

    assert(root1->parent == NULL && root2->parent == NULL);

    if (root1 == root2) {
        return;
    }
    // keep the tree relatively balanced by checking rank
    if (root1->rank > root2->rank) {
        root2->parent = root1;
    } else if (root1->rank < root2->rank) {
        root1->parent = root2;
    } else {
        root1->parent = root2;
        root2->rank++;
    }
}

static stUnionFindEntry *find(stUnionFindEntry *entry) {
    // Recursively compress the path to the root by making every node
    // visited point directly to the root.
    if (entry->parent != NULL) {
        entry->parent = find(entry->parent);
    }
    return entry->parent == NULL ? entry : entry->parent;
}

void *stUnionFind_find(stUnionFind *unionFind, void *object) {
    stUnionFindEntry *entry = stHash_search(unionFind->objectToEntry, object);
    assert(entry != NULL);
    return find(entry);
}

stUnionFindIt *stUnionFind_getIterator(stUnionFind *unionFind) {
    stHash *rootToSet = stHash_construct();
    stList *entries = stHash_getValues(unionFind->objectToEntry);
    for (int64_t i = 0; i < stList_length(entries); i++) {
        stUnionFindEntry *entry = stList_get(entries, i);
        stUnionFindEntry *root = find(entry);
        stSet *set = stHash_search(rootToSet, root);
        if (set == NULL) {
            set = stSet_construct();
            stHash_insert(rootToSet, root, set);
        }
        stSet_insert(set, entry->object);
    }
    stUnionFindIt *it = st_malloc(sizeof(stUnionFindIt));
    it->sets = stHash_getValues(rootToSet);
    stList_setDestructor(it->sets, (void (*)(void *)) stSet_destruct);
    it->curIndex = 0;

    stList_destruct(entries);
    stHash_destruct(rootToSet);
    return it;
}

stSet *stUnionFindIt_getNext(stUnionFindIt *it) {
    if (it->curIndex < stList_length(it->sets)) {
        return stList_get(it->sets, it->curIndex++);
    } else {
        return NULL;
    }
}

void stUnionFind_destructIterator(stUnionFindIt *it) {
    stList_destruct(it->sets);
    free(it);
}
