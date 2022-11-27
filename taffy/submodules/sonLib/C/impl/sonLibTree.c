/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * tree.c
 *
 *  Created on: 21 May 2010
 *      Author: benedictpaten
 */

#include "sonLibGlobalsInternal.h"
#include <stdbool.h>
#include "sonLibString.h"

/*
 * The actual datastructure.
 */

struct _stTree {
    double branchLength;
    stList *nodes;
    char *label;
    void *clientData;
    stTree *parent;
};

/*
 * The functions..
 */

stTree *stTree_construct(void) {
    stTree *tree = st_malloc(sizeof(stTree));
    tree->branchLength = INFINITY;
    tree->nodes = stList_construct3(0, (void (*)(void *))stTree_destruct);
    tree->label = NULL;
    tree->parent = NULL;
    tree->clientData = NULL;
    return tree;
}

void stTree_destruct(stTree *tree) {
    stList_destruct(tree->nodes);
    if(tree->label != NULL) {
        free(tree->label);
    }
    free(tree);
}

/* clone a node */
stTree *stTree_cloneNode(stTree *node) {
    stTree *node2 = stTree_construct();
    stTree_setBranchLength(node2, stTree_getBranchLength(node));
    stTree_setClientData(node2, stTree_getClientData(node));
    stTree_setLabel(node2, stTree_getLabel(node));
    return node2;
}

/* recursively clone a tree */
static stTree *tree_clonetree(stTree *node, stTree *parent2) {
    stTree *node2 = stTree_cloneNode(node);
    stTree_setParent(node2, parent2);
    for (int i = 0; i < stTree_getChildNumber(node); i++) {
        tree_clonetree(stTree_getChild(node, i), node2);
    }
    return node2;
}

stTree *stTree_clone(stTree *root) {
    return tree_clonetree(root, NULL);
}

// Set client data to NULL (optionally recursively).
static void stTree_clearClientData(stTree *tree, bool recursive) {
    stTree_setClientData(tree, NULL);
    if (recursive) {
        for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
            stTree_clearClientData(stTree_getChild(tree, i), true);
        }
    }
}

// clone this node, make its supertree a child, and clone all children
// other than oldNode, leaving this node as a child of nodeToAddTo
static void tree_cloneFlippedTree(stTree *node, stTree *oldNode,
                                  stTree *nodeToAddTo,
                                  double branchLength) {
    if(stTree_getParent(node) != NULL || stTree_getChildNumber(node) > 2) {
        // This node isn't the root
        stTree *clonedNode = stTree_cloneNode(node);
        stTree_setParent(clonedNode, nodeToAddTo);
        stTree_setBranchLength(clonedNode, branchLength);
        // Clone its children (other than oldNode) and their subtrees
        for(int64_t i = 0; i < stTree_getChildNumber(node); i++) {
            stTree *child = stTree_getChild(node, i);
            if(child != oldNode) {
                stTree *clonedChild = stTree_clone(child);
                stTree_setParent(clonedChild, clonedNode);
            }
        }

        if (stTree_getParent(node) != NULL) {
            // Recurse on the parent of this node.
            tree_cloneFlippedTree(stTree_getParent(node), node,
                                  clonedNode, stTree_getBranchLength(node));
        }
    } else if (stTree_getParent(node) == NULL) {
        // We have to treat the root specially in binary trees,
        // because we're going to eliminate it. Just add all the other
        // children of the root as children of nodeToAddTo.
        for(int64_t i = 0; i < stTree_getChildNumber(node); i++) {
            stTree *child = stTree_getChild(node, i);
            if(child != oldNode) {
                stTree *clonedChild = stTree_clone(child);
                stTree_setParent(clonedChild, nodeToAddTo);
                stTree_setBranchLength(clonedChild, stTree_getBranchLength(child) + branchLength);
            }
        }
    }
}

// Return a new tree rooted a given distance above the given node.
static stTree *stTree_reRootP(stTree *node, double distanceAbove, bool clearClientData) {
    if(stTree_getParent(node) == NULL) {
        // This node is already the root.
        stTree *newRoot = stTree_clone(node);
        if (clearClientData) {
            stTree_clearClientData(newRoot, true);
        }
        return newRoot;
    }

    assert(stTree_getBranchLength(node) >= distanceAbove);

    stTree *newRoot = stTree_construct();
    // This node and its children (if any) are fine already.
    stTree *clonedNode = stTree_clone(node);
    stTree_setParent(clonedNode, newRoot);
    stTree_setBranchLength(clonedNode, distanceAbove);
    tree_cloneFlippedTree(stTree_getParent(node), node, newRoot,
                          stTree_getBranchLength(node) - distanceAbove);
    // Having the same client data can be a problem
    if (clearClientData) {
        stTree_clearClientData(newRoot, true);
    }
    return newRoot;
}

stTree *stTree_reRoot(stTree *node, double distanceAbove) {
    return stTree_reRootP(node, distanceAbove, true);
}

stTree *stTree_reRootAndKeepClientData(stTree *node, double distanceAbove) {
    return stTree_reRootP(node, distanceAbove, false);
}

stTree *stTree_getParent(stTree *tree) {
    return tree->parent;
}

void stTree_setParent(stTree *tree, stTree *parent) {
    if(stTree_getParent(tree) != NULL) {
        stList_removeItem(stTree_getParent(tree)->nodes, tree);
    }
    tree->parent = parent;
    if(parent != NULL) {
        stList_append(parent->nodes, tree);
    }
}

int64_t stTree_getChildNumber(stTree *tree) {
    return stList_length(tree->nodes);
}

stTree *stTree_getChild(stTree *tree, int64_t i) {
    return stList_get(tree->nodes, i);
}

stTree *stTree_findChild(stTree *tree, const char *label) {
    if (tree->label != NULL && strcmp(tree->label, label) == 0) {
        return tree;
    }
    for (int i = 0; i < stList_length(tree->nodes); i++) {
        stTree *node = (stTree *)stList_get(tree->nodes, i);
        if ((node->label != NULL) && (strcmp(node->label, label) == 0)) {
            return node;
        }
        stTree *hit = stTree_findChild(node, label);
        if (hit != NULL) {
            return hit;
        }
    }
    return NULL;
}

double stTree_getBranchLength(stTree *tree) {
    return tree->branchLength;
}

void stTree_setBranchLength(stTree *tree, double distance) {
    tree->branchLength = distance;
}

void *stTree_getClientData(stTree *tree) {
    return tree->clientData;
}

void stTree_setClientData(stTree *tree, void *clientData) {
    tree->clientData = clientData;
}


const char *stTree_getLabel(stTree *tree) {
    return tree->label;
}

void stTree_setLabel(stTree *tree, const char *label) {
    if(tree->label != NULL) {
        free(tree->label);
    }
    tree->label = label == NULL ? NULL : stString_copy(label);
}

int stTree_getNumNodes(stTree *root) {
    int cnt = 1; // this node
    for (int i = 0; i < stTree_getChildNumber(root); i++) {
        cnt += stTree_getNumNodes(stTree_getChild(root, i));
    }
    return cnt;
}

/////////////////////////////
//Newick tree parser
/////////////////////////////

/*
 * Gets the next token from the list.
 */
static void tree_parseNewickTreeString_getNextToken(char **token, char **newickTreeString) {
    assert(*token != NULL);
    free(*token);
    *token = stString_getNextWord(newickTreeString);
    if (*token == NULL) {
        st_errAbort("invalid Newick tree, doesn't end with ';': %s", **newickTreeString);
    }
}

/*
 * Sets the label, if the token is a label and updates the token.
 */
static void tree_parseNewickString_getLabel(char **token, char **newickTreeString, stTree *tree) {
    if(**token != ':' && **token != ',' && **token != ';' && **token != ')') {
        stTree_setLabel(tree, *token);
        tree_parseNewickTreeString_getNextToken(token, newickTreeString);
    }
}

/*
 * Parses any available branch length and updates the token.
 */
static void tree_parseNewickString_getBranchLength(char **token, char **newickTreeString, stTree *tree) {
    if (**token == ':') {
        tree_parseNewickTreeString_getNextToken(token, newickTreeString);
        double distance;
        int j = sscanf(*token, "%lf", &distance);
        if (j != 1) {
            st_errAbort("Invalid Newick tree");
        }
        stTree_setBranchLength(tree, distance);
        tree_parseNewickTreeString_getNextToken(token, newickTreeString);
    }
}

static stTree *tree_parseNewickStringP(char **token, char **newickTreeString) {
    stTree *tree = stTree_construct();
    if((*token)[0] == '(') {
        if (strlen(*token) != 1) {
            st_errAbort("invalid Newick tree");
        }
        tree_parseNewickTreeString_getNextToken(token, newickTreeString);
        while(1) {
            stTree_setParent(tree_parseNewickStringP(token, newickTreeString), tree);
            if (strlen(*token) != 1) {
                st_errAbort("invalid Newick tree");
            }
            if ((*token)[0] == ',') {
                tree_parseNewickTreeString_getNextToken(token, newickTreeString);
            } else {
                break;
            }
        }
        if ((*token)[0] != ')') {
            st_errAbort("Invalid Newick tree, unbalanced parenthesis");
        }
        tree_parseNewickTreeString_getNextToken(token, newickTreeString);
    }
    tree_parseNewickString_getLabel(token, newickTreeString, tree);
    tree_parseNewickString_getBranchLength(token, newickTreeString, tree);
    if (!((**token == ',') || (**token == ';') || (**token == ')'))) {
        st_errAbort("Invalid Newick tree, missing ',', ';', or ')'");
    }
    return tree;
}

stTree *stTree_parseNewickString(const char *string) {
    //lax newick tree parser
    char *cA = stString_replace(string, "(", " ( ");
    char *cA2 = stString_replace(cA, ")", " ) ");
    free(cA);
    cA = cA2;
    cA2 = stString_replace(cA, ":", " : ");
    free(cA);
    cA = cA2;
    cA2 = stString_replace(cA, ",", " , ");
    free(cA);
    cA = cA2;
    cA2 = stString_replace(cA, ";", " ; ");
    free(cA);
    cA = cA2;
    char *token = stString_getNextWord(&cA);
    if (token == NULL) {
        st_errAbort("Invalid Newick tree");
    }
    stTree *tree = tree_parseNewickStringP(&token, &cA);
    if (*token != ';') {
        st_errAbort("Invalid Newick tree, missing ';'");
    }
    free(cA2);
    free(token);
    return tree;
}

/////////////////////////////
//Newick tree writer
/////////////////////////////

static char *tree_getNewickTreeStringP(stTree *tree) {
    char *cA, *cA2;
    if(stTree_getChildNumber(tree) > 0) {
        int64_t i;
        cA = stString_copy("(");
        for(i=0; i<stTree_getChildNumber(tree); i++) {
            cA2 = tree_getNewickTreeStringP(stTree_getChild(tree, i));
            char *cA3 = stString_print((i+1 < stTree_getChildNumber(tree) ? "%s%s," : "%s%s"), cA, cA2);
            free(cA);
            free(cA2);
            cA = cA3;
        }
        cA2 = stString_print("%s)", cA);
        free(cA);
        cA = cA2;
    }
    else {
        cA = stString_copy("");
    }
    if(stTree_getLabel(tree) != NULL) {
        cA2 = stString_print("%s%s", cA, stTree_getLabel(tree));
        free(cA);
        cA = cA2;
    }
    if(stTree_getBranchLength(tree) != INFINITY) {
        char *cA2 = stString_print("%s:%g", cA, stTree_getBranchLength(tree));
        free(cA);
        cA = cA2;
    }
    return cA;
}

char *stTree_getNewickTreeString(stTree *tree) {
    char *cA = tree_getNewickTreeStringP(tree), *cA2;
    cA2 = stString_print("%s;", cA);
    free(cA);
    return cA2;
}

bool stTree_equals(stTree *tree1, stTree *tree2) {
    if (stTree_getBranchLength(tree1) != stTree_getBranchLength(tree2)) {
        return false;
    }
    if (!stString_eq(stTree_getLabel(tree1), stTree_getLabel(tree2))) {
        return false;
    }
    int numChildren = stTree_getChildNumber(tree1);
    if (stTree_getChildNumber(tree2) != numChildren) {
        return false;
    }
    for (int i = 0; i < numChildren; i++) {
        if (!stTree_equals(stTree_getChild(tree1, i), stTree_getChild(tree2, i))) {
            return false;
        }
    }
    return true;
}

/* must warp function pointer in data, as ISO C doesn't allow conversions
 * between void* and function pointers */
struct sortFuncArgs {
    int (*cmpFn)(stTree *a, stTree *b);
};

static int sortChildrenListCmpFn(const void *a, const void *b, void* args) {
    return ((struct sortFuncArgs*)args)->cmpFn((stTree*)a, (stTree*)b);
}

void stTree_sortChildren(stTree *root, int (*cmpFn)(stTree *a, stTree *b)) {
    struct sortFuncArgs args = {cmpFn};
    stList_sort2(root->nodes, sortChildrenListCmpFn, &args);
    for (int i = 0; i < stTree_getChildNumber(root); i++) {
        stTree_sortChildren(stTree_getChild(root, i), cmpFn);
    }
}

void stTree_setChild(stTree *tree, int64_t i, stTree *child) {
    assert(i >= 0);
    assert(i < stTree_getChildNumber(tree));
    stList_set(tree->nodes, i, child);
}

stTree *stTree_getMRCA(stTree *node1, stTree *node2) {
    // Find all of node 1's parents (inclusive of node 1)
    stSet *parents = stSet_construct();
    stTree *curNode = node1;
    do {
        stSet_insert(parents, curNode);
    } while ((curNode = stTree_getParent(curNode)) != NULL);

    // Find the first parent of node 2 that is a parent of node 1
    stTree *ret = NULL;
    curNode = node2;
    do {
        if (stSet_search(parents, curNode) != NULL) {
            ret = curNode;
            break;
        }
    } while ((curNode = stTree_getParent(curNode)) != NULL);

    stSet_destruct(parents);
    return ret;
}
