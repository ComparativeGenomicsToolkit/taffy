#include <stdlib.h>
#include "sonLib.h"
#include "stPhylogeny.h"
// QuickTree includes
#include "cluster.h"
#include "tree.h"
#include "buildtree.h"
// Spimap C/C++ translation layer include
#include "stSpimapLayer.h"

// Free stReconciliationInfo properly.
static void stReconciliationInfo_destruct(stReconciliationInfo *info) {
    free(info);
}

static stReconciliationInfo *stReconciliationInfo_clone(stReconciliationInfo *info) {
    stReconciliationInfo *ret = st_malloc(sizeof(stReconciliationInfo));
    ret->species = info->species;
    ret->event = info->event;
    return ret;
}

static void stIndexedTreeInfo_destruct(stIndexedTreeInfo *info) {
    assert(info != NULL);
    free(info->leavesBelow);
    free(info);
}

// Free a stPhylogenyInfo struct
void stPhylogenyInfo_destruct(stPhylogenyInfo *info) {
    if (info->index != NULL) {
        stIndexedTreeInfo_destruct(info->index);
    }
    if (info->recon != NULL) {
        stReconciliationInfo_destruct(info->recon);
    }
    free(info);
}

// Clone a stIndexedTreeInfo struct
static stIndexedTreeInfo *stIndexedTreeInfo_clone(stIndexedTreeInfo *info) {
    stIndexedTreeInfo *ret = st_malloc(sizeof(stIndexedTreeInfo));
    memcpy(ret, info, sizeof(stIndexedTreeInfo));
    ret->leavesBelow = malloc(ret->totalNumLeaves * sizeof(char));
    memcpy(ret->leavesBelow, info->leavesBelow, ret->totalNumLeaves * sizeof(char));
    return ret;
}

// Clone a stPhylogenyInfo struct
stPhylogenyInfo *stPhylogenyInfo_clone(stPhylogenyInfo *info) {
    stPhylogenyInfo *ret = st_malloc(sizeof(stPhylogenyInfo));
    if (info->index != NULL) {
        ret->index = stIndexedTreeInfo_clone(info->index);
    }
    if (info->recon != NULL) {
        ret->recon = stReconciliationInfo_clone(info->recon);
    } else {
        ret->recon = NULL;
    }
    return ret;
}

// Free the stPhylogenyInfo struct for this node and all nodes below it.
void stPhylogenyInfo_destructOnTree(stTree *tree) {
    stPhylogenyInfo_destruct(stTree_getClientData(tree));
    stTree_setClientData(tree, NULL);
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        stPhylogenyInfo_destructOnTree(stTree_getChild(tree, i));
    }
}

// Helper function to add the stIndexedTreeInfo that is normally
// generated during neighbor-joining to a tree that has leaf-labels 0,
// 1, 2, etc.
void addStIndexedTreeInfoR(stTree *tree)
{
    stPhylogenyInfo *info = stTree_getClientData(tree);
    if (info == NULL) {
        info = st_calloc(1, sizeof(stPhylogenyInfo));
        stTree_setClientData(tree, info);
    }
    if (info->index != NULL) {
        stIndexedTreeInfo_destruct(info->index);
    }
    stIndexedTreeInfo *indexInfo = st_calloc(1, sizeof(stIndexedTreeInfo));
    info->index = indexInfo;
    if(stTree_getChildNumber(tree) == 0) {
        int ret;
        ret = sscanf(stTree_getLabel(tree), "%" PRIi64, &indexInfo->matrixIndex);
        (void) ret;
        assert(ret == 1);
    } else {
        indexInfo->matrixIndex = -1;
        for(int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
            addStIndexedTreeInfoR(stTree_getChild(tree, i));
        }
    }
}

int64_t stTree_getNumLeaves(stTree *tree) {
    int64_t numLeaves = 0;
    stList *stack = stList_construct();
    stList_append(stack, tree);
    while (stList_length(stack) != 0) {
        tree = stList_pop(stack);
        for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
            stList_append(stack, stTree_getChild(tree, i));
        }
        if (stTree_getChildNumber(tree) == 0) {
            numLeaves++;
        }
    }
    stList_destruct(stack);
    return numLeaves;
}

void stPhylogeny_addStIndexedTreeInfo(stTree *tree) {
    addStIndexedTreeInfoR(tree);
    stPhylogeny_setLeavesBelow(tree, stTree_getNumLeaves(tree));
}


// Set (and allocate) the leavesBelow and totalNumLeaves attribute in
// the phylogenyInfo for the given tree and all subtrees. The
// phylogenyInfo structure (in the clientData field) must already be
// allocated!
void stPhylogeny_setLeavesBelow(stTree *tree, int64_t totalNumLeaves)
{
    int64_t i, j;
    assert(stTree_getClientData(tree) != NULL);
    stPhylogenyInfo *info = stTree_getClientData(tree);
    assert(info->index != NULL);
    stIndexedTreeInfo *indexInfo = info->index;
    for (i = 0; i < stTree_getChildNumber(tree); i++) {
        stPhylogeny_setLeavesBelow(stTree_getChild(tree, i), totalNumLeaves);
    }

    indexInfo->totalNumLeaves = totalNumLeaves;
    if (indexInfo->leavesBelow != NULL) {
        // leavesBelow has already been allocated somewhere else, free it.
        free(indexInfo->leavesBelow);
    }
    indexInfo->leavesBelow = st_calloc(totalNumLeaves, sizeof(char));
    if (stTree_getChildNumber(tree) == 0) {
        assert(indexInfo->matrixIndex < totalNumLeaves);
        assert(indexInfo->matrixIndex >= 0);
        indexInfo->leavesBelow[indexInfo->matrixIndex] = 1;
    } else {
        for (i = 0; i < totalNumLeaves; i++) {
            for (j = 0; j < stTree_getChildNumber(tree); j++) {
                stPhylogenyInfo *childInfo = stTree_getClientData(stTree_getChild(tree, j));
                stIndexedTreeInfo *childIndexInfo = childInfo->index;
                indexInfo->leavesBelow[i] |= childIndexInfo->leavesBelow[i];
            }
        }
    }
}

static stTree *quickTreeToStTreeR(struct Tnode *tNode) {
    stTree *ret = stTree_construct();
    bool hasChild = false;
    if (tNode->left != NULL) {
        stTree *left = quickTreeToStTreeR(tNode->left);
        stTree_setParent(left, ret);
        hasChild = true;
    }
    if (tNode->right != NULL) {
        stTree *right = quickTreeToStTreeR(tNode->right);
        stTree_setParent(right, ret);
        hasChild = true;
    }

    // Allocate the phylogenyInfo for this node.
    stPhylogenyInfo *info = st_calloc(1, sizeof(stPhylogenyInfo));
    stIndexedTreeInfo *indexInfo = st_calloc(1, sizeof(stIndexedTreeInfo));
    info->index = indexInfo;
    if (!hasChild) {
        indexInfo->matrixIndex = tNode->nodenumber;
    } else {
        indexInfo->matrixIndex = -1;
    }
    stTree_setClientData(ret, info);

    stTree_setBranchLength(ret, tNode->distance);

    // Can remove if needed, probably not useful except for testing.
    char *label = stString_print_r("%u", tNode->nodenumber);
    stTree_setLabel(ret, label);
    free(label);

    return ret;
}

// Helper function for converting an unrooted QuickTree Tree into an
// stTree. The tree is rooted halfway along the longest branch if
// outgroups is NULL, otherwise it's rooted halfway along the longest
// branch to an outgroup.
static stTree *quickTreeToStTree(struct Tree *tree, stList *outgroups) {
    struct Tree *rootedTree = get_root_Tnode(tree);
    stTree *ret = quickTreeToStTreeR(rootedTree->child[0]);
    stPhylogeny_setLeavesBelow(ret, (stTree_getNumNodes(ret) + 1) / 2);
    if(outgroups != NULL && stList_length(outgroups) != 0) {
        // Find the longest branch to an outgroup and root the tree
        // there. If we're desperate for speedups, we can just root
        // the tree once instead of rooting once and then re-rooting.
        assert(stList_length(outgroups) != 0);
        double maxLength = -1;
        stTree *maxNode = NULL;
        for(int64_t i = 0; i < stList_length(outgroups); i++) {
            int64_t outgroupIndex = stIntTuple_get(stList_get(outgroups, i), 0);
            stTree *outgroup = stPhylogeny_getLeafByIndex(ret, outgroupIndex);
            if(stTree_getBranchLength(outgroup) > maxLength) {
                maxLength = stTree_getBranchLength(outgroup);
                maxNode = outgroup;
            }
        }
        assert(maxNode != NULL);
        stTree *reRooted = stTree_reRoot(maxNode, maxLength/2);

        // Get rid of the old tree
        stPhylogenyInfo_destructOnTree(ret);
        stTree_destruct(ret);

        stPhylogeny_addStIndexedTreeInfo(reRooted);
        ret = reRooted;
    }
    free_Tree(tree);
    free_Tree(rootedTree);
    return ret;
}

// Compare a single partition to a single bootstrap partition and
// update its support if they are identical.
static void updatePartitionSupportFromPartition(stTree *partitionToScore,
                                                stTree *originalPartition,
                                                stTree *bootstrap) {
    stPhylogenyInfo *partitionInfo, *bootstrapInfo;

    partitionInfo = stTree_getClientData(partitionToScore);
    bootstrapInfo = stTree_getClientData(bootstrap);
    assert(partitionInfo != NULL && partitionInfo->index != NULL);
    assert(partitionInfo != NULL && bootstrapInfo->index != NULL);
    assert(partitionInfo->index->totalNumLeaves == bootstrapInfo->index->totalNumLeaves);
    // Check if the set of leaves is equal in both partitions. If not,
    // the partitions can't be equal.
    if (memcmp(partitionInfo->index->leavesBelow, bootstrapInfo->index->leavesBelow,
            partitionInfo->index->totalNumLeaves * sizeof(char))) {
        return;
    }
    // The partitions are equal, increase the support 
    partitionInfo->index->numBootstraps++;
}

// Increase a partition's bootstrap support by supplying the given
// function with the closest possible bootstrap candidate from this
// tree (i.e. it may not be exactly the same branch, if that branch is
// not in the bootstrap).
static void updateSupportFromTree(stTree *partitionToScore,
                                  stTree *originalPartition,
                                  stTree *bootstrapTree,
                                  void (*updateAgainstBootstrapCandidate)(stTree *, stTree *, stTree *))
{
    stPhylogenyInfo *partitionInfo = stTree_getClientData(partitionToScore);
    stPhylogenyInfo *bootstrapInfo = stTree_getClientData(bootstrapTree);
    assert(partitionInfo != NULL);
    assert(bootstrapInfo != NULL);
    stIndexedTreeInfo *partitionIndex = partitionInfo->index;
    stIndexedTreeInfo *bootstrapIndex = bootstrapInfo->index;
    assert(partitionIndex != NULL);
    assert(bootstrapIndex != NULL);
    // This partition should be updated against the bootstrap
    // partition only if none of the bootstrap's children have a leaf
    // set that is a superset of the partition's leaf set.
    bool checkThisPartition = TRUE;

    (void) bootstrapIndex;
    assert(partitionIndex->totalNumLeaves == bootstrapIndex->totalNumLeaves);
    // Check that the leaves under the partition are a subset of the
    // leaves under the current bootstrap node. This should always be
    // true, since that's checked before running this function
    for (int64_t i = 0; i < partitionIndex->totalNumLeaves; i++) {
        if (partitionIndex->leavesBelow[i]) {
            assert(bootstrapIndex->leavesBelow[i]);
        }
    }
    
    for(int64_t i = 0; i < stTree_getChildNumber(bootstrapTree); i++) {
        stTree *bootstrapChild = stTree_getChild(bootstrapTree, i);
        stPhylogenyInfo *bootstrapChildInfo = stTree_getClientData(bootstrapChild);
        assert(bootstrapChildInfo != NULL);
        stIndexedTreeInfo *bootstrapChildIndex = bootstrapChildInfo->index;
        // If any of the bootstrap's children has a leaf set that is a
        // superset of the partition's leaf set, we should update
        // against that child instead.
        bool isSuperset = TRUE;
        for (int64_t j = 0; j < partitionIndex->totalNumLeaves; j++) {
            if (partitionIndex->leavesBelow[j]) {
                if (!bootstrapChildIndex->leavesBelow[j]) {
                    isSuperset = FALSE;
                    break;
                }
            }
        }
        if (isSuperset) {
            updateSupportFromTree(partitionToScore, originalPartition,
                                  bootstrapChild,
                                  updateAgainstBootstrapCandidate);
            checkThisPartition = FALSE;
            break;
        }
    }

    if (checkThisPartition) {
        // This bootstrap partition is the closest candidate. Check
        // the partition against this node.
        updateAgainstBootstrapCandidate(partitionToScore, originalPartition,
                                        bootstrapTree);
    }
}

// Return a new tree which has its partitions scored by how often they
// appear in the bootstrap. This fills in the numBootstraps and
// bootstrapSupport fields of each node. All trees must have valid
// stPhylogenyInfo.
stTree *stPhylogeny_scoreFromBootstrap(stTree *tree, stTree *bootstrap)
{
    stList *list = stList_construct();
    stList_append(list, bootstrap);
    stTree *ret = stPhylogeny_scoreFromBootstraps(tree, list);
    stList_destruct(list);
    return ret;
}

// Return a new tree which has its partitions scored by how often they
// appear in the bootstraps. This fills in the numBootstraps and
// bootstrapSupport fields of each node. All trees must have valid
// stPhylogenyInfo.
stTree *stPhylogeny_scoreFromBootstraps(stTree *tree, stList *bootstraps)
{
    int64_t i;
    stTree *ret = stTree_cloneNode(tree);
    stPhylogenyInfo *info = stPhylogenyInfo_clone(stTree_getClientData(tree));
    stTree_setClientData(ret, info);
    // Update child partitions (if any)
    for(i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree_setParent(stPhylogeny_scoreFromBootstraps(stTree_getChild(tree, i), bootstraps), ret);
    }

    // Check the current partition against all bootstraps
    for(i = 0; i < stList_length(bootstraps); i++) {
        updateSupportFromTree(ret, tree, stList_get(bootstraps, i),
                              updatePartitionSupportFromPartition);
    }

    info->index->bootstrapSupport = ((double) info->index->numBootstraps) / stList_length(bootstraps);
    return ret;
}

void updateReconciliationSupportFromPartition(stTree *partitionToScore,
                                              stTree *originalPartition,
                                              stTree *bootstrapPartition) {
    stPhylogenyInfo *partitionInfo, *bootstrapInfo;

    // First, check that the branches are identical.
    partitionInfo = stTree_getClientData(partitionToScore);
    bootstrapInfo = stTree_getClientData(bootstrapPartition);
    assert(partitionInfo != NULL);
    assert(bootstrapInfo != NULL);
    stIndexedTreeInfo *partitionIndex = partitionInfo->index;
    stIndexedTreeInfo *bootstrapIndex = bootstrapInfo->index;
    assert(partitionIndex != NULL);
    assert(bootstrapIndex != NULL);
    assert(partitionIndex->totalNumLeaves == bootstrapIndex->totalNumLeaves);
    // Check if the set of leaves is equal in both partitions. If not,
    // the partitions can't be equal.
    if (memcmp(partitionIndex->leavesBelow, bootstrapIndex->leavesBelow,
            partitionIndex->totalNumLeaves * sizeof(char))) {
        return;
    }

    // Now check the reconciliation of the parents.
    stTree *partitionParent = stTree_getParent(originalPartition);
    stTree *bootstrapParent = stTree_getParent(bootstrapPartition);
    if (partitionParent == NULL && bootstrapParent == NULL) {
        // We count this case as having identical reconciliation.
        partitionIndex->numBootstraps++;
        return;
    } else if (partitionParent == NULL || bootstrapParent == NULL) {
        // If only one is the root, we don't consider them to have the
        // same reconciliation.
        return;
    }
    stPhylogenyInfo *partitionParentInfo = stTree_getClientData(partitionParent);
    stPhylogenyInfo *bootstrapParentInfo = stTree_getClientData(bootstrapParent);
    assert(partitionParentInfo != NULL);
    assert(bootstrapParentInfo != NULL);
    stReconciliationInfo *partitionParentRecon = partitionParentInfo->recon;
    stReconciliationInfo *bootstrapParentRecon = bootstrapParentInfo->recon;
    assert(partitionParentRecon != NULL);
    assert(bootstrapParentRecon != NULL);
    if (partitionParentRecon->event != bootstrapParentRecon->event ||
        partitionParentRecon->species != bootstrapParentRecon->species) {
        // Not the same reconciliation / duplication labeling.
        return;
    }
    // The partitions are equal and they have the same reconciliation,
    // increase the support
    partitionIndex->numBootstraps++;
}

stTree *stPhylogeny_scoreReconciliationFromBootstrap(stTree *tree,
                                                     stTree *bootstrap) {
    stList *list = stList_construct();
    stList_append(list, bootstrap);
    stTree *ret = stPhylogeny_scoreReconciliationFromBootstraps(tree, list);
    stList_destruct(list);
    return ret;
}

// Return a new tree which has its partitions scored by how often they
// appear with the same reconciliation in the bootstraps. Obviously,
// this is always less than or equal to the normal bootstrap
// score. This fills in the numBootstraps and bootstrapSupport fields
// of each node. All trees must have valid stPhylogenyInfo.
stTree *stPhylogeny_scoreReconciliationFromBootstraps(stTree *tree,
                                                      stList *bootstraps)
{
    int64_t i;
    stTree *ret = stTree_cloneNode(tree);
    stPhylogenyInfo *info = stPhylogenyInfo_clone(stTree_getClientData(tree));
    stTree_setClientData(ret, info);
    // Update child partitions (if any)
    for(i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree_setParent(stPhylogeny_scoreReconciliationFromBootstraps(stTree_getChild(tree, i), bootstraps), ret);
    }
    // Check the current partition against all bootstraps
    for(i = 0; i < stList_length(bootstraps); i++) {
        updateSupportFromTree(ret, tree, stList_get(bootstraps, i),
                              updateReconciliationSupportFromPartition);
    }

    info->index->bootstrapSupport = ((double) info->index->numBootstraps) / stList_length(bootstraps);
    return ret;
}

// Only one half of the distanceMatrix is used, distances[i][j] for which i > j
// Tree returned is labeled by the indices of the distance matrix. The
// tree is rooted halfway along the longest branch if outgroups is
// NULL; otherwise, it's rooted halfway along the longest outgroup
// branch.
stTree *stPhylogeny_neighborJoin(stMatrix *distances, stList *outgroups) {
    struct DistanceMatrix *distanceMatrix;
    struct Tree *tree;
    int64_t i, j;
    int64_t numSequences = stMatrix_n(distances);
    assert(numSequences > 2);
    assert(distances != NULL);
    assert(stMatrix_n(distances) == stMatrix_m(distances));
    // Set up the basic QuickTree data structures to represent the sequences.
    // The data structures are only filled in as much as absolutely
    // necessary, so they will probably be invalid for anything but
    // running neighbor-joining.
    struct ClusterGroup *clusterGroup = empty_ClusterGroup();
    struct Cluster **clusters = st_malloc(numSequences * sizeof(struct Cluster *));
    for (i = 0; i < numSequences; i++) {
        struct Sequence *seq = empty_Sequence();
        seq->name = stString_print_r("%" PRIi64, i);
        clusters[i] = single_Sequence_Cluster(seq);
    }
    clusterGroup->clusters = clusters;
    clusterGroup->numclusters = numSequences;
    // Fill in the QuickTree distance matrix
    distanceMatrix = empty_DistanceMatrix(numSequences);
    for (i = 0; i < numSequences; i++) {
        for (j = 0; j <= i; j++) {
            distanceMatrix->data[i][j] = *stMatrix_getCell(distances, i, j);
        }
    }
    clusterGroup->matrix = distanceMatrix;
    // Finally, run the neighbor-joining algorithm.
    tree = neighbour_joining_buildtree(clusterGroup, 0);
    free_ClusterGroup(clusterGroup);
    return quickTreeToStTree(tree, outgroups);
}

// Get the distance to a leaf from an internal node
static double stPhylogeny_distToLeaf(stTree *tree, int64_t leafIndex) {
    int64_t i;
    stPhylogenyInfo *info = stTree_getClientData(tree);
    (void)info;
    assert(info->index->leavesBelow[leafIndex]);
    if(stTree_getChildNumber(tree) == 0) {
        return 0.0;
    }
    for(i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        stPhylogenyInfo *childInfo = stTree_getClientData(child);
        if(childInfo->index->leavesBelow[leafIndex]) {
            return stTree_getBranchLength(child) + stPhylogeny_distToLeaf(child, leafIndex);
        }
    }
    // We shouldn't've gotten here--none of the children have the
    // leaf under them, but this node claims to have the leaf under it!
    assert(false);
    return 0.0/0.0;
}

// Get the distance to a node from an internal node above it. Will
// fail if the target node is not below.
static double stPhylogeny_distToChild(stTree *tree, stTree *target) {
    int64_t i, j;
    stPhylogenyInfo *treeInfo, *targetInfo;
    treeInfo = stTree_getClientData(tree);
    targetInfo = stTree_getClientData(target);
    assert(treeInfo != NULL);
    assert(targetInfo != NULL);
    stIndexedTreeInfo *treeIndex = treeInfo->index;
    stIndexedTreeInfo *targetIndex = targetInfo->index;
    assert(treeIndex != NULL);
    assert(targetIndex != NULL);
    assert(treeIndex->totalNumLeaves == targetIndex->totalNumLeaves);
    if(memcmp(treeIndex->leavesBelow, targetIndex->leavesBelow,
              treeIndex->totalNumLeaves) == 0) {
        // This node is the target node
        return 0.0;
    }
    for(i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        stPhylogenyInfo *childInfo = stTree_getClientData(child);
        assert(childInfo != NULL);
        stIndexedTreeInfo *childIndex = childInfo->index;
        assert(childIndex != NULL);
        bool childIsSuperset = true;
        // Go through all the children and find one which is above (or
        // is) the target node. (Any node above the target will have a
        // leaf set that is a superset of the target's leaf set.)
        assert(childIndex->totalNumLeaves == treeIndex->totalNumLeaves);
        for(j = 0; j < childIndex->totalNumLeaves; j++) {
            if(targetIndex->leavesBelow[j] && !childIndex->leavesBelow[j]) {
                childIsSuperset = false;
                break;
            }
        }
        if(childIsSuperset) {
            return stPhylogeny_distToChild(child, target) + stTree_getBranchLength(child);
            break;
        }
    }
    // We shouldn't've gotten here--none of the children have the
    // target under them, but this node claims to have the target under it!
    assert(false);
    return 0.0/0.0;
}


// Return the MRCA of the given leaves.
stTree *stPhylogeny_getMRCA(stTree *tree, int64_t leaf1, int64_t leaf2) {
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        stPhylogenyInfo *childInfo = stTree_getClientData(child);
        assert(childInfo != NULL);
        stIndexedTreeInfo *childIndex = childInfo->index;
        assert(childIndex != NULL);
        if(childIndex->leavesBelow[leaf1] && childIndex->leavesBelow[leaf2]) {
            return stPhylogeny_getMRCA(child, leaf1, leaf2);
        }
    }

    // If we've gotten to this point, then this is the MRCA of the leaves.
    return tree;
}

// Find the distance between leaves (given by their index in the
// distance matrix.)
double stPhylogeny_distanceBetweenLeaves(stTree *tree, int64_t leaf1,
                                         int64_t leaf2) {
    stTree *mrca = stPhylogeny_getMRCA(tree, leaf1, leaf2);
    return stPhylogeny_distToLeaf(mrca, leaf1) + stPhylogeny_distToLeaf(mrca, leaf2);
}

// Find the distance between two arbitrary nodes (which must be in the
// same tree, with stPhylogenyInfo attached properly).
double stPhylogeny_distanceBetweenNodes(stTree *node1, stTree *node2) {
    stPhylogenyInfo *info1, *info2;
    if(node1 == node2) {
        return 0.0;
    }
    info1 = stTree_getClientData(node1);
    info2 = stTree_getClientData(node2);
    stIndexedTreeInfo *index1 = info1->index;
    stIndexedTreeInfo *index2 = info2->index;
    assert(index1->totalNumLeaves == index2->totalNumLeaves);
    // Check if node1 is under node2, vice versa, or if they aren't on
    // the same path to the root
    bool oneAboveTwo = false, twoAboveOne = false, differentSubsets = true;
    for(int64_t i = 0; i < index1->totalNumLeaves; i++) {
        if(index1->leavesBelow[i] && index2->leavesBelow[i]) {
            differentSubsets = false;
        } else if(index2->leavesBelow[i] && !index1->leavesBelow[i]) {
            twoAboveOne = true;
            // Technically we can break here, but it's cheap to
            // double-check that everything is correct.
            assert(differentSubsets || oneAboveTwo == false);
        } else if(index1->leavesBelow[i] && !index2->leavesBelow[i]) {
            oneAboveTwo = true;
            assert(differentSubsets || twoAboveOne == false);
        }
    }
    // If differentSubsets is true, then the values of oneAboveTwo and
    // twoAboveOne don't matter; if differentSubsets is false, exactly
    // one should be true.
    (void) twoAboveOne; // stop gcc complaining about this value,
                        // which is used in asserts
    assert(differentSubsets || (oneAboveTwo ^ twoAboveOne));

    if(differentSubsets) {
        stTree *parent = node1;
        for(;;) {
            parent = stTree_getParent(parent);
            assert(parent != NULL);
            stPhylogenyInfo *parentInfo = stTree_getClientData(parent);
            assert(parentInfo != NULL && parentInfo->index != NULL);
            stIndexedTreeInfo *parentIndex = parentInfo->index;
            assert(parentIndex->totalNumLeaves == index1->totalNumLeaves);
            bool isCommonAncestor = true;
            for(int64_t i = 0; i < parentIndex->totalNumLeaves; i++) {
                if((index1->leavesBelow[i] || index2->leavesBelow[i]) &&
                   !parentIndex->leavesBelow[i]) {
                    isCommonAncestor = false;
                }
            }
            if(isCommonAncestor) {
                // Found the MRCA of both nodes
                break;
            }
        }
        return stPhylogeny_distToChild(parent, node1) +
            stPhylogeny_distToChild(parent, node2);
    } else if(oneAboveTwo) {
        return stPhylogeny_distToChild(node1, node2);
    } else {
        assert(twoAboveOne);
        return stPhylogeny_distToChild(node2, node1);
    }
}

// Gets the (leaf) node corresponding to an index in the distance matrix.
stTree *stPhylogeny_getLeafByIndex(stTree *tree, int64_t leafIndex) {
    int64_t i;
    stPhylogenyInfo *info = stTree_getClientData(tree);
    assert(info != NULL);
    stIndexedTreeInfo *index = info->index;
    assert(index != NULL);
    assert(leafIndex < index->totalNumLeaves);
    if(index->matrixIndex == leafIndex) {
        return tree;
    }
    for(i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        stPhylogenyInfo *childInfo = stTree_getClientData(child);
        assert(childInfo != NULL);
        stIndexedTreeInfo *childIndex = childInfo->index;
        assert(childIndex != NULL);
        assert(index->totalNumLeaves == childIndex->totalNumLeaves);
        if(childIndex->leavesBelow[leafIndex]) {
            return stPhylogeny_getLeafByIndex(child, leafIndex);
        }
    }

    // Shouldn't get here if the stIndexedTreeInfo is set properly
    return NULL;
}

// Helper function for computeJoinCosts. Populates a map from a tree to a unique int.
static void populateSpeciesToIndex(stTree *speciesTree, stHash *speciesToIndex) {
    stList *bfQueue = stList_construct();
    stList_append(bfQueue, speciesTree);
    int64_t curIdx = 0;
    while (stList_length(bfQueue) != 0) {
        stTree *node = stList_pop(bfQueue);
        for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
            stList_append(bfQueue, stTree_getChild(node, i));
        }
        stHash_insert(speciesToIndex, node, stIntTuple_construct1(curIdx++));
    }
    stList_destruct(bfQueue);
}

// Get the number of nodes between a descendant and its ancestor that
// could cause losses, i.e. that have more than one child. (Exclusive
// of both the ancestor and its descendant, so if the descendant is a
// direct child of the ancestor, that counts as 0.)
static int64_t numSkipsToAncestor(stTree *descendant, stTree *ancestor) {
    if (descendant == ancestor) {
        return 0;
    }
    stTree *curNode = descendant;
    int64_t ret = 0;
    while ((curNode = stTree_getParent(curNode)) != NULL) {
        if (curNode == ancestor) {
            return ret;
        }
        if (stTree_getChildNumber(curNode) != 1) {
            ret++;
        }
    }

    // Shouldn't get here.
    assert(false);
    return -INT64_MAX;
}

// Get the number of losses in a (left,right)parent; subtree
static int64_t lossesInSubtree(stTree *parent, stTree *left, stTree *right) {
    int64_t ret = 0;
    ret += numSkipsToAncestor(left, parent) + numSkipsToAncestor(right, parent);
    if ((left == parent || right == parent) && (left != right)) {
        ret++;
    }
    return ret;
}

// Compute join costs for a species tree for use in guided
// neighbor-joining. These costs are calculated by penalizing
// according to the number of dups and losses implied by the
// reconciliation when joining two genes reconciled to certain nodes
// of the species tree.
// speciesToIndex (a blank hash) will be populated with stIntTuples
// corresponding to each species' index into the join cost matrix.
// NB: the species tree must be binary.
// TODO: take branch lengths into account
stMatrix *stPhylogeny_computeJoinCosts(stTree *speciesTree, stHash *speciesToIndex, double costPerDup, double costPerLoss) {
    assert(stHash_size(speciesToIndex) == 0);

    int64_t numSpecies = stTree_getNumNodes(speciesTree);

    populateSpeciesToIndex(speciesTree, speciesToIndex);
    assert(stHash_size(speciesToIndex) == numSpecies);

    // Fill in the join cost matrix.
    stMatrix *ret = stMatrix_construct(numSpecies, numSpecies);
    stHash *indexToSpecies = stHash_invert(speciesToIndex, (uint64_t (*)(const void *)) stIntTuple_hashKey,
                                           (int (*)(const void *, const void *)) stIntTuple_equalsFn, NULL, NULL);
    for (int64_t i = 0; i < numSpecies; i++) {
        // get the species node for this index (the intTuple is
        // necessary to query the hash properly)
        stIntTuple *query_i = stIntTuple_construct1(i);
        stTree *species_i = stHash_search(indexToSpecies, query_i);
        assert(species_i != NULL);
        for (int64_t j = i; j < numSpecies; j++) {
            // get the species node for this index (the intTuple is
            // necessary to query the hash properly)
            stIntTuple *query_j = stIntTuple_construct1(j);
            stTree *species_j = stHash_search(indexToSpecies, query_j);
            assert(species_j != NULL);

            // Can't use stPhylogeny_getMRCA as that is only defined for leaves.
            stTree *mrca = stTree_getMRCA(species_i, species_j);

            // Calculate the number of dups implied when joining species i and j.
            if (species_i == mrca || species_j == mrca) {
                // One species is the ancestor of the other, or they
                // are equal. This implies one dup.
                *stMatrix_getCell(ret, i, j) += costPerDup;
                if (j != i) {
                    *stMatrix_getCell(ret, j, i) += costPerDup;
                }
            }

            // Calculate the minimum number of losses implied when
            // joining species i and j.
            int64_t numLosses = lossesInSubtree(mrca, species_i, species_j);
            *stMatrix_getCell(ret, i, j) += costPerLoss * numLosses;
            if (j != i) {
                *stMatrix_getCell(ret, j, i) += costPerLoss * numLosses;
            }

            stIntTuple_destruct(query_j);
        }
        stIntTuple_destruct(query_i);
    }

    stHash_destruct(indexToSpecies);
    return ret;
}

int64_t **stPhylogeny_getMRCAMatrix(stTree *speciesTree, stHash *speciesToIndex) {
    int64_t numSpecies = stTree_getNumNodes(speciesTree);
    int64_t **ret = st_calloc(numSpecies, sizeof(int64_t *));
    for (int64_t i = 0; i < numSpecies; i++) {
        ret[i] = st_calloc(numSpecies, sizeof(int64_t));
    }
    stHash *indexToSpecies = stHash_invert(speciesToIndex, (uint64_t (*)(const void *)) stIntTuple_hashKey,
                                           (int (*)(const void *, const void *)) stIntTuple_equalsFn, NULL, NULL);
    for (int64_t i = 0; i < numSpecies; i++) {
        stIntTuple *query_i = stIntTuple_construct1(i);
        for (int64_t j = i; j < numSpecies; j++) {
            stTree *node_i = stHash_search(indexToSpecies, query_i);
            assert(node_i != NULL);
            stIntTuple *query_j = stIntTuple_construct1(j);
            stTree *node_j = stHash_search(indexToSpecies, query_j);
            assert(node_j != NULL);
            stTree *mrca = stTree_getMRCA(node_i, node_j);
            stIntTuple *mrcaIndex = stHash_search(speciesToIndex, mrca);
            assert(mrcaIndex != NULL);
            ret[i][j] = stIntTuple_get(mrcaIndex, 0);
            ret[j][i] = ret[i][j];
            stIntTuple_destruct(query_j);
        }
        stIntTuple_destruct(query_i);
    }
    stHash_destruct(indexToSpecies);
    return ret;
}

stTree *stPhylogeny_guidedNeighborJoining(stMatrix *distanceMatrix,
                                          stMatrix *similarityMatrix,
                                          stMatrix *joinCosts,
                                          stHash *matrixIndexToJoinCostIndex,
                                          stHash *speciesToJoinCostIndex,
                                          int64_t **speciesMRCAMatrix,
                                          stTree *speciesTree) {
    int64_t numLeaves = stMatrix_n(similarityMatrix);
    assert(numLeaves == stMatrix_m(similarityMatrix));
    assert(numLeaves >= 3);

    // Stores the reconciliation (in join cost matrix indices) for
    // each leaf/node.
    int64_t *recon = st_calloc(numLeaves, sizeof(int64_t));

    // Fill in the initial reconciliation array.
    stHashIterator *hashIt = stHash_getIterator(matrixIndexToJoinCostIndex);
    stIntTuple *matrixIndex;
    while ((matrixIndex = stHash_getNext(hashIt)) != NULL) {
        int64_t i = stIntTuple_get(matrixIndex, 0);
        stIntTuple *joinCostIndex = stHash_search(matrixIndexToJoinCostIndex, matrixIndex);
        assert(joinCostIndex != NULL);
        recon[i] = stIntTuple_get(joinCostIndex, 0);
    }
    stHash_destructIterator(hashIt);

    // Distance matrix. We clone the one fed to us since we need to
    // modify it during the process. Note: only valid for i < j.
    double **distances = st_calloc(numLeaves, sizeof(double *));
    for (int64_t i = 0; i < numLeaves; i++) {
        distances[i] = st_calloc(numLeaves, sizeof(double));
    }

    // Confidence matrix (i.e. count of similarities + differences for
    // i,j). Note: only valid for i < j.
    double **confidences = st_calloc(numLeaves, sizeof(double *));
    for (int64_t i = 0; i < numLeaves; i++) {
        confidences[i] = st_calloc(numLeaves, sizeof(double));
    }

    // Fill in our copy of the distance matrix.
    for (int64_t i = 0; i < numLeaves; i++) {
        for (int64_t j = i + 1; j < numLeaves; j++) {
            distances[i][j] = *stMatrix_getCell(distanceMatrix, i, j);
        }
    }

    // Initial "r" cost (the average of the distances from each node
    // to all others) to weight with.
    double *r = st_calloc(numLeaves, sizeof(double));
    for (int64_t i = 0; i < numLeaves; i++) {
        assert(r[i] == 0.0);
        for (int64_t j = 0; j < numLeaves; j++) {
            if (i == j) {
                continue;
            }
            if (i < j) {
                r[i] += distances[i][j];
            } else {
                r[i] += distances[j][i];
            }
        }
        r[i] /= numLeaves - 2;
    }

    // Initialize the join distances (i.e. distance-relative rather
    // than the similarity/difference-relative join costs.) Only
    // defined for row < col.
    double **joinDistances = st_calloc(numLeaves, sizeof(double *));
    for (int64_t i = 0; i < numLeaves; i++) {
        joinDistances[i] = st_calloc(numLeaves, sizeof(double));
        for (int64_t j = i + 1; j < numLeaves; j++) {
            if (confidences[i][j] != 0) {
                joinDistances[i][j] = *stMatrix_getCell(joinCosts, recon[i], recon[j]) / confidences[i][j];
            } else {
                // the distance will already be INT64_MAX, just try to
                // bias it toward the species tree (although we may
                // run into floating-point resolution errors here)
                joinDistances[i][j] = *stMatrix_getCell(joinCosts, recon[i], recon[j]);
            }
        }
    }

    // The actual neighbor-joining process.
    stTree **nodes = st_calloc(numLeaves, sizeof(stTree *));
    for (int64_t i = 0; i < numLeaves; i++) {
        // Initialize nodes.
        nodes[i] = stTree_construct();
        char *name = stString_print_r("%" PRIi64, i);
        stTree_setLabel(nodes[i], name);
        free(name);
    }
    int64_t numJoinsLeft = numLeaves - 1;
    while (numJoinsLeft > 0) {
        // Find the lowest distance between any two roots in the
        // forest.
        double minDist = DBL_MAX;
        int64_t mini = -1, minj = -1;
        for (int64_t i = 0; i < numLeaves; i++) {
            if (recon[i] == -1) {
                // Signals that this node has been joined and its
                // index is abandoned.
                continue;
            }
            for (int64_t j = i + 1; j < numLeaves; j++) {
                if (recon[j] == -1) {
                    // Signals that this node has been joined and its
                    // index is abandoned.
                    continue;
                }
                double dist = distances[i][j] + joinDistances[i][j] - r[i] - r[j];
                if (dist < minDist) {
                    mini = i;
                    minj = j;
                    minDist = dist;
                }
            }
        }
        assert(mini != -1);
        assert(minj != -1);

        double dist_mini_minj = distances[mini][minj];

        // Get the branch lengths for the children of the new node.
        double branchLength_mini = (dist_mini_minj + r[mini] - r[minj]) / 2;
        double branchLength_minj = dist_mini_minj - branchLength_mini;
        // Fix the distances in case of negative branch length.
        if ((branchLength_mini <= 0 || branchLength_minj <= 0) && dist_mini_minj < 0) {
            branchLength_mini = 0;
            branchLength_minj = 0;
        } else if (branchLength_mini < 0) {
            branchLength_mini = 0;
            branchLength_minj = dist_mini_minj;
        } else if (branchLength_minj < 0) {
            branchLength_mini = dist_mini_minj;
            branchLength_minj = 0;
        }

        // Join the nodes.
        stTree *nodei = nodes[mini];
        assert(nodei != NULL);
        stTree *nodej = nodes[minj];
        assert(nodej != NULL);
        stTree *joined = stTree_construct();
        stTree_setParent(nodei, joined);
        stTree_setParent(nodej, joined);
        stTree_setBranchLength(nodei, branchLength_mini);
        stTree_setBranchLength(nodej, branchLength_minj);
        nodes[mini] = joined;
        // Not strictly necessary.
        nodes[minj] = NULL;

        // Update the reconciliation of our new joined node to be the
        // MRCA of its children's reconciliations.
        int64_t recon_i = recon[mini];
        int64_t recon_j = recon[minj];
        recon[mini] = speciesMRCAMatrix[recon_i][recon_j];
        recon[minj] = -1;

        // Update the new row of the distance matrix.
        for (int64_t k = 0; k < numLeaves; k++) {
            if (recon[k] == -1) {
                // Node is gone, don't need to update.
                continue;
            }
            if (k == mini) {
                // Skip distance calculation with self
                continue;
            }

            // calculate the proper indices for accessing i<->k and
            // j<->k distances and confidences (the matrix is
            // symmetric, so only one half is stored)
            int64_t mini_kRow = mini;
            int64_t mini_kCol = k;
            int64_t minj_kRow = minj;
            int64_t minj_kCol = k;
            if (mini > k) {
                mini_kRow = k;
                mini_kCol = mini;
            }
            if (minj > k) {
                minj_kRow = k;
                minj_kCol = minj;
            }

            confidences[mini_kRow][mini_kCol] = (confidences[mini_kRow][mini_kCol] + confidences[minj_kRow][minj_kCol]) / 2;

            // Update the distances.
            double dist_mini_k = distances[mini_kRow][mini_kCol];
            double dist_minj_k = distances[minj_kRow][minj_kCol];
            distances[mini_kRow][mini_kCol] = (dist_mini_k + dist_minj_k - dist_mini_minj) / 2;

            // Update the join distance.
            if (confidences[mini_kRow][mini_kCol] != 0) {
                joinDistances[mini_kRow][mini_kCol] = *stMatrix_getCell(joinCosts, recon[mini], recon[k]) / confidences[mini_kRow][mini_kCol];
            } else {
                joinDistances[mini_kRow][mini_kCol] = *stMatrix_getCell(joinCosts, recon[mini], recon[k]);
            }

            // Update r[k].
            if (numJoinsLeft > 2) {
                r[k] = ((r[k] * (numJoinsLeft - 1)) - dist_mini_k - dist_minj_k + distances[mini_kRow][mini_kCol]) / (numJoinsLeft - 2);
            } else {
                r[k] = 0.0;
            }
        }

        // Set r for the new column.
        r[mini] = 0.0;
        if (numJoinsLeft > 2) {
            for (int64_t k = 0; k < numLeaves; k++) {
                if (recon[k] == -1) {
                    continue;
                }
                if (k < mini) {
                    r[mini] += distances[k][mini];
                } else {
                    r[mini] += distances[mini][k];
                }
            }
            r[mini] /= numJoinsLeft - 2;
        } else {
            r[mini] = 0.0;
        }
        numJoinsLeft--;
    }

    stTree *ret = nodes[0];
    assert(ret != NULL);

    // Clean up.
    free(recon);
    free(r);
    free(nodes);

    for (int64_t i = 0; i < numLeaves; i++) {
        free(distances[i]);
    }
    free(distances);

    for (int64_t i = 0; i < numLeaves; i++) {
        free(confidences[i]);
    }
    free(confidences);

    for (int64_t i = 0; i < numLeaves; i++) {
        free(joinDistances[i]);
    }
    free(joinDistances);

    assert(stTree_getNumNodes(ret) == numLeaves * 2 - 1);

    stPhylogeny_addStIndexedTreeInfo(ret);
    return ret;
}

// Fills in stReconciliationInfo, creating the containing
// stPhylogenyInfo if necessary.
static void fillInReconciliationInfo(stTree *gene, stTree *recon,
                                     stReconciliationEvent event,
                                     bool relabelAncestors) {
    stPhylogenyInfo *info = stTree_getClientData(gene);
    if (info == NULL) {
        info = st_calloc(1, sizeof(stPhylogenyInfo));
        stTree_setClientData(gene, info);
    }
    stReconciliationInfo *reconInfo = info->recon;
    if (reconInfo == NULL) {
        reconInfo = st_calloc(1, sizeof(stReconciliationInfo));
        info->recon = reconInfo;
    }
    reconInfo->species = recon;
    reconInfo->event = event;
    if (stTree_getChildNumber(gene) != 0 && relabelAncestors) {
        stTree_setLabel(gene, stTree_getLabel(recon));
    }
}

static stTree *stPhylogeny_reconcileAtMostBinary_R(stTree *gene,
                                                   stHash *leafToSpecies,
                                                   bool relabelAncestors) {
    stTree *recon;
    stReconciliationEvent event;
    if (stTree_getChildNumber(gene) == 0) {
        // Leaves are already reconciled.
        recon = stHash_search(leafToSpecies, gene);
        assert(recon != NULL);
        event = LEAF;
    } else {
        event = SPECIATION;
        recon = stPhylogeny_reconcileAtMostBinary_R(
            stTree_getChild(gene, 0), leafToSpecies, relabelAncestors);
        for (int64_t i = 1; i < stTree_getChildNumber(gene); i++) {
            stTree *childRecon = stPhylogeny_reconcileAtMostBinary_R(
                stTree_getChild(gene, i), leafToSpecies, relabelAncestors);
            recon = stTree_getMRCA(childRecon, recon);
        }
        for (int64_t i = 0; i < stTree_getChildNumber(gene); i++) {
            stPhylogenyInfo *childInfo = stTree_getClientData(stTree_getChild(gene, i));
            stTree *childRecon = childInfo->recon->species;
            if (childRecon == recon) {
                event = DUPLICATION;
            }
        }
    }
    fillInReconciliationInfo(gene, recon, event, relabelAncestors);
    return recon;
}

// Reconcile a gene tree (without rerooting), set the proper
// stReconcilationInfo (as an entry of stPhylogenyInfo) as client data
// on all nodes, and optionally set the labels of the ancestors to the
// labels of their reconciliation in the species tree.
//
// The gene tree must be binary, and the species tree must be
// "at-most-binary", i.e. it must have no nodes with more than 3
// children, but may have nodes with only one child.
void stPhylogeny_reconcileAtMostBinary(stTree *geneTree, stHash *leafToSpecies,
                                       bool relabelAncestors) {
    stPhylogeny_reconcileAtMostBinary_R(geneTree, leafToSpecies,
                                        relabelAncestors);
}

static bool getLinkedSpeciesTree_R(stTree *speciesNode, stTree *polytomy, stHash *speciesToNumGenes, stTree *linkedNode) {
    bool hasGene = false;
    for (int64_t i = 0; i < stTree_getChildNumber(polytomy); i++) {
        stTree *child = stTree_getChild(polytomy, i);
        stPhylogenyInfo *childInfo = stTree_getClientData(child);
        stReconciliationInfo *childRecon = childInfo->recon;
        stTree *childSpecies = childRecon->species;
        if (childSpecies == speciesNode) {
            hasGene = true;
            int64_t *numGenes = stHash_search(speciesToNumGenes, linkedNode);
            if (numGenes == NULL) {
                numGenes = calloc(1, sizeof(int64_t));
                stHash_insert(speciesToNumGenes, linkedNode, numGenes);
            }
            (*numGenes)++;
        }
    }
    if (stTree_getChildNumber(speciesNode) == 0 && !hasGene) {
        return false;
    }
    bool descendantHasGene = false;
    for (int64_t i = 0; i < stTree_getChildNumber(speciesNode); i++) {
        if (getLinkedSpeciesTree_R(stTree_getChild(speciesNode, i), polytomy,
                                   speciesToNumGenes, stTree_getChild(linkedNode, i))) {
            descendantHasGene = true;
        }
    }

    if (!descendantHasGene) {
        // have to be careful to always destruct the 0th child here:
        // we can't iterate the list of children while destroying the
        // children.
        int64_t numChildren = stTree_getChildNumber(linkedNode);
        for (int64_t i = 0; i < numChildren; i++) {
            stTree *linkedChild = stTree_getChild(linkedNode, 0);
            stTree_setParent(linkedChild, NULL);
            stTree_destruct(linkedChild);
        }
    }
    return descendantHasGene || hasGene;
}

// Get the "linked species tree" for a polytomy: i.e. the part of the
// species tree corresponding to the polytomy.
stTree *getLinkedSpeciesTree(stTree *speciesTree, stTree *polytomy, stHash **speciesToNumGenes) {
    *speciesToNumGenes = stHash_construct2(NULL, free);
    stTree *linkedTree = stTree_clone(speciesTree);
    getLinkedSpeciesTree_R(speciesTree, polytomy, *speciesToNumGenes, linkedTree);
    return linkedTree;
}

// Calculate the minimum cost for an x-partial resolution of a
// polytomy.
static int64_t minCost(stTree *species, int64_t x, stHash *cupValues) {
    stIntTuple *cupParameters = stHash_search(cupValues, species);
    assert(cupParameters != NULL);
    int64_t m1 = stIntTuple_get(cupParameters, 0);
    int64_t m2 = stIntTuple_get(cupParameters, 1);
    int64_t gamma = stIntTuple_get(cupParameters, 2);

    if (x < m1) {
        return gamma + m1 - x;
    } else if (x >= m1 && x <= m2) {
        return gamma;
    } else {
        assert (x > m2);
        return gamma + x - m2;
    }
}

static void getCupValues_R(stTree *species, stHash *speciesToNumGenes, stHash *cupValues) {
    for (int64_t i = 0; i < stTree_getChildNumber(species); i++) {
        getCupValues_R(stTree_getChild(species, i), speciesToNumGenes, cupValues);
    }

    int64_t *numGenesPtr = stHash_search(speciesToNumGenes, species);
    int64_t numGenes = numGenesPtr ? *numGenesPtr : 0;
    int64_t m1, m2, gamma;
    if (stTree_getChildNumber(species) == 0) {
        if (numGenes == 0) {
            m1 = 1;
            m2 = 1;
            gamma = 1;
        } else {
            m1 = numGenes;
            m2 = numGenes;
            gamma = 0;
        }
    } else {
        assert(stTree_getChildNumber(species) == 2);
        stTree *child1 = stTree_getChild(species, 0);
        stTree *child2 = stTree_getChild(species, 1);
        stIntTuple *child1Tuple = stHash_search(cupValues, child1);
        stIntTuple *child2Tuple = stHash_search(cupValues, child2);
        int64_t child1m1 = stIntTuple_get(child1Tuple, 0);
        int64_t child1m2 = stIntTuple_get(child1Tuple, 1);
        int64_t child1gamma = stIntTuple_get(child1Tuple, 2);
        int64_t child2m1 = stIntTuple_get(child2Tuple, 0);
        int64_t child2m2 = stIntTuple_get(child2Tuple, 1);
        int64_t child2gamma = stIntTuple_get(child2Tuple, 2);

        if (child1m1 < child2m1 && child1m2 < child2m1) {
            m1 = child1m2;
            m2 = child2m1;
            gamma = child1gamma + child2gamma + child2m1 - child1m2;
        } else if (child1m1 < child2m1 && child1m2 >= child2m1 && child1m2 <= child2m2) {
            m1 = child2m1;
            m2 = child1m2;
            gamma = child1gamma + child2gamma;
        } else if (child1m1 < child2m1 && child1m2 > child2m2) {
            m1 = child2m1;
            m2 = child2m2;
            gamma = child1gamma + child2gamma;
        } else if (child2m1 <= child1m1 && child1m1 <= child2m2 && child2m1 <= child1m2 && child1m2 <= child2m2) {
            // The first part of this if condition is duplicated in
            // the paper, but it is clear that they mean the case
            // where child1's breakpoints are totally enclosed within
            // child2's breakpoints.
            m1 = child1m1;
            m2 = child1m2;
            gamma = child1gamma + child2gamma;
        } else if (child2m1 <= child1m1 && child1m1 <= child2m2 && child1m2 > child2m2) {
            m1 = child1m2;
            m2 = child2m2;
            gamma = child1gamma + child2gamma;
        } else {
            m1 = child2m2;
            m2 = child1m1;
            gamma = child1gamma + child2gamma + child1m1 - child2m2;
        }

        m1 = m1 + numGenes;
        m2 = m2 + numGenes;
        if (m1 <= numGenes) {
            m1 = numGenes + 1;
        }
        if (m2 <= numGenes) {
            m2 = numGenes + 1;
            gamma = minCost(child1, 1, cupValues) + minCost(child2, 1, cupValues);
        }
    }
    stHash_insert(cupValues, species, stIntTuple_construct3(m1, m2, gamma));
}

// Get the "cup values" (i.e. parameters for the minCost function) for
// a polytomy.
static stHash *getCupValues(stTree *speciesTree, stHash *speciesToNumGenes) {
    stHash *cupValues = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
    getCupValues_R(speciesTree, speciesToNumGenes, cupValues);
    return cupValues;
}

// Get the number of duplications and losses in a minimum-cost
// resolution of a polytomy at each node in the species tree.
static void dupLoss(stTree *species, int64_t k, stHash *speciesToNumGenes, stHash *cupValues, stHash *dups, stHash *losses) {
    int64_t *numGenesPtr = stHash_search(speciesToNumGenes, species);
    int64_t numGenes = numGenesPtr ? *numGenesPtr : 0;

    stIntTuple *cupParameters = stHash_search(cupValues, species);
    assert(cupParameters != NULL);
    int64_t m1 = stIntTuple_get(cupParameters, 0);
    int64_t m2 = stIntTuple_get(cupParameters, 1);

    stTree *child1 = stTree_getChildNumber(species) > 0 ? stTree_getChild(species, 0) : NULL;
    stTree *child2 = stTree_getChildNumber(species) > 1 ? stTree_getChild(species, 1) : NULL;

    int64_t *numDups = stHash_search(dups, species);
    if (numDups == NULL) {
        numDups = st_calloc(1, sizeof(int64_t));
        stHash_insert(dups, species, numDups);
    }
    int64_t *numLosses = stHash_search(losses, species);
    if (numLosses == NULL) {
        numLosses = st_calloc(1, sizeof(int64_t));
        stHash_insert(losses, species, numLosses);
    }

    if (stTree_getChildNumber(species) == 0) {
        if (k >= numGenes) {
            *numDups = 0;
            *numLosses = k - numGenes;
        } else {
            *numDups = numGenes - k;
            *numLosses = 0;
        }
    } else if (k - numGenes > 0 && minCost(species, k, cupValues) == minCost(child1, k - numGenes, cupValues) + minCost(child2, k - numGenes, cupValues)) {
        *numDups = 0;
        *numLosses = 0;
        dupLoss(child1, k - numGenes, speciesToNumGenes, cupValues, dups, losses);
        dupLoss(child2, k - numGenes, speciesToNumGenes, cupValues, dups, losses);
    } else if (k < m1) {
        *numDups = m1 - k;
        *numLosses = 0;
        dupLoss(child1, m1 - numGenes, speciesToNumGenes, cupValues, dups, losses);
        dupLoss(child2, m1 - numGenes, speciesToNumGenes, cupValues, dups, losses);
    } else if (k > m2) {
        *numDups = 0;
        *numLosses = k - m2;
        dupLoss(child1, m2 - numGenes, speciesToNumGenes, cupValues, dups, losses);
        dupLoss(child2, m2 - numGenes, speciesToNumGenes, cupValues, dups, losses);
    }
}

// For a tree that has already been reconciled by
// reconcileAtMostBinary, calculates the number of dups and losses
// implied by the reconciliation. dups and losses must be set to 0
// before calling.
void stPhylogeny_reconciliationCostAtMostBinary(stTree *reconciledTree,
                                                int64_t *dups,
                                                int64_t *losses) {
    stPhylogenyInfo *info = stTree_getClientData(reconciledTree);
    assert(info != NULL);
    stReconciliationInfo *recon = info->recon;
    assert(recon != NULL);
    stTree *species = recon->species;
    if (stTree_getChildNumber(reconciledTree) == 2) {
        if (recon->event == DUPLICATION) {
            (*dups)++;
        }
        // Count losses.
        // Look at all this node's children. If the children's recons
        // aren't direct children of this node's recon, then count N
        // losses, where N is the number of nodes "skipped" on the way to
        // this node's recon (plus one if this is a dup and the children's
        // recons are not equal). Nodes that have only one child should
        // not count as "skipped".
        stTree *leftChild = stTree_getChild(reconciledTree, 0);
        stPhylogenyInfo *leftInfo = stTree_getClientData(leftChild);
        assert(leftInfo != NULL);
        stReconciliationInfo *leftRecon = leftInfo->recon;
        assert(leftRecon != NULL);
        stTree *leftSpecies = leftRecon->species;

        stTree *rightChild = stTree_getChild(reconciledTree, 1);
        stPhylogenyInfo *rightInfo = stTree_getClientData(rightChild);
        assert(rightInfo != NULL);
        stReconciliationInfo *rightRecon = rightInfo->recon;
        assert(rightRecon != NULL);
        stTree *rightSpecies = rightRecon->species;

        *losses += lossesInSubtree(species, leftSpecies, rightSpecies);
    } else if (stTree_getChildNumber(reconciledTree) > 2) {
        // We follow the algorithm for finding the minimum mutation
        // cost in an apparent polytomy given by Lafond, Swenson,
        // El-Mabrouk, 2012.
        stHash *speciesToNumGenes;
        stTree *linkedSpeciesTree = getLinkedSpeciesTree(species, reconciledTree, &speciesToNumGenes);
        stHash *cupValues = getCupValues(linkedSpeciesTree, speciesToNumGenes);
        stHash *dupsPerSpecies = stHash_construct2(NULL, free);
        stHash *lossesPerSpecies = stHash_construct2(NULL, free);
        dupLoss(linkedSpeciesTree, 1, speciesToNumGenes, cupValues, dupsPerSpecies, lossesPerSpecies);
        stList *dupValues = stHash_getValues(dupsPerSpecies);
        for (int64_t i = 0; i < stList_length(dupValues); i++) {
            int64_t *dupsInSpecies = stList_get(dupValues, i);
            *dups += *dupsInSpecies;
        }

        stList *lossValues = stHash_getValues(lossesPerSpecies);
        for (int64_t i = 0; i < stList_length(lossValues); i++) {
            int64_t *lossesInSpecies = stList_get(lossValues, i);
            *losses += *lossesInSpecies;
        }

        stList_destruct(dupValues);
        stList_destruct(lossValues);
        stHash_destruct(lossesPerSpecies);
        stHash_destruct(dupsPerSpecies);
        stHash_destruct(cupValues);
        stTree_destruct(linkedSpeciesTree);
        stHash_destruct(speciesToNumGenes);
    }

    for (int64_t i = 0; i < stTree_getChildNumber(reconciledTree); i++) {
        stPhylogeny_reconciliationCostAtMostBinary(stTree_getChild(reconciledTree, i),
                                                   dups, losses);
    }
}

// Recurse down a tree testing roots to see which would give the
// lowest recon cost if the tree was rooted at that position.
// curRoot is the child of the branch to root on.
void stPhylogeny_rootByReconciliationAtMostBinary_R(stTree *curRoot,
                                                    stTree *prevRootParentSpecies,
                                                    int64_t prevRootDups,
                                                    int64_t prevRootLosses,
                                                    int64_t *bestDups,
                                                    int64_t *bestLosses,
                                                    stTree **bestRoot) {
    // The difference between the tree rooted at this branch and the
    // one that was rooted previous to this is just the reconciliation
    // of the parent of this branch and the new root. So recalculate
    // the reconciliation of these nodes, pretending that the tree has
    // already been rerooted.
    stTree *parent = stTree_getParent(curRoot);
    stPhylogenyInfo *parentInfo = stTree_getClientData(parent);
    assert(parentInfo != NULL && parentInfo->recon != NULL);
    stTree *parentOldSpecies = parentInfo->recon->species;

    // First, the parent of our branch should just be reconciled to
    // the MRCA of prev root's parent reconciliation and our sibling's
    // reconciliation (which has not changed).
    stTree *sibling = NULL;
    for (int64_t i = 0; i < stTree_getChildNumber(parent); i++) {
        if (stTree_getChild(parent, i) != curRoot) {
            sibling = stTree_getChild(parent, i);
        }
    }
    stPhylogenyInfo *siblingInfo = stTree_getClientData(sibling);
    assert(siblingInfo != NULL && siblingInfo->recon != NULL);
    stTree *siblingSpecies = siblingInfo->recon->species;
    // Call this the parent's new species for consistency, although
    // now it's the sibling in the new tree. The name is confusing
    // either way.
    stTree *parentNewSpecies = stTree_getMRCA(prevRootParentSpecies, siblingSpecies);

    // Next, the new root's recon is just the MRCA of our parent's
    // recon in the rerooted tree, and the recon of this node (which
    // would stay the same).
    stPhylogenyInfo *curInfo = stTree_getClientData(curRoot);
    assert(curInfo != NULL && curInfo->recon != NULL);
    stTree *curSpecies = curInfo->recon->species;
    stTree *newRootSpecies = stTree_getMRCA(curSpecies, parentNewSpecies);

    // Find the new cost in dups. This is just (# of old dups) - (old root
    // was dup? 1 : 0) - (parent used to be dup? 1 : 0) + (new root is
    // dup? 1 : 0) + (parent (now sibling) is now dup? 1 : 0)
    int64_t curRootDups = prevRootDups;
    if (parentOldSpecies == curSpecies || parentOldSpecies == siblingSpecies) {
        curRootDups--;
    }
    stTree *oldRootSpecies = stTree_getMRCA(prevRootParentSpecies,
                                            parentOldSpecies);
    if (oldRootSpecies == prevRootParentSpecies || oldRootSpecies == parentOldSpecies) {
        curRootDups--;
    }
    if (newRootSpecies == curSpecies || newRootSpecies == parentNewSpecies) {
        curRootDups++;
    }
    if (parentNewSpecies == siblingSpecies || parentNewSpecies == prevRootParentSpecies) {
        curRootDups++;
    }

    // Now find the new cost in losses.
    int64_t curRootLosses = prevRootLosses;
    curRootLosses -= lossesInSubtree(parentOldSpecies, curSpecies, siblingSpecies);
    curRootLosses -= lossesInSubtree(oldRootSpecies, prevRootParentSpecies, parentOldSpecies);
    assert(curRootLosses >= 0);
    curRootLosses += lossesInSubtree(newRootSpecies, curSpecies, parentNewSpecies);
    curRootLosses += lossesInSubtree(parentNewSpecies, siblingSpecies, prevRootParentSpecies);

    if (curRootDups < *bestDups || (curRootDups == *bestDups && curRootLosses < *bestLosses)) {
        *bestDups = curRootDups;
        *bestLosses = curRootLosses;
        *bestRoot = curRoot;
    }
    for (int64_t i = 0; i < stTree_getChildNumber(curRoot); i++) {
        stPhylogeny_rootByReconciliationAtMostBinary_R(stTree_getChild(curRoot, i),
                                                       parentNewSpecies,
                                                       curRootDups,
                                                       curRootLosses, bestDups,
                                                       bestLosses, bestRoot);
    }
}

// Return a copy of geneTree that is rooted to minimize duplications,
// with the amount of losses as tiebreaker between roots.
// NOTE: the returned tree does *not* have reconciliation info set,
// and this function will reconcile geneTree, resetting any
// reconciliation information that potentially already exists.
stTree *stPhylogeny_rootByReconciliationAtMostBinary(stTree *geneTree,
                                                     stHash *leafToSpecies) {
    stPhylogeny_reconcileAtMostBinary(geneTree, leafToSpecies, false);

    // Find the root which has the lowest reconciliation cost.
    int64_t dups = 0, losses = 0;
    stPhylogeny_reconciliationCostAtMostBinary(geneTree, &dups, &losses);
    stTree *bestRoot = geneTree;
    int64_t bestDups = dups;
    int64_t bestLosses = losses;
    if (stTree_getChildNumber(geneTree) == 0) {
        return stTree_clone(geneTree);
    } else {
        assert(stTree_getChildNumber(geneTree) == 2);
        stTree *leftChild = stTree_getChild(geneTree, 0);
        stTree *rightChild = stTree_getChild(geneTree, 1);
        stPhylogenyInfo *leftChildInfo = stTree_getClientData(leftChild);
        stTree *leftChildSpecies = leftChildInfo->recon->species;
        stPhylogenyInfo *rightChildInfo = stTree_getClientData(rightChild);
        stTree *rightChildSpecies = rightChildInfo->recon->species;
        for (int64_t i = 0; i < stTree_getChildNumber(leftChild); i++) {
            stPhylogeny_rootByReconciliationAtMostBinary_R(stTree_getChild(leftChild, i),
                                                           rightChildSpecies,
                                                           dups, losses,
                                                           &bestDups,
                                                           &bestLosses,
                                                           &bestRoot);
        }
        for (int64_t i = 0; i < stTree_getChildNumber(rightChild); i++) {
            stPhylogeny_rootByReconciliationAtMostBinary_R(stTree_getChild(rightChild, i),
                                                           leftChildSpecies,
                                                           dups, losses,
                                                           &bestDups,
                                                           &bestLosses,
                                                           &bestRoot);
        }
        return stTree_reRoot(bestRoot, stTree_getBranchLength(bestRoot)/2);
    }
}

static void getNewLeafToSpecies_R(stTree *node, stHash *leafToSpecies) {
    if (stTree_getChildNumber(node) == 0) {
        stPhylogenyInfo *info = stTree_getClientData(node);
        assert(info != NULL);
        stReconciliationInfo *recon = info->recon;
        assert(recon != NULL);
        stHash_insert(leafToSpecies, node, recon->species);
    }
    for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
        getNewLeafToSpecies_R(stTree_getChild(node, i), leafToSpecies);
    }
    stTree_setClientData(node, NULL);
}

static stHash *getNewLeafToSpecies(stTree *rerooted) {
    stHash *leafToSpecies = stHash_construct();
    getNewLeafToSpecies_R(rerooted, leafToSpecies);
    return leafToSpecies;
}

stTree *stPhylogeny_rootByReconciliationNaive(stTree *tree, stHash *leafToSpecies) {
    stPhylogeny_reconcileAtMostBinary(tree, leafToSpecies, false);
    stList *stack = stList_construct();
    stList_append(stack, tree);
    stTree *bestTree = NULL;
    int64_t bestDups = INT64_MAX;
    int64_t bestLosses = INT64_MAX;
    while (stList_length(stack) > 0) {
        stTree *node = stList_pop(stack);
        for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
            stList_append(stack, stTree_getChild(node, i));
        }

        stTree *curTree;
        if (node != tree) {
            curTree = stTree_reRootAndKeepClientData(node, stTree_getBranchLength(node)/2);
        } else {
            curTree = stTree_clone(tree);
        }
        stHash *newLeafToSpecies = getNewLeafToSpecies(curTree);
        stPhylogeny_reconcileAtMostBinary(curTree, newLeafToSpecies, false);
        stHash_destruct(newLeafToSpecies);
        int64_t dups = 0, losses = 0;
        stPhylogeny_reconciliationCostAtMostBinary(curTree, &dups, &losses);
        if (dups < bestDups || (dups == bestDups && losses < bestLosses)) {
            if (bestTree != NULL) {
                stPhylogenyInfo_destructOnTree(bestTree);
                stTree_destruct(bestTree);
            }
            bestTree = curTree;
            bestDups = dups;
            bestLosses = losses;
        } else {
            stPhylogenyInfo_destructOnTree(curTree);
            stTree_destruct(curTree);
        }
    }
    stList_destruct(stack);
    return bestTree;
}

static stSet *climb(stTree *childGene, stTree *childLCARecon, stTree *parentGene,
                    stTree *parentLCARecon, stHash *N) {
    stSet *childN = stHash_search(N, childGene);
    assert(childN != NULL);
    if (stSet_search(childN, parentLCARecon)) {
        return childN;
    }

    for (int64_t i = 0; i < stTree_getChildNumber(parentLCARecon); i++) {
        if (stSet_search(childN, stTree_getChild(parentLCARecon, i))) {
            return childN;
        }
    }

    stSetIterator *setIt = stSet_getIterator(childN);
    stTree *x = stSet_getNext(setIt);
    stSet_destructIterator(setIt);

    while (stTree_getParent(x) != parentLCARecon) {
        x = stTree_getParent(x);
    }

    stSet *ret = stSet_construct();
    stSet_insert(ret, x);
    stSet_destruct(childN);
    return ret;
}

static stTree *stPhylogeny_reconcileNonBinary_R(stTree *gene, stHash *leafToSpecies,
                                                stHash *N, bool relabelAncestors) {
    stTree *LCARecon;
    stReconciliationEvent event;
    if (stTree_getChildNumber(gene) == 0) {
        // Leaves are already reconciled.
        LCARecon = stHash_search(leafToSpecies, gene);
        assert(LCARecon != NULL);
        stSet *myN = stSet_construct();
        stSet_insert(myN, LCARecon);
        stHash_insert(N, gene, myN);
        event = LEAF;
    } else {
        // Internal node
        // Calculate the LCA mapping
        stTree *leftLCARecon = stPhylogeny_reconcileNonBinary_R(
            stTree_getChild(gene, 0), leafToSpecies, N, relabelAncestors);
        stTree *rightLCARecon = stPhylogeny_reconcileNonBinary_R(
            stTree_getChild(gene, 1), leafToSpecies, N, relabelAncestors);
        LCARecon = stTree_getMRCA(leftLCARecon, rightLCARecon);
        // Calculate if this is a required duplication. We don't
        // really care if it's a conditional duplication.
        stSet *leftN = climb(stTree_getChild(gene, 0), leftLCARecon,
                             gene, LCARecon, N);
        stSet *rightN = climb(stTree_getChild(gene, 1), rightLCARecon,
                             gene, LCARecon, N);
        stSet *myN = stSet_getUnion(leftN, rightN);
        stHash_insert(N, gene, myN);
        stSet *intersect = stSet_getIntersection(leftN, rightN);
        if (stSet_size(intersect) != 0) {
            // This node is a required duplication
            event = DUPLICATION;
        } else {
            // This event is a conditional duplication or a
            // speciation.
            event = SPECIATION;
        }
        stSet_destruct(leftN);
        stSet_destruct(rightN);
        stSet_destruct(intersect);
        if (stTree_getParent(gene) == NULL) {
            stSet_destruct(myN);
        }
    }
    fillInReconciliationInfo(gene, LCARecon, event, relabelAncestors);
    return LCARecon;
}

void stPhylogeny_reconcileNonBinary(stTree *geneTree, stHash *leafToSpecies, bool relabelAncestors) {
    // TODO: this hash is likely unnecessary and values could probably
    // be passed up along the tree by stPhylogeny_reconcile_R.
    stHash *N = stHash_construct();
    stPhylogeny_reconcileNonBinary_R(geneTree, leafToSpecies, N, relabelAncestors);
    stHash_destruct(N);
}

void stPhylogeny_nni(stTree *anc, stTree **tree1, stTree **tree2) {
    if (stTree_getChildNumber(anc) == 0 || stTree_getParent(anc) == NULL) {
        // The branch we will be NNI'ing (the branch above anc) isn't
        // an internal branch.
        *tree1 = NULL;
        *tree2 = NULL;
        return;
    }

    // Get the root of this tree, and a trail of breadcrumbs back to
    // the position of this node (child of the branch to be NNI'd).
    int64_t distToRoot = 0;
    stTree *root = anc;
    while (stTree_getParent(root) != NULL) {
        root = stTree_getParent(root);
        distToRoot++;
    }

    bool pathFromRoot[distToRoot];
    int64_t i = distToRoot - 1;
    root = anc;
    stTree *prev;
    while (stTree_getParent(root) != NULL) {
        prev = root;
        root = stTree_getParent(root);
        if (prev == stTree_getChild(root, 0)) {
            pathFromRoot[i--] = 0;
        } else {
            assert(prev == stTree_getChild(root, 1));
            pathFromRoot[i--] = 1;
        }
    }    

    *tree1 = stTree_clone(root);
    *tree2 = stTree_clone(root);

    // Traverse down to the right node in tree1 and tree2.
    stTree *anc1 = *tree1;
    for (i = 0; i < distToRoot; i++) {
        if (pathFromRoot[i]) {
            anc1 = stTree_getChild(anc1, 1);
        } else {
            anc1 = stTree_getChild(anc1, 0);
        }
    }

    stTree *anc2 = *tree2;
    for (i = 0; i < distToRoot; i++) {
        if (pathFromRoot[i]) {
            anc2 = stTree_getChild(anc2, 1);
        } else {
            anc2 = stTree_getChild(anc2, 0);
        }
    }

    if (stTree_getParent(stTree_getParent(anc)) != NULL) {
        /*
         * Not a branch off the root node.
         * Original tree:
         *     \
         *     /\
         *    4 /\
         *     3 /\
         *      1  2
         * Tree 1:
         *     \
         *     /\
         *    2 /\
         *     3 /\
         *      1  4
         * Tree 2:
         *     \
         *     /\
         *    4 /\
         *     2 /\
         *      1  3
         */
        // Swap 2 and 4 in tree 1.
        stTree *two = stTree_getChild(anc1, 1);
        stTree *four;
        if (stTree_getChild(stTree_getParent(stTree_getParent(anc1)), 0) == stTree_getParent(anc1)) {
            four = stTree_getChild(stTree_getParent(stTree_getParent(anc1)), 1);
        } else {
            assert(stTree_getChild(stTree_getParent(stTree_getParent(anc1)), 1) == stTree_getParent(anc1));
            four = stTree_getChild(stTree_getParent(stTree_getParent(anc1)), 0);
        }
        stTree_setParent(two, stTree_getParent(four));
        stTree_setParent(four, anc1);
        // Swap 2 and 3 in tree 2.
        two = stTree_getChild(anc2, 1);
        stTree *three;
        if (stTree_getChild(stTree_getParent(anc2), 0) == anc2) {
            three = stTree_getChild(stTree_getParent(anc2), 1);
        } else {
            assert(stTree_getChild(stTree_getParent(anc2), 1) == anc2);
            three = stTree_getChild(stTree_getParent(anc2), 0);
        }
        stTree_setParent(two, stTree_getParent(three));
        stTree_setParent(three, anc2);
    } else {
        /*
         * A branch off the root node.
         * Original tree:
         *      /\
         *     /  \
         *    /    \
         *   /\    /\
         *  1  2  3  4
         * Tree 1:
         *      /\
         *     /  \
         *    /    \
         *   /\    /\
         *  1  4  3  2
         * Tree 2:
         *      /\
         *     /  \
         *    /    \
         *   /\    /\
         *  1  3  2  4
         */
        stTree *two = stTree_getChild(anc1, 1);
        stTree *three, *four;
        if (stTree_getChild(root, 0) == anc) {
            three = stTree_getChild(stTree_getChild(*tree2, 1), 0);
            four = stTree_getChild(stTree_getChild(*tree1, 1), 1);
        } else {
            assert(stTree_getChild(root, 1) == anc);
            three = stTree_getChild(stTree_getChild(*tree2, 0), 0);
            four = stTree_getChild(stTree_getChild(*tree1, 0), 1);
        }
        // Swap 2 and 4 in tree 1.
        stTree_setParent(two, stTree_getParent(four));
        stTree_setParent(four, anc1);
        // Swap 2 and 3 in tree 2.
        two = stTree_getChild(anc2, 1);
        stTree_setParent(two, stTree_getParent(three));
        stTree_setParent(three, anc2);
    }
}

// Determines whether a split satisfies the four-point criterion of
// Bandelt and Dress 1992. The "relaxed" parameter, if true, uses the
// condition stated in the paper (where the intra-split distance must
// not be larger than *both* inter-split distances), but if false,
// uses a stricter condition (that the intra-split distance must be
// smaller than *both* inter-split distances).
static bool satisfiesFourPoint(stMatrix *distanceMatrix, stList *leftSplitIndices, stList *rightSplitIndices, bool relaxed) {
    // This is a bit convoluted, but generates all possible
    // unordered combinations of indices i, j in the left side of the split. i,j
    // are distance matrix indices, not indices in the split list!
    for (int64_t left_i = 0; left_i < stList_length(leftSplitIndices); left_i++) {
        for (int64_t left_j = left_i + 1; left_j < stList_length(leftSplitIndices); left_j++) {
            int64_t i = stIntTuple_get(stList_get(leftSplitIndices, left_i), 0);
            int64_t j = stIntTuple_get(stList_get(leftSplitIndices, left_j), 0);
            // Similarly, generate all possible unordered combinations
            // k, l from the right side of the split.
            for (int64_t right_i = 0; right_i < stList_length(rightSplitIndices); right_i++) {
                for (int64_t right_j = right_i + 1; right_j < stList_length(rightSplitIndices); right_j++) {
                    int64_t k = stIntTuple_get(stList_get(rightSplitIndices, right_i), 0);
                    int64_t l = stIntTuple_get(stList_get(rightSplitIndices, right_j), 0);
                    // Do the check.
                    double intra = *stMatrix_getCell(distanceMatrix, i, j) + *stMatrix_getCell(distanceMatrix, k, l);
                    double inter1 = *stMatrix_getCell(distanceMatrix, i, k) + *stMatrix_getCell(distanceMatrix, j, l);
                    double inter2 = *stMatrix_getCell(distanceMatrix, i, l) + *stMatrix_getCell(distanceMatrix, j, k);
                    if (relaxed) {
                        if (intra >= inter1 && intra >= inter2) {
                            return false;
                        }
                    } else {
                        if (intra >= inter1 || intra >= inter2) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    // If we're here, then we've checked all the possible quartets without failing.
    return true;
}

static stSplit *stSplit_construct(stList *leftSplit, stList *rightSplit, double isolationIndex) {
    stSplit *ret = st_malloc(sizeof(stSplit));
    ret->leftSplit = leftSplit;
    ret->rightSplit = rightSplit;
    ret->isolationIndex = isolationIndex;
    return ret;
}

static void stSplit_destruct(stSplit *split) {
    stList_destruct(split->leftSplit);
    stList_destruct(split->rightSplit);
    free(split);
}

// Compare two d-splits by their isolation indexes.
static int stSplit_cmp(stSplit *split1, stSplit *split2) {
    if (split1->isolationIndex < split2->isolationIndex) {
        return -1;
    } else if (split1->isolationIndex > split2->isolationIndex) {
        return 1;
    } else {
        return 0;
    }
}

static void assignIsolationIndex(stMatrix *distanceMatrix, stSplit *split) {
    // We want to find the minimum of (maximum of inter-split distances - intra-split distance) / 2
    // from all cross-split quartets.
    double min_isolation = DBL_MAX;
    for (int64_t left_i = 0; left_i < stList_length(split->leftSplit); left_i++) {
        for (int64_t left_j = left_i + 1; left_j < stList_length(split->leftSplit); left_j++) {
            int64_t i = stIntTuple_get(stList_get(split->leftSplit, left_i), 0);
            int64_t j = stIntTuple_get(stList_get(split->leftSplit, left_j), 0);

            for (int64_t right_i = 0; right_i < stList_length(split->rightSplit); right_i++) {
                for (int64_t right_j = right_i + 1; right_j < stList_length(split->rightSplit); right_j++) {
                    int64_t k = stIntTuple_get(stList_get(split->rightSplit, right_i), 0);
                    int64_t l = stIntTuple_get(stList_get(split->rightSplit, right_j), 0);
                    double intra = *stMatrix_getCell(distanceMatrix, i, j) + *stMatrix_getCell(distanceMatrix, k, l);
                    double inter1 = *stMatrix_getCell(distanceMatrix, i, k) + *stMatrix_getCell(distanceMatrix, j, l);
                    double inter2 = *stMatrix_getCell(distanceMatrix, i, l) + *stMatrix_getCell(distanceMatrix, j, k);
                    double max_dist = intra;
                    if (inter1 > max_dist) {
                        max_dist = inter1;
                    }
                    if (inter2 > max_dist) {
                        max_dist = inter2;
                    }
                    max_dist -= intra;
                    if (max_dist < min_isolation) {
                        min_isolation = max_dist;
                    }
                }
            }
        }
    }
    split->isolationIndex = min_isolation / 2;
}

stList *stPhylogeny_getSplits(stMatrix *distanceMatrix, bool relaxed) {
    assert(stMatrix_m(distanceMatrix) == stMatrix_n(distanceMatrix));
    stList *splits = stList_construct3(0, (void (*)(void *)) stSplit_destruct);
    for (int64_t i = 1; i < stMatrix_m(distanceMatrix); i++) {
        stList *singletonSplitLeft = stList_construct3(0, free);
        stList_append(singletonSplitLeft, stIntTuple_construct1(i));
        stList *singletonSplitRight = stList_construct3(0, free);
        for (int64_t j = 0; j < i; j++) {
            stList_append(singletonSplitRight, stIntTuple_construct1(j));
        }
        stList *newSplits = stList_construct3(0, (void (*)(void *)) stSplit_destruct);
        stList_append(newSplits, stSplit_construct(singletonSplitLeft, singletonSplitRight, 0.0));
        while (stList_length(splits) > 0) {
            stSplit *split = stList_pop(splits);
            stIntTuple *iTuple = stIntTuple_construct1(i);
            stList_append(split->leftSplit, iTuple);
            bool addToLeft = satisfiesFourPoint(distanceMatrix, split->leftSplit, split->rightSplit, relaxed);
            stList_pop(split->leftSplit);
            stList_append(split->rightSplit, iTuple);
            bool addToRight = satisfiesFourPoint(distanceMatrix, split->leftSplit, split->rightSplit, relaxed);
            stList_pop(split->rightSplit);
            if (addToRight && addToLeft) {
                // We are making two new splits, so have to clone the
                // lists and their elements. For no particular reason,
                // the cloned one becomes the one with i added to the
                // right.

                // First we handle the split with i added to its right
                // split, cloning the existing list, and add that to
                // the new split list.
                stList *addedToRight_rightSplit = stList_construct3(0, free);
                for (int64_t j = 0; j < stList_length(split->rightSplit); j++) {
                    stIntTuple *new = stIntTuple_construct1(stIntTuple_get(stList_get(split->rightSplit, j), 0));
                    stList_append(addedToRight_rightSplit, new);
                }
                stList_append(addedToRight_rightSplit, stIntTuple_construct1(i));
                stList *addedToRight_leftSplit = stList_construct3(0, free);
                for (int64_t j = 0; j < stList_length(split->leftSplit); j++) {
                    stIntTuple *new = stIntTuple_construct1(stIntTuple_get(stList_get(split->leftSplit, j), 0));
                    stList_append(addedToRight_leftSplit, new);
                }
                stSplit *addedToRight = stSplit_construct(addedToRight_leftSplit, addedToRight_rightSplit, 0.0);
                stList_append(newSplits, addedToRight);

                // Now add i to the left split of the existing split
                // and add that to the new split list.
                stList_append(split->leftSplit, iTuple);
                stList_append(newSplits, split);
            } else if (addToRight) {
                stList_append(split->rightSplit, iTuple);
                stList_append(newSplits, split);
            } else if (addToLeft) {
                stList_append(split->leftSplit, iTuple);
                stList_append(newSplits, split);
            } else {
                stSplit_destruct(split);
            }
        }
        stList_destruct(splits);
        splits = newSplits;
    }

    // Remove the remaining trivial splits and assign isolation indexes.
    for (int64_t i = 0; i < stList_length(splits); i++) {
        stSplit *split = stList_get(splits, i);
        assignIsolationIndex(distanceMatrix, split);
        if (stList_length(split->leftSplit) == 1 || stList_length(split->rightSplit) == 1) {
            // inefficient
            stList_remove(splits, i);
            stSplit_destruct(split);
            // adjusting index to account for removed item
            i--;
        }
    }

    // Sort by isolation index in descending order.
    stList_sort(splits, (int (*)(const void *, const void *)) stSplit_cmp);
    stList_reverse(splits);
    return splits;
}

static bool isCompatibleSplit(stList *splitIndices, stHash *indexToLeaf) {
    stTree *parent = stTree_getParent(stHash_search(indexToLeaf, stList_get(splitIndices, 0)));
    assert(parent != NULL);
    for (int64_t i = 1; i < stList_length(splitIndices); i++) {
        stTree *leaf = stHash_search(indexToLeaf, stList_get(splitIndices, i));
        if (stTree_getParent(leaf) != parent) {
            return false;
        }
    }
    return true;
}

static void applyCompatibleSplit(stList *splitIndices, stHash *indexToLeaf) {
    stTree *parent = stTree_getParent(stHash_search(indexToLeaf, stList_get(splitIndices, 0)));
    stTree *newNode = stTree_construct();
    stTree_setParent(newNode, parent);
    // Branch lengths are arbitrarily set to 1.0.
    stTree_setBranchLength(newNode, 1.0);
    for (int64_t i = 0; i < stList_length(splitIndices); i++) {
        stTree *leaf = stHash_search(indexToLeaf, stList_get(splitIndices, i));
        stTree_setParent(leaf, newNode);
    }
}

stTree *stPhylogeny_greedySplitDecomposition(stMatrix *distanceMatrix, bool relaxed) {
    assert(stMatrix_m(distanceMatrix) == stMatrix_n(distanceMatrix));
    stHash *indexToLeaf = stHash_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey, (int (*)(const void *, const void *)) stIntTuple_equalsFn, (void (*)(void *)) stIntTuple_destruct, NULL);
    // We start out with a complete star phylogeny.
    stTree *root = stTree_construct();
    for (int64_t i = 0; i < stMatrix_m(distanceMatrix); i++) {
        stTree *leaf = stTree_construct();
        stHash_insert(indexToLeaf, stIntTuple_construct1(i), leaf);
        char *label = stString_print_r("%" PRIi64, i);
        stTree_setLabel(leaf, label);
        free(label);
        stTree_setParent(leaf, root);
        // The branch lengths are all arbitrarily set to 1.0 to avoid
        // infinite branch lengths setting up a minefield for any
        // arithmetic later on.
        stTree_setBranchLength(leaf, 1.0);
    }

    stList *splits = stPhylogeny_getSplits(distanceMatrix, relaxed);
    // Start adding compatible splits to the tree, creating a new
    // internal node for each split which groups together one of its
    // sides.
    for (int64_t i = 0; i < stList_length(splits); i++) {
        stSplit *split = stList_get(splits, i);
        if (isCompatibleSplit(split->leftSplit, indexToLeaf)) {
            applyCompatibleSplit(split->leftSplit, indexToLeaf);
        } else if (isCompatibleSplit(split->rightSplit, indexToLeaf)) {
            applyCompatibleSplit(split->rightSplit, indexToLeaf);
        }
    }
    stList_destruct(splits);
    stHash_destruct(indexToLeaf);
    stPhylogeny_addStIndexedTreeInfo(root);
    return root;
}

void stPhylogeny_applyJukesCantorCorrection(stMatrix *distanceMatrix) {
    for (int64_t i = 0; i < stMatrix_m(distanceMatrix); i++) {
        for (int64_t j = 0; j < stMatrix_n(distanceMatrix); j++) {
            if (*stMatrix_getCell(distanceMatrix, i, j) < 0.75) {
                *stMatrix_getCell(distanceMatrix, i, j) = -0.75 * log(1 - 4 * (*stMatrix_getCell(distanceMatrix, i, j)) / 3);
            } else {
                // Having <25% identity isn't valid under the JC
                // model, so we just set the distance to something
                // higher than any realistic distance (not infinity as
                // that may break some arithmetic down the road).
                *stMatrix_getCell(distanceMatrix, i, j) = 10000.0;
            }
        }
    }
}
