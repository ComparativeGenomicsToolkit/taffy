/*
 * Copyright (C) 2006-2014 by Joel Armstrong & Benedict Paten
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef st_phylogeny_h_
#define st_phylogeny_h_

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

// Used for identifying the event underlying a node in a reconciled tree.
typedef enum {
    DUPLICATION,
    SPECIATION,
    LEAF
} stReconciliationEvent;

// Reconciliation information. Only defined on internal nodes.
typedef struct {
    stTree *species; // The node in the species tree that this node maps to.
    stReconciliationEvent event; // Duplication or speciation node.
} stReconciliationInfo;

// Data structure for storing information about a node in a
// neighbor-joined tree.
typedef struct {
    int64_t matrixIndex;         // = -1 if an internal node, index
                                 // into distance matrix if a leaf.
    char *leavesBelow;           // leavesBelow[i] = 1 if leaf i is
                                 // present below this node, 0
                                 // otherwise. Could be a bit array,
                                 // which would make things much
                                 // faster.
    int64_t numBootstraps;       // Number of bootstraps that support
                                 // this split.
    double bootstrapSupport;     // Fraction of bootstraps that
                                 // support this split.
    int64_t totalNumLeaves;      // Total number of leaves overall in
                                 // the tree (which is the size of
                                 // leavesBelow). Not strictly
                                 // necessary, but convenient
} stIndexedTreeInfo;

// Used as client-data on nodes in trees created/used by the
// functions in this file.
typedef struct {
    stReconciliationInfo *recon; // Reconciliation info. Can be NULL.
    stIndexedTreeInfo *index;    // Index, in case this tree has a
                                 // previous numbering of its
                                 // leaves. Allows more efficient traversals.
                                 // Can be NULL.
} stPhylogenyInfo;

// Represents a d-split as characterized by Bandelt and Dress, 1992.
typedef struct {
    stList *leftSplit;
    stList *rightSplit;
    double isolationIndex;
} stSplit;

// Free a stPhylogenyInfo struct.
void stPhylogenyInfo_destruct(stPhylogenyInfo *info);

// Free stPhylogenyInfo data on all nodes in the tree.
void stPhylogenyInfo_destructOnTree(stTree *tree);

// Clones stPhylogenyInfo and all sub-info.
stPhylogenyInfo *stPhylogenyInfo_clone(stPhylogenyInfo *info);

// Add valid stIndexedTreeInfo (as a subfield of stPhylogenyInfo) to a
// tree with leaves labeled 0, 1, 2, 3, etc. and with internal nodes
// unlabeled. This function will fail on trees labeled any other way.
void stPhylogeny_addStIndexedTreeInfo(stTree *tree);

// Set (and allocate) the leavesBelow and totalNumLeaves attribute in
// the stIndexedTreeInfo for the given tree and all subtrees.
// The stIndexedTreeInfo must already be allocated!
void stPhylogeny_setLeavesBelow(stTree *tree, int64_t totalNumLeaves);

// Return a new tree which has its partitions scored by how often they
// appear in the bootstrap. This fills in the numBootstraps and
// bootstrapSupport fields of each node. All trees must have valid
// stPhylogenyInfo.
stTree *stPhylogeny_scoreFromBootstrap(stTree *tree, stTree *bootstrap);

// Return a new tree which has its partitions scored by how often they
// appear in the bootstraps. This fills in the numBootstraps and
// bootstrapSupport fields of each node. All trees must have valid
// stPhylogenyInfo.
stTree *stPhylogeny_scoreFromBootstraps(stTree *tree, stList *bootstraps);

stTree *stPhylogeny_scoreReconciliationFromBootstrap(stTree *tree,
                                                     stTree *bootstrap);

// Return a new tree which has its partitions scored by how often they
// appear with the same reconciliation in the bootstraps. Obviously,
// this is always less than or equal to the normal bootstrap
// score. This fills in the numBootstraps and bootstrapSupport fields
// of each node. All trees must have valid stPhylogenyInfo.
stTree *stPhylogeny_scoreReconciliationFromBootstraps(stTree *tree,
                                                      stList *bootstraps);

// Only one half of the distanceMatrix is used, distances[i][j] for which i > j
// Tree returned is labeled by the indices of the distance matrix. The
// tree is rooted halfway along the longest branch if outgroups is
// NULL; otherwise, it's rooted halfway along the longest outgroup
// branch.
stTree *stPhylogeny_neighborJoin(stMatrix *distances, stList *outgroups);

// Gets the (leaf) node corresponding to an index in the distance matrix.
// Requires an indexed tree (which has stPhylogenyInfo with non-null
// stIndexedTreeInfo.)
stTree *stPhylogeny_getLeafByIndex(stTree *tree, int64_t leafIndex);

// Find the distance between two arbitrary nodes (which must be in the
// same tree).
// Requires an indexed tree (which has stPhylogenyInfo with non-null
// stIndexedTreeInfo.)
double stPhylogeny_distanceBetweenNodes(stTree *node1, stTree *node2);

// Find the distance between leaves (given by their index in the
// distance matrix.)
// Requires an indexed tree (which has stPhylogenyInfo with non-null
// stIndexedTreeInfo.)
double stPhylogeny_distanceBetweenLeaves(stTree *tree, int64_t leaf1,
                                         int64_t leaf2);

// Return the MRCA of the given leaves. More efficient than
// stTree_getMRCA.
// Requires an indexed tree (which has stPhylogenyInfo with non-null
// stIndexedTreeInfo.)
stTree *stPhylogeny_getMRCA(stTree *tree, int64_t leaf1, int64_t leaf2);

// Compute join costs for a species tree for use in guided
// neighbor-joining. These costs are calculated by penalizing
// according to the number of dups and losses implied by the
// reconciliation when joining two genes reconciled to certain nodes
// of the species tree.
// speciesToIndex (a blank hash) will be populated with stIntTuples
// corresponding to each species' index into the join cost matrix.
// NB: the species tree must be binary.
stMatrix *stPhylogeny_computeJoinCosts(stTree *speciesTree,
                                       stHash *speciesToIndex,
                                       double costPerDup, double costPerLoss);

// Precompute an MRCA matrix to pass to guided neighbor-joining.
int64_t **stPhylogeny_getMRCAMatrix(stTree *speciesTree, stHash *speciesToIndex);

// Neighbor joining guided by a species tree. The similarity matrix
// has entries for each i and j such that s_ij if i > j is #
// differences between i and j, i < j is # similarities between i and
// j. Join costs should be precomputed by
// stPhylogeny_computeJoinCosts. indexToSpecies is a map from matrix
// index (of the similarity matrix) to species leaves.
stTree *stPhylogeny_guidedNeighborJoining(stMatrix *distanceMatrix,
                                          stMatrix *similarityMatrix,
                                          stMatrix *joinCosts,
                                          stHash *matrixIndexToJoinCostIndex,
                                          stHash *speciesToJoinCostIndex,
                                          int64_t **speciesMRCAMatrix,
                                          stTree *speciesTree);

// Reconcile a gene tree (without rerooting), set the proper
// stReconcilationInfo (as an entry of stPhylogenyInfo) as client data
// on all nodes, and optionally set the labels of the ancestors to the
// labels of their reconciliation in the species tree.
//
// The species tree must be "at-most-binary", i.e. it must have no
// nodes with more than 3 children, but may have nodes with only one
// child.
void stPhylogeny_reconcileAtMostBinary(stTree *geneTree, stHash *leafToSpecies,
                                       bool relabelAncestors);

// For a tree that has already been reconciled by
// reconcileAtMostBinary, calculates the number of dups and losses
// implied by the reconciliation. dups and losses must be set to 0
// before calling. The gene-tree may have polytomies.
void stPhylogeny_reconciliationCostAtMostBinary(stTree *reconciledTree,
                                                int64_t *dups,
                                                int64_t *losses);

// Return a copy of geneTree that is rooted to minimize duplications.
// The gene tree must not have any polytomies.
// NOTE: the returned tree does *not* have reconciliation info set,
// and this function will reconcile geneTree, resetting any
// reconciliation information that potentially already exists.
stTree *stPhylogeny_rootByReconciliationAtMostBinary(stTree *geneTree,
                                                     stHash *leafToSpecies);

// Return a copy of geneTree that is rooted to minimize duplications.
// This function is slower than rootByReconciliationAtMostBinary, but
// supports gene trees with polytomies.
stTree *stPhylogeny_rootByReconciliationNaive(stTree *geneTree, stHash *leafToSpecies);

// Reconcile a binary gene tree to a species tree that may include
// polytomies. Based on the method used by NOTUNG, described in
// Vernot, Stolzer, Goldman, Durand, J Comput Biol 2008.
void stPhylogeny_reconcileNonBinary(stTree *geneTree, stHash *leafToSpecies,
                                    bool relabelAncestors);

// Nearest-neighbor interchange (for strictly binary trees). Returns
// two new neighboring trees in the parameters tree1 and tree2. Does
// not modify the original tree.
void stPhylogeny_nni(stTree *anc, stTree **tree1, stTree **tree2);

// Get the non-trivial d-splits as described in Bandelt and Dress,
// 1992. The "relaxed" parameter, if true, uses the condition stated
// in the paper (where the intra-split distance must not be larger
// than *both* inter-split distances), but if false, uses a stricter
// condition (that the intra-split distance must be smaller than
// *both* inter-split distances).
stList *stPhylogeny_getSplits(stMatrix *distanceMatrix, bool relaxed);

// Build a tree greedily using the d-splits from stPhylogeny_getSplits.
stTree *stPhylogeny_greedySplitDecomposition(stMatrix *distanceMatrix, bool relaxed);

// Apply the Jukes-Cantor distance correction to the input distance matrix.
void stPhylogeny_applyJukesCantorCorrection(stMatrix *distanceMatrix);

#ifdef __cplusplus
}
#endif
#endif
