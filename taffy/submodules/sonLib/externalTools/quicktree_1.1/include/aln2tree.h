/*  Last edited: Aug 24 14:45 1999 (klh) */
/***********************************************************************
 ** FILE: aln2tree.h
 ** NOTES:
 **  This program constructs a phylogenetic tree from a given align-
 **  ment file. The process comprises of two phases: 
 **    1. Construction of a distance matrix. Here, the optional -k
 **       is used to determine whether Kimura's tranlation is used to
 **       hopefully obtain more accurate distance
 **    2. The construction of the tree. It is printed in New-Hampshire
 **       format to standard output
 ***********************************************************************/

#ifndef _ALN2TREE
#define _ALN2TREE

#include "options.h"
#include "util.h"
#include "align.h"
#include "cluster.h"
#include "buildtree.h"

/********************** function prototypes ***************************/


/*********************************************************************
 FUNCTION: maketree_treemaker
 DESCRIPTION: 
   Given a fule handle to an alignment, creates and prints a neighbour
   joining phlogenetic tree for the sequences in the alignment
 RETURNS: void
 ARGS: 
   A file handle to an alignment (should be in MUL format);
 NOTES: 
 *********************************************************************/
void maketree_treemaker( FILE * );


#endif
