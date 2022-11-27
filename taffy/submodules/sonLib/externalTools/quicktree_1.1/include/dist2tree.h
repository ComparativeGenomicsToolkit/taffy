/*  Last edited: Aug 25 15:33 1999 (klh) */
/**********************************************************************
 ** FILE: dist2tree.h
 ** NOTES:
 ** This program can be used to generate neighbour-joining
 ** trees given a distance matrix file.
 ** The matrix file should be in phylip format; the application was
 ** originally coded to compare the efficiency of the writer's 
 ** implementation of neighbourjoining with that of phylip's 
 ** 'neighbor' program. The tree is writtento stdout in New Hants 
 ** format.
 **********************************************************************/


#ifndef _DIST2TREE
#define _DIST2TREE

#include "align.h"
#include "options.h"
#include "util.h"
#include "cluster.h"
#include "buildtree.h"


/********************** function prototypes ***************************/

/*********************************************************************
 FUNCTION: maketree_neighbourjoin
 DESCRIPTION: 
   Given a file handle to a distance matrix, creates and prints a neighbour
   joining phlogenetic tree for the sequences in the alignment
 RETURNS: void
 ARGS: 
   A file handle to a distance matrix (in 'phylip' format)
 NOTES: 
 *********************************************************************/
void maketree_neighbourjoin( FILE * );


#endif
