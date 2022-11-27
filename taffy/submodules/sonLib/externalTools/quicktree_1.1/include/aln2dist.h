/*  Last edited: Aug 24 14:46 1999 (klh) */
/***********************************************************************
 ** FILE: aln2dist.h
 ** NOTES:
 **  This program produces a distance matrix in phylip format for the 
 **  given alignment file. The optional -k flag is used to determine 
 **  whether Kimura's tranlation is used to hopefully obtain more 
 **  accurate distances. The matrix is written to stdout
 ***********************************************************************/

#ifndef _ALN2DIST
#define _ALN2DIST

#include "options.h"
#include "util.h"
#include "align.h"
#include "cluster.h"


/********************** function prototypes ***************************/


/*********************************************************************
 FUNCTION: generate_matrix_matrixmaker
 DESCRIPTION: 
   Given a fule handle to an alignment, creates and prints a 
   distance matrix for the sequences in the alignment in phylip format
   to stdout
 RETURNS: void
 ARGS: 
   An open alignment FILE
 NOTES: 
 *********************************************************************/
void generate_DistanceMatrix_matrixmaker( FILE * );

#endif
