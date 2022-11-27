/*  Last edited: Jun 21 11:57 2001 (klh) */
/***********************************************************************
 ** FILE: aln2tree.c
 ** NOTES:
 **  This program constructs a phylogenetic tree from a given align-
 **  ment file. The process comprises of two phases: 
 **    1. Construction of a distance matrix. Here, the optional -k
 **       is used to determine whether Kimura's tranlation is used to
 **       hopefully obtain more accurate distance
 **    2. The construction of the tree. It is printed in New-Hampshire
 **       format to standard output
 ***********************************************************************/

#include "aln2tree.h"


/********************* command-line options * ************************/

static char usage[] = "\
Usage: treemaker [-options] <alignfile>\n\
Options are:\n\
-kimura          : Use the kimura translation for pairwise distances\n\
-boot <n>        : Calcuate bootstrap values with n iterations\n\
";

static struct Option options[] = {
    { "-kimura", NO_ARGS },
    { "-boot",   INT_ARG },
    { "-h",      NO_ARGS }
};

static unsigned int use_kimura = 0; 
static unsigned int calc_bootstraps = 0;


/*********************************************************************
 FUNCTION: maketree_treemaker
 DESCRIPTION: 
   Given a fule handle to an alignment, creates and prints a neighbour
   joining phlogenetic tree for the sequences in the alignment
 RETURNS: void
 ARGS: 
   void
 NOTES: 
 *********************************************************************/

void maketree_treemaker( FILE *input ) {
  struct Alignment *aln, *cons_aln;
  struct Tree *myTree, *testTree;
  struct ClusterGroup *group;
  unsigned int trial; 

  /******* METHOD ***************************************
    1. Read the alignment.
    2. Create a ClusterGroup from the alignment,
         merging all identical sequences
    3. Grab the alignment of consensi from the 
         ClusterGroup for later bootstrapping
    4. Build the tree from the cluster group
    5. Bootstrap if required
    6. Print the tree
  *******************************************************/

  if ( ( aln = read_Stockholm_Alignment( input )) == NULL)
    fatal_util( "Parse error in alignment file: Must be in MUL format");
  fclose( input );

  /*
  write_MUL_Alignment(stderr, aln);
  exit(0);
  */

  group = alignment_to_ClusterGroup( aln, TRUE );
  aln = free_Alignment( aln );
  cons_aln = consensus_aln_from_ClusterGroup( group );
  group->matrix = empty_DistanceMatrix( group->numclusters );
  calc_DistanceMatrix( group->matrix, cons_aln, FALSE, use_kimura );

  myTree = neighbour_joining_buildtree( group, calc_bootstraps );

  if (calc_bootstraps) {
    for (trial=0; trial < calc_bootstraps; trial++) {
      calc_DistanceMatrix( group->matrix, cons_aln, TRUE, use_kimura );

      testTree = neighbour_joining_buildtree( group, calc_bootstraps );
      update_bootstraps_Tree( myTree, testTree, group->numclusters );
      testTree = free_Tree( testTree );
    }
    scale_bootstraps_Tree( myTree, calc_bootstraps );
    write_newhampshire_Tree( stdout, myTree, TRUE );
  }
  else {
    write_newhampshire_Tree( stdout, myTree, FALSE );
  }

  cons_aln = free_Alignment( cons_aln );
  group = free_ClusterGroup( group );
  myTree = free_Tree( myTree );

  exit(0);
}





/*********************************************************************
		        MAIN
 *********************************************************************/

void main (int argc, char *argv[]) {

  FILE *input;
  char *fname;
  
  char *optname;                /* name of option found by getoption */
  char *optarg;                 /* argument found by getoption       */
  unsigned int optindex;        /* index in argv[]             */


  while (get_option( argc, 
		     argv, 
		     options,
		     sizeof(options) / sizeof( struct Option ),
		     usage,
		     &optindex,
		     &optname,
		     &optarg )){
    if (strcmp(optname, "-kimura") == 0) use_kimura = 1;
    else if (strcmp(optname, "-boot") == 0) calc_bootstraps = atoi( optarg );
    else if (strcmp(optname, "-h") == 0) {
      fprintf( stderr, "%s", usage );
      exit(0);
    }
  }

  if (argc - optindex != 1) {
    fatal_util("Fatal error: Incorrect number of arguements.\n%s\n", usage );
  }

  fname = argv[optindex];

  if ((input = fopen( fname, "r" )) == NULL)
    fatal_util( "Could not open file %s for reading", fname );
  
  maketree_treemaker( input );
}
