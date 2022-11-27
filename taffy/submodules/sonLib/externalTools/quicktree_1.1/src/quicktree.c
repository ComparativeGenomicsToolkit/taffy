/*  Last edited: Feb  1 18:22 2002 (klh) */
/***********************************************************************
 ** FILE: quicktree.c
 ** NOTES:
 **  This program does one of three things:
 ** 1. Given an aligment, construct a distance matrix.
 ** 2. Given a distance matrix, build a tree
 ** 3. Given an aligment, build a tree
 ** Which of these modes is active is dependent upon the command options
 ** 
 ***********************************************************************/

#include "options.h"
#include "util.h"
#include "align.h"
#include "cluster.h"
#include "buildtree.h"
#include "distancemat.h"
#include "tree.h"


/********************* command-line options * ************************/

static char usage[] = "\
Usage: quicktree [-options] <inputfile>\n\
Options:\n\
-in <m|a>        : input file is a distance matrix in phylip format (m)\n\
                   or an alignment in stockholm format* (a, default)\n\
-out <m|t>       : output is a distance matrix in phylip format (m) or\n\
                   a tree in New Hampshire format\n\n\
Advanced options :\n\
-upgma           : Use the UPGMA method to construct the tree\n\
                     (ignored for distance matrix outputs)\n\
-kimura          : Use the kimura translation for pairwise distances\n\
                     (ignored for distance matrix inputs)\n\
-boot <n>        : Calcuate bootstrap values with n iterations\n\
                     (ignored for distance matrix outputs)\n\n\
*Use sreformat, part of the HMMer package to convert your alignment to Stockholm format\n";

static struct Option options[] = {
    { "-in", CHAR_ARG },
    { "-out", CHAR_ARG },
    { "-upgma", NO_ARGS },
    { "-kimura", NO_ARGS },
    { "-boot",   INT_ARG },
    { "-h",      NO_ARGS }
};

static unsigned int use_kimura = 0; 
static unsigned int calc_bootstraps = 0;
static unsigned int use_upgma = 0;
static unsigned int input_is_matrix = 0;
static unsigned int output_is_matrix = 0;


static void quicktree( FILE *input) {
  unsigned int trial;
  struct Alignment *aln, *cons_aln = NULL;
  struct DistanceMatrix *mat;
  struct ClusterGroup *group;
  struct Tree *myTree, *testTree;

  struct Tree *(*tree_func)(struct ClusterGroup *, unsigned int);
  /* step 1: get the distance matrix (either build or read it) */


  if (input_is_matrix) {

    mat = read_phylip_DistanceMatrix( input, &aln  );

    fclose( input );

    /* Note that the alignment here contains no sequences, just names,
       so there is no point in trying to merge identical sequences, hence 
       the FALSE argument 

    */

    group = alignment_to_ClusterGroup( aln, FALSE );
    group->matrix = mat;

  }
  else {
    /* input must be an alignment */

    if ( ( aln = read_Stockholm_Alignment( input )) == NULL)
      fatal_util( "Parse error in alignment file: Must be in MUL format");
    fclose( input );
    
    if (output_is_matrix) {
      mat = empty_DistanceMatrix( aln->numseqs );
      calc_DistanceMatrix( mat, aln, FALSE, use_kimura );
      write_phylip_DistanceMatrix( stdout, mat, aln );
      mat = free_DistanceMatrix( mat );
      return;
    }

    group = alignment_to_ClusterGroup( aln, TRUE );
    aln = free_Alignment( aln );
    cons_aln = consensus_aln_from_ClusterGroup( group );
    group->matrix = empty_DistanceMatrix( group->numclusters );
    calc_DistanceMatrix( group->matrix, cons_aln, FALSE, use_kimura );

  }

  /* step 2 produce tree */

  tree_func = use_upgma ? &UPGMA_buildtree : &neighbour_joining_buildtree;
  
  myTree = (*tree_func)( group, calc_bootstraps ); 
  
  if (calc_bootstraps) {
    for (trial=0; trial < calc_bootstraps; trial++) {
      calc_DistanceMatrix( group->matrix, cons_aln, TRUE, use_kimura );
      testTree = (*tree_func)( group, calc_bootstraps ); 
      update_bootstraps_Tree( myTree, testTree, group->numclusters );
      testTree = free_Tree( testTree );
    }
    scale_bootstraps_Tree( myTree, calc_bootstraps );
    cons_aln = free_Alignment( cons_aln );
  }
  
  write_newhampshire_Tree( stdout, myTree, calc_bootstraps );

  aln = free_Alignment( aln );
  group = free_ClusterGroup( group );
  myTree = free_Tree( myTree );

}



/*********************************************************************
		        MAIN
 *********************************************************************/

int main (int argc, char *argv[]) {

  FILE *input;
  char *fname;
  
  char *optname;                /* name of option found by getoption */
  char *optarg;                 /* argument found by getoption       */
  unsigned int optindex;        /* index in argv[]             */

  /**********************************/
  /** process command line options **/				    
  /**********************************/

  while (get_option( argc, 
		     argv, 
		     options,
		     sizeof(options) / sizeof( struct Option ),
		     usage,
		     &optindex,
		     &optname,
		     &optarg )){
    if (strcmp(optname, "-upgma") == 0) use_upgma = 1;
    else if (strcmp(optname, "-kimura") == 0) use_kimura = 1;
    else if (strcmp(optname, "-boot") == 0) calc_bootstraps = atoi( optarg );
    else if (strcmp(optname, "-in") == 0) {
      if (optarg[0] != 'm' && optarg[0] != 'a')
	fatal_util("Error: Incorrect use of '-in' optio\n%s\n", usage );
      
      if (optarg[0] == 'm')
	input_is_matrix = 1;
      else if (optarg[0] == 'a')
	input_is_matrix = 0;
    }
    else if (strcmp(optname, "-out") == 0) {
      if (optarg[0] != 'm' && optarg[0] != 't')
	fatal_util("Error: Incorrect use of '-out' option\n%s\n", usage );
      
      if (optarg[0] == 'm')
	output_is_matrix = 1;
      else if (optarg[0] == 't')
	output_is_matrix = 0;

    } 
    else if (strcmp(optname, "-h") == 0) {
      fprintf( stderr, "%s", usage );
      exit(0);
    }
  }

  if (argc - optindex != 1) {
    fatal_util("Fatal error: Incorrect number of arguements.\n%s\n", usage );
  }

  fname = argv[optindex];

  /* do some options combination checking - some combinations do not make sense */

  if (input_is_matrix && use_kimura)
    fatal_util("You can only use kimura when giving an aligment\n%s\n", usage );

  if ( (input_is_matrix || output_is_matrix) && calc_bootstraps)
    fatal_util("You can't do bootstrapping in matrix mode. You need an alignment\n%s\n", usage );


  if ((input = fopen( fname, "r" )) == NULL)
    fatal_util( "Could not open file %s for reading", fname );

  quicktree( input );

  exit(0);
}
