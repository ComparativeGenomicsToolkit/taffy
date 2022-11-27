/*  Last edited: Aug 25 15:29 1999 (klh) */
/***********************************************************************
 ** FILE: aln2dist.c
 ** NOTES:
 **  This program constructs a DistanceMatrix from the given 
 **  multiple sequence alignment
 ***********************************************************************/

#include "aln2dist.h"


/********************* command-line options * ************************/

static char usage[] = "\
Usage: treemaker [-options] <alignfile>\n\
Options are:\n\
-kimura          : Use the kimura translation for pairwise distances\n\
";

static struct Option options[] = {
    { "-kimura", NO_ARGS },
    { "-h",      NO_ARGS }
};

static unsigned int use_kimura = 0; 


/*********************************************************************
 FUNCTION: produce_matrix_aln2dist
 DESCRIPTION: 
   Given a fule handle to an alignment, creates and prints a neighbour
   joining phlogenetic tree for the sequences in the alignment
 RETURNS: void
 ARGS: 
   void
 NOTES: 
 *********************************************************************/

void maketree_treemaker( FILE *input ) {
  struct Alignment *aln;
  struct DistanceMatrix *mat;

  /******* METHOD ***************************************
    1. Read the alignment.
    2. Use it to create the distance matrix
    3. Print the matrix in phylip format
  *******************************************************/

  if ( ( aln = read_MUL_Alignment( input )) == NULL)
    fatal_util( "Parse error in alignment file: Must be in MUL format");
  fclose( input );

  mat = empty_DistanceMatrix( aln->numseqs );
  calc_DistanceMatrix( mat, aln, FALSE, use_kimura );
  write_phylip_DistanceMatrix( stdout, mat, aln );

  /* before freeing the cluster, nullify the members field, which is actually
     pointing to part of the alignment */

  mat = free_DistanceMatrix( mat ); 
  aln = free_Alignment( aln );

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
