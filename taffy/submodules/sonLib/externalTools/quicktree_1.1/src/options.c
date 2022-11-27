/*  Last edited: Feb  1 15:45 2002 (klh) */
/**********************************************************************
 ** FILE: options.c
 ** DESCRIPTION:
 **  Rudimentary provision of command-line options
 **********************************************************************/

#include "options.h"



/*********************************************************************
 FUNCTION: get_option
 DESCRIPTION: 
   Gets an option from the given command line
 RETURNS:
   1, if a valid option was found
   0, if no valid option was found and option parsing is therefore 
      complete
 ARGS: 
 NOTES: 
*********************************************************************/
unsigned int get_option(int argc, 
			char **argv, 
			struct Option *opt, 
			unsigned int num_opts, 
			char *usage,
			unsigned int *ret_optindex, 
			char **ret_optname, 
			char **ret_optarg) {

  static int optindex   = 1;        /* init to 1 on first call  */
  static char *optptr = NULL;     /* ptr to next valid switch */

  unsigned int i, arglen, matches, opti = 0;

  /* Check to see if we've run out of options. '-' on its own is not an option */

  if (optindex >= argc || argv[optindex][0] != '-' || strcmp(argv[optindex], "-") == 0) { 
      *ret_optindex  = optindex; 
      *ret_optarg  = NULL; 
      *ret_optname = NULL; 
      return 0; 
  }

  /* We have an option. All options are treated as 'full' optiions,
     including single letter ones, which means you can't join
     single letter options with -abc; use -a -b -c */

  if (optptr == NULL && (argv[optindex][0] == '-')) {
      if ((optptr = strchr(argv[optindex], '=')) != NULL) { 
	*optptr = '\0'; 
	optptr++; 
      }

      arglen = strlen(argv[optindex]);
      matches = 0;
      for (i = 0; i < num_opts; i++) {
	if (strncmp(opt[i].name, argv[optindex], arglen) == 0) { 
	    matches++;
	    opti = i;
	}
      }
      if (matches > 1 && arglen != strlen(opt[opti].name)) 
	fatal_util("Option \'%s\' is ambiguous.\n%s", argv[optindex], usage);
      if (matches == 0)
	fatal_util("No such option as \'%s\'.\n%s", argv[optindex], usage);

      *ret_optname = opt[opti].name;

      if (opt[opti].type != NO_ARGS) {
	if (optptr != NULL) {
	  *ret_optarg = optptr;
	  optptr = NULL;
	  optindex++;
	}
	else if (optindex+1 >= argc) {
	  fatal_util("Option %s needs an argument\n%s", opt[opti].name, usage);
	}
	else {
	  *ret_optarg = argv[optindex+1];
	  optindex += 2;
	}
      }
      else {
	if (optptr != NULL) {
	  fatal_util("Option %s does not take any arguments\n%s", opt[optindex].name, usage);
	}
	*ret_optarg = NULL;
	optindex++;
      }
  }

  *ret_optindex = optindex;
  return 1;
}
