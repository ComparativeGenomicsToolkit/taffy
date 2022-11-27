/*  Last edited: Feb  1 15:25 2002 (klh) */
/**********************************************************************
 ** FILE: options.h
 ** DESCRIPTION:
 **  Rudimentary provision of command-line options
 **********************************************************************/

#ifndef _GETOPTIONS
#define _GETOPTIONS 


#include <string.h>
#include <stdlib.h>

#include "util.h"


/************************** constants *********************************/

#define NO_ARGS    0 
#define INT_ARG    1
#define FLOAT_ARG  2
#define CHAR_ARG   3
#define STRING_ARG 4

/******************* structure definitions ****************************/

struct Option {
  char *name;         /* name of option, e.g. "-option" */
  unsigned int type;  /* for typechecking, e.g. INT_ARG     */
};


/******************* function prototypes ****************************/

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
unsigned int get_option(int, char **, struct Option *, unsigned int,
			char *, unsigned int *, char **, char **);


#endif

