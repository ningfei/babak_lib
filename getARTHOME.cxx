#include <stdlib.h>
#include <stdio.h>

#define _getARTHOME

char *ARTHOME;

///////////////////////////////////////////////////////////////////////////////////////////////
// get the value of the ARTHOME environment variable
// The getenv() function searches the environment list for a string that matches "ARTHOME".
// It returns a pointer to the value in the environment, or NULL if there is no match.
///////////////////////////////////////////////////////////////////////////////////////////////
void getARTHOME()
{
   ARTHOME=getenv("ARTHOME");

   if(ARTHOME == NULL)
   {
      printf("The ARTHOME environment variable is not defined.\n");
      exit(0);
   }
}
