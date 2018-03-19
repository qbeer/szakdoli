/*-------------------------------------------------------------*
     time.c  -   function to print out time
     Peter A. Henning, June 1994  otime_(int *retval) or otime(int *retval)
 *-------------------------------------------------------------*/
#include <time.h>
#include <stdlib.h>


otime(int *retval)
{  clock_t now;
   now = clock();
   *retval = (int)(100.0*now/CLOCKS_PER_SEC);
/*   printf(" CPU time used %g seconds\n",(1.0*now)/CLOCKS_PER_SEC) */;
}
