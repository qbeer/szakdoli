/*
Write out buffer to stdout and stderr immediately.
*/

#include <stdio.h>

f77flush_()
{
  fflush(stdout);
  fflush(stderr);
}
