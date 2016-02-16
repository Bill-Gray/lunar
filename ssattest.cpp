#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "watdefs.h"
#include "lunar.h"

int main( const int argc, const char **argv)
{
   int i;
   double loc[3], jd;

   if( argc != 2)
      {
      printf( "'ssattest' takes a JD on the command line,  and outputs\n");
      printf( "coordinates for the eight main satellites of Saturn.\n");
      return( -1);
      }
   jd = atof( argv[1]);
   printf( "Date: %.5f\n", jd);
   for( i = 0; i < 8; i++)
      {
      jd = atof( argv[1]);
      calc_ssat_loc( jd, loc, i, 0L);
      printf( "%d: %9.6f %9.6f %9.6f\n", i, loc[0] * 100., loc[1] * 100.,
                  loc[2] * 100.);
      }
   return( 0);
}
