#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "watdefs.h"
#include "mpc_func.h"

/* Run as,  e.g.,  './ll_test "n23 45 12.3, w31.4159, 200ft"' to test out
the lat/lon parsing code in 'mpc_code.cpp.'  Compile with

gcc -Wall -Wextra -pedantic -o ll_test ll_test.c liblunar.a -lm      */


#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int main( const int argc, const char **argv)
{
   if( argc > 1)
      {
      mpc_code_t cinfo;
      const int rval = get_lat_lon_info( &cinfo, argv[1]);

      if( rval)
         printf( "Failed\n");
      else
         {
         printf( "Lat %lf\nLon %lf\n", cinfo.lat * 180. / PI, cinfo.lon * 180. / PI);
         printf( "Alt %lf metres\n", cinfo.alt);
         }
      }
   else
      printf( "Code to test lat/lon parsing\n");
   return( 0);
}
