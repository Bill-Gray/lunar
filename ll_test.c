#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include "watdefs.h"
#include "mpc_func.h"

/* Run as,  e.g.,  './ll_test "n23 45 12.3, w31.4159, 200ft"' to test out
the lat/lon parsing code in 'mpc_code.cpp.'  Compile with

gcc -Wall -Wextra -pedantic -o ll_test ll_test.c liblunar.a -lm      */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

static void run_tests( void)
{
   FILE *ifile = fopen( "ll_test.txt", "rb");
   char buff[200];

   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      if( *buff != '#')
         {
         double lat, lon, alt;
         int expected_rval, rval;
         mpc_code_t cinfo;

         if( !memcmp( buff + 32, "COM Long.", 9))
            rval = get_xxx_location_info( &cinfo, buff + 32);
         else
            rval = get_lat_lon_info( &cinfo, buff + 32);
         assert( 4 == sscanf( buff, "%d %lf %lf %lf", &expected_rval, &lon, &lat, &alt));
         if( rval != expected_rval)
            fprintf( stderr, "Error with string '%s'; got %d\n", buff + 32, rval);
         else if( !rval)         /* should have succeeded */
            {
            const double lat1 = cinfo.lat * 180. / PI;
            const double lon1 = cinfo.lon * 180. / PI;
            const double rounding_tolerance = 1e-8;

            if( alt - cinfo.alt < -rounding_tolerance
                        || alt - cinfo.alt > rounding_tolerance
                        || lat - lat1 < -rounding_tolerance
                        || lat - lat1 > rounding_tolerance
                        || lon - lon1 < -rounding_tolerance
                        || lon - lon1 > rounding_tolerance)
               fprintf( stderr, "Error with string '%s'; got %.8lf %.8lf %.8lf\n",
                       buff + 32,
                       cinfo.lat * 180. / PI, cinfo.lon * 180. / PI, cinfo.alt);
            }
         }
   fclose( ifile);
}

int main( const int argc, const char **argv)
{
   if( argc > 1)
      {
      mpc_code_t cinfo;
      int rval;

      if( !memcmp( argv[1], "COM Long.", 9))
         rval = get_xxx_location_info( &cinfo, argv[1]);
      else
         rval = get_lat_lon_info( &cinfo, argv[1]);

      if( rval)
         printf( "Failed %d\n", rval);
      else
         {
         printf( "Lat %lf\nLon %lf\n", cinfo.lat * 180. / PI, cinfo.lon * 180. / PI);
         printf( "Alt %lf metres\n", cinfo.alt);
         }
      }
   else
      run_tests( );
   return( 0);
}
