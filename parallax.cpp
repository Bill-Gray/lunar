#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "watdefs.h"
#include "mpc_func.h"

/* Test program I run when I encounter a latitude/altitude and want to turn
it into parallax constants,  or vice versa.  See following 'error_exit'
function,  or run without command line arguments,  for usage. */

#define EARTH_MAJOR_AXIS_IN_METERS    6378137.
#define EARTH_MINOR_AXIS_IN_METERS    6356752.
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

static void error_exit( void)
{
   printf( "Run 'parallax' with three arguments (latitude, longitude, altitude)\n"
           "and they will be converted to parallax constants.  Run with two arguments\n"
           "(rho_cos_phi, rho_sin_phi) and the corresponding latitude and altitude\n"
           "will be computed and shown.\n");
   exit( -1);
}

static double get_angle( const char *ibuff)
{
   int place = 0;
   bool is_negative = false;
   double rval;

   if( *ibuff == '-')
      {
      ibuff++;
      is_negative = true;
      }
   while( ibuff[place] && ibuff[place] != '.')
      place++;
   if( !ibuff[place])         /* gotta have a decimal point */
      error_exit( );
   if( place < 6)
      rval = atof( ibuff);
   else
      {
      const long seconds = atol( ibuff);

      rval = (seconds / 10000L) * 3600L + ((seconds / 100L) % 100L) * 60L
                     + seconds % 100L;
      rval += atof( ibuff + place);
      rval /= 3600.;       /* cvt arcsec to degrees */
      }
   if( is_negative)
      rval = -rval;
   return( rval);
}

static char *show_angle( char *buff, double angle)
{
   const long microarcsec = (long)( fabs( angle) * 3600e+3);

   *buff = (angle < 0 ? '-' : '+');
   sprintf( buff + 1, "%02ld %02ld %02ld.%03ld",
               (microarcsec / 3600000L),         /* degrees */
               (microarcsec / 60000L) % 60L,     /* arcminutes */
               (microarcsec / 1000L) % 60L,     /* arcseconds */
                microarcsec % 1000L);           /* milliarcseconds */
   return( buff);
}

int main( const int argc, const char **argv)
{
   double lat, lon, alt, rho_cos_phi, rho_sin_phi;
   char buff[80];

   if( argc < 3 || argc > 4)
      error_exit( );
   else if( argc == 3)          /* parallax constants provided */
      {
      rho_cos_phi = atof( argv[1]);
      rho_sin_phi = atof( argv[2]);
      lat = point_to_ellipse( 1., EARTH_MINOR_AXIS_IN_METERS / EARTH_MAJOR_AXIS_IN_METERS,
                 rho_cos_phi, rho_sin_phi, &alt);
      alt *= EARTH_MAJOR_AXIS_IN_METERS;
      lon = 0.;
      }
   else                         /* lat/lon/alt provided */
      {
      lat = get_angle( argv[1]) * PI / 180.;
      lon = get_angle( argv[2]) * PI / 180.;
      alt = atof( argv[3]);
      lat_alt_to_parallax( lat, alt, &rho_cos_phi, &rho_sin_phi,
               EARTH_MAJOR_AXIS_IN_METERS, EARTH_MINOR_AXIS_IN_METERS);
      }
   lat *= 180. / PI;
   lon *= 180. / PI;
   if( lon)
      printf( "Longitude %11.7f = %s\n", lon, show_angle( buff, lon));
   printf( "Latitude  %11.7f = %s\n", lat, show_angle( buff, lat));
   printf( "Altitude %.3f meters\n", alt);
   printf( "Parallax constants %.7f %+.7f\n", rho_cos_phi, rho_sin_phi);
   return( 0);
}
