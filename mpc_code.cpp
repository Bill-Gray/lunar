#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "watdefs.h"
#include "mpc_func.h"
#include "lunar.h"

#define SUN_RADIUS          696000e+3
#define MERCURY_RADIUS      2439700.
#define VENUS_RADIUS        6051800.
#define EARTH_MAJOR_AXIS    6378137.
#define EARTH_MINOR_AXIS    6356752.
#define MARS_MAJOR_AXIS     3396190.
#define MARS_MINOR_AXIS     3376200.
#define MOON_RADIUS         1737400.
#define JUPITER_MAJOR_AXIS   71492e+3
#define JUPITER_MINOR_AXIS   66854e+3
#define SATURN_MAJOR_AXIS    60268e+3
#define SATURN_MINOR_AXIS    54364e+3
#define URANUS_MAJOR_AXIS    25559e+3
#define URANUS_MINOR_AXIS    24973e+3
#define NEPTUNE_MAJOR_AXIS   24764e+3
#define NEPTUNE_MINOR_AXIS   24341e+3
#define PLUTO_RADIUS          1195e+3
#define IO_MEAN_RADIUS        1821.49e+3
#define EUROPA_MEAN_RADIUS    1560.8e+3
#define GANYMEDE_MEAN_RADIUS  2631.2e+3
#define CALLISTO_MEAN_RADIUS  2410.3e+3

      /* Earth dimensions are WGS84 constants */
      /* Other sizes are from http://adsabs.harvard.edu/abs/2011CeMDA.109..101A */
      /* or http://astropedia.astrogeology.usgs.gov/alfresco/d/d/workspace/SpacesStore/28fd9e81-1964-44d6-a58b-fbbf61e64e15/WGCCRE2009reprint.pdf */

static const double equatorial_radii[15] = {
      SUN_RADIUS, MERCURY_RADIUS, VENUS_RADIUS, EARTH_MAJOR_AXIS,
      MARS_MAJOR_AXIS, JUPITER_MAJOR_AXIS, SATURN_MAJOR_AXIS,
      URANUS_MAJOR_AXIS, NEPTUNE_MAJOR_AXIS, PLUTO_RADIUS,
      MOON_RADIUS, IO_MEAN_RADIUS, EUROPA_MEAN_RADIUS,
      GANYMEDE_MEAN_RADIUS, CALLISTO_MEAN_RADIUS };

double planet_radius_in_meters( const int planet_idx)
{
   assert( planet_idx >= 0 && planet_idx < 15);
   return( equatorial_radii[planet_idx]);
}

/* NOTE that for the Earth,  I'm using the WGS84 ellipsoid.  The
IAU1976 ellipsoid has a major axis of 6378.1377
  with
major axis 6378.140 km (see above) and flattening 1/f = 298.257.
The GRS80 and WGS84 ellipsoids have a major axis of 6378.137 meters
(three meters less than the IAU1976 value) and 1/f = 298.257223563
and 298.257222101 respectively (making for semimajor axes also about
three meters less than the IAU1976 value).      */

#define N_POLAR_RADII      9

const double polar_radii[N_POLAR_RADII] = {
      SUN_RADIUS, MERCURY_RADIUS, VENUS_RADIUS, EARTH_MINOR_AXIS,
      MARS_MINOR_AXIS, JUPITER_MINOR_AXIS, SATURN_MINOR_AXIS,
      URANUS_MINOR_AXIS, NEPTUNE_MINOR_AXIS };

double planet_axis_ratio( const int planet_idx)
{
   return( planet_idx >= N_POLAR_RADII ?
            1 : polar_radii[planet_idx] / equatorial_radii[planet_idx]);
}

const double PI =
   3.1415926535897932384626433832795028841971693993751058209749445923;

int lat_alt_to_parallax( const double lat, const double ht_in_meters,
            double *rho_cos_phi, double *rho_sin_phi,
            const double major_axis_in_meters,
            const double minor_axis_in_meters)
{
   const double axis_ratio = minor_axis_in_meters / major_axis_in_meters;
   const double u = atan2( sin( lat) * axis_ratio, cos( lat));

   *rho_sin_phi = axis_ratio * sin( u) +
                            (ht_in_meters / major_axis_in_meters) * sin( lat);
   *rho_cos_phi = cos( u) + (ht_in_meters / major_axis_in_meters) * cos( lat);
   return( 0);
}

/* Given an ellipse with semimajor axis a,  semiminor axis b,  centered
at the origin,  and an arbitrary point (x, y),  point_to_ellipse() will
compute the closest distance between that point and the ellipse,  and
the angle to the ellipse.

   This is an exact method from _Explanatory Supplement to the Astronomical
Almanac_, pgs 206-207,  in turn from K. M. Borkowski (1989), "Accurate
Algorithms to Transform Geocentric to Geodetic Coordinates",  _Bulletin
Geodesique_ 63, no. 1, 50-56,  modified slightly to handle the possibilities
of negative x and/or y.  It is also described at

http://www.astro.uni.torun.pl/~kb/Papers/ASS/Geod-ASS.htm

This reduces the problem to finding the roots of a quartic polynomial,
but does so in a form that is somewhat straightforward,  with unit
leading and trailing coefficients and a zero quadratic coefficient.

References are to the _Explanatory Supplement_ and then the above URL.
For example,  the equation for 'e' is given at 4.22-12 in the ES and
as equation (6) at the above URL.    */

double point_to_ellipse( const double a, const double b,
                         const double x, const double y, double *dist)
{
   const double fy = fabs( y), fx = fabs( x);
   double lat;

   if( x == 0.)
      {
      lat = PI / 2.;
      if( dist)
         *dist = fy - a;
      }
   else
      {
      const double c_squared = a * a - b * b;
      const double e = (b * fy - c_squared) / (a * fx);      /* 4.22-12/6 */
      const double f = (b * fy + c_squared) / (a * fx);      /* 4.22-13/7 */
      const double p = (4. / 3.) * (e * f + 1.);             /* 4.22-14/9 */
      const double q = 2. * (e * e - f * f);                 /* 4.22-15/10 */
      const double d = p * p * p + q * q;                    /* 4.22-16/12 */
      double v, g, t;

      if( d >= 0.)
         {
         const double sqrt_d = sqrt( d);

         v = cbrt( sqrt_d - q) - cbrt( sqrt_d + q);         /* 4.22-17/11a */
         }
      else
         {
         const double sqp = sqrt( -p);
         const double temp_ang = acos( q / (sqp * p));

         v = 2. * sqp * cos( temp_ang / 3.);                      /* 11b */
         }
      g = (sqrt( e * e + v) + e) * .5;                     /* 4.22-18/14 */
      t = sqrt( g * g + (f - v * g) / (2. * g - e)) - g;   /* 4.22-19/13 */
      lat = atan2( a * (1. - t * t), 2. * b * t);          /* 4.22-20/15a */
      if( dist)                                            /* 4.22-21/15b */
         *dist = (fx - a * t) * cos( lat) + (fy - b) * sin( lat);
      }
   if( x < 0.)
      lat = PI - lat;
   if( y < 0.)
      lat = -lat;
   return( lat);
}

/* You can store locations in 'rovers.txt' in base-60 form,  with the
lat and longitude smashed together;  e.g.,  19 13' 33.1" would be
stored as 191333.1.  The following code would take 191331.1 as input
and return 19 + 13/60 + 33.1/3600 = 19.22586111.   */

static double convert_base_60_to_decimal( const double ival)
{
   const int deg = (int)( ival / 10000.);
   const int min = (int)( fmod( ival / 100., 100.));
   const double sec = fmod( ival, 100.);
   const double rval = (double)deg + (double)min / 60. + sec / 3600.;

   return( rval);
}

int get_mpc_code_info( mpc_code_t *cinfo, const char *buff)
{
   int i = 0, rval = -1;

   while( buff[i] > ' ' && buff[i] <= '~')
      i++;

   cinfo->lat = cinfo->lon = cinfo->alt
                  = cinfo->rho_sin_phi = cinfo->rho_cos_phi = 0.;
   if( i == 3 && strlen( buff) > 33 && buff[3] == ' ')
      {
      rval = 3;         /* assume earth */

      while( buff[i] == ' ')
         i++;
      if( buff[4] == '!')     /* rovers.txt format */
         {
         if( sscanf( buff + 5, "%lf%lf%lf", &cinfo->lon,
                               &cinfo->lat, &cinfo->alt) != 3)
            rval = -1;
         else
            {
            const char *tptr = strchr( buff, '@');

            if( fabs( cinfo->lon) > 361. || fabs( cinfo->lat) > 91.)
               {
               cinfo->lon = convert_base_60_to_decimal( cinfo->lon);
               cinfo->lat = convert_base_60_to_decimal( cinfo->lat);
               }
            cinfo->lon *= PI / 180.;
            cinfo->lat *= PI / 180.;
            cinfo->name = buff + 47;
            cinfo->format = MPC_CODE_LAT_LON_ALT;
            if( tptr)                     /* non-earth location */
               rval = atoi( tptr + 1);
            if( cinfo->lat && cinfo->lon)   /* i.e.,  topocentric */
               {
               const double major = planet_radius_in_meters( rval);
               const double minor = major * planet_axis_ratio( rval);

               lat_alt_to_parallax( cinfo->lat, cinfo->alt,
                        &cinfo->rho_cos_phi, &cinfo->rho_sin_phi,
                        major, minor);
               }
            }
         }
      else if( buff[7] == '.' && (buff[21] == '+' || buff[21] == '-')
               && buff[14] == '.' && buff[23] == '.')
         {                 /* 'standard' MPC format */
         if( sscanf( buff + 3, "%lf%lf%lf", &cinfo->lon,
                  &cinfo->rho_cos_phi, &cinfo->rho_sin_phi) != 3)
            rval = -1;
         else
            {
            cinfo->name = buff + 30;
            cinfo->format = MPC_CODE_PARALLAXES;
            cinfo->lon *= PI / 180.;
            cinfo->lat = point_to_ellipse( 1., EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS,
                    cinfo->rho_cos_phi, cinfo->rho_sin_phi, &cinfo->alt);
            cinfo->alt *= EARTH_MAJOR_AXIS;
            }
         }
      else if( i == 30)
         {
         cinfo->name = buff + 30;
         cinfo->format = MPC_CODE_SATELLITE;
         }
      else
         rval = -1;
      }
   if( rval != -1)
      {
      cinfo->planet = rval;
      memcpy( cinfo->code, buff, 3);
      cinfo->code[3] = '\0';
      if( cinfo->lon < 0.)
         cinfo->lon += PI + PI;
      }
   return( rval);
}

#ifdef TEST_CODE

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( (argc == 1 ? "ObsCodes.htm" : argv[1]), "rb");
   char buff[200];
   mpc_code_t code;

   const char *header =
            "Cod Longitude  Latitude     Altitude rho_cos   rho_sin_phi";

   if( !ifile)
      ifile = fopen( "ObsCodes.html", "rb");
   if( !ifile)
      {
      printf( "ObsCodes not opened\n");
      return( -1);
      }
   printf( "%s\n", header);
   while( fgets( buff, sizeof( buff), ifile))
      if( get_mpc_code_info( &code, buff) >= 0)
         printf( "%s %10.6f %+10.6f %10.3f %9.7f %+10.7f %s",
                  code.code, code.lon * 180. / PI, code.lat * 180. / PI,
                  code.alt, code.rho_cos_phi, code.rho_sin_phi, code.name);
   printf( "%s\n", header);
   fclose( ifile);
   return( 0);
}
#endif
