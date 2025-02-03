/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include "watdefs.h"
#include "mpc_func.h"
#include "lunar.h"
#include "stringex.h"

#ifdef _MSC_VER
     /* Microsoft Visual C/C++ has no strncasecmp.  strncmp will do.  */
#define strncasecmp strncmp
#endif

#define SUN_RADIUS          695700e+3
#define MERCURY_MAJOR_AXIS  2440530.
#define MERCURY_MINOR_AXIS  2438260.
#define VENUS_RADIUS        6051800.
#define EARTH_MAJOR_AXIS    6378137.
#define EARTH_MINOR_AXIS    6356752.314140347
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
#define PLUTO_RADIUS          1188300.
#define IO_MEAN_RADIUS        1821.49e+3
#define EUROPA_MEAN_RADIUS    1560.8e+3
#define GANYMEDE_MEAN_RADIUS  2631.2e+3
#define CALLISTO_MEAN_RADIUS  2410.3e+3

      /* Earth dimensions are WGS84 constants */
      /* Other sizes are from http://adsabs.harvard.edu/abs/2011CeMDA.109..101A */
      /* or https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf  */

#define N_EQUATORIAL_RADII 15

static const double equatorial_radii[N_EQUATORIAL_RADII] = {
      SUN_RADIUS, MERCURY_MAJOR_AXIS, VENUS_RADIUS, EARTH_MAJOR_AXIS,
      MARS_MAJOR_AXIS, JUPITER_MAJOR_AXIS, SATURN_MAJOR_AXIS,
      URANUS_MAJOR_AXIS, NEPTUNE_MAJOR_AXIS, PLUTO_RADIUS,
      MOON_RADIUS, IO_MEAN_RADIUS, EUROPA_MEAN_RADIUS,
      GANYMEDE_MEAN_RADIUS, CALLISTO_MEAN_RADIUS };

double planet_radius_in_meters( const int planet_idx)
{
   if( planet_idx >= 0 && planet_idx < N_EQUATORIAL_RADII)
      return( equatorial_radii[planet_idx]);
   else
      return( 0.);
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
      SUN_RADIUS, MERCURY_MINOR_AXIS, VENUS_RADIUS, EARTH_MINOR_AXIS,
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

/* MS only got around to adding cbrt in VS2013 : */

#if (defined( _MSC_VER) && (_MSC_VER < 1800)) || defined( __WATCOMC__)

static double cbrt( const double z)
{
   double rval;

   if( z > 0.)
      rval = pow( z, 1. / 3.);
   else if( z < 0.)
      rval = -pow( -z, 1. / 3.);
   else
      rval = 0.;
   return( rval);
}
#endif

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
as equation (6) at the above URL.

The same point-to-ellipse problem can come up in computing MOIDs, which
is why the function is not of type static : it's used in moid.cpp. */

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

/* atof( ) and strtod( ) have to handle a lot of odd cases,  such as
exponentials.  The requirements for rounding are very strict.  As a
result,  it is quite slow.  If you are sure your input string has no
odd cases (and does not need real error checking), this is faster and
"good enough". It only handles strings with whitespace,  digits,  and
optional decimal point and more digits.  See 'f_strtod.cpp' in the
Bill-Gray/jpl_eph repository for further discussion,  and a more
robust solution.        */

double quick_strtod( const char *ibuff, const char **endptr)
{
   bool is_negative = false, got_a_digit = false;
   int integer_part = 0;
   double rval;

   if( endptr)
      *endptr = ibuff;
   while( *ibuff == ' ')
      ibuff++;
   if( *ibuff == '-')
      {
      is_negative = true;
      ibuff++;
      }
   else if( *ibuff == '+')
      ibuff++;
   if( *ibuff >= '0' && *ibuff <= '9')
      got_a_digit = true;
   while( *ibuff >= '0' && *ibuff <= '9')
      integer_part = integer_part * 10 + (*ibuff++ - '0');
   if( *ibuff != '.')
      rval = (double)integer_part;
   else
      {
      int frac_part = 0, divisor = 1;

      ibuff++;
      if( *ibuff < '0' || *ibuff > '9')
         if( !got_a_digit)       /* it's not a number */
            return( 0.);
      while( *ibuff >= '0' && *ibuff <= '9')
         {
         frac_part = frac_part * 10 + (*ibuff++ - '0');
         divisor *= 10;
         }
      rval = (double)integer_part + (double)frac_part / (double)divisor;
      }
   if( endptr)
      *endptr = ibuff;
   if( is_negative)
      rval = -rval;
   return( rval);
}

double quick_atof( const char *ibuff)
{
   return( quick_strtod( ibuff, NULL));
}

/* Useful if floating-point text may be immediately followed by more
floating-point text,  and you want scanning to stop after 'len'
digits.  This happens in ObsCodes.htm.    */

static double limited_atof( const char *ibuff, const size_t len)
{
   char tbuff[80];

   assert( len < sizeof( tbuff) - 1);
   memcpy( tbuff, ibuff, len);
   tbuff[len] = '\0';
   return( quick_strtod( tbuff, NULL));
}

/* You can store locations in 'rovers.txt' in base-60 form,  with the
degrees/minutes/seconds smashed together;  e.g.,  19 13' 33.1" would be
stored as 191333.1.  The following code would take 191331.1 as input
and return 19 + 13/60 + 33.1/3600 = 19.22586111.   */

static double convert_base_60_to_decimal( const double ival)
{
   const int secs = (int)ival;
   const double rval = (double)( secs / 10000)
                     + (double)((secs / 100) % 100) / 60.
                     + (double)( secs % 100) / 3600.
                     + (ival - (double)secs) / 3600.;

   return( rval);
}

static void _set_parallax_constants( mpc_code_t *cinfo)
{
   const double major = planet_radius_in_meters( cinfo->planet);
   const double minor = major * planet_axis_ratio( cinfo->planet);

   lat_alt_to_parallax( cinfo->lat, cinfo->alt,
                        &cinfo->rho_cos_phi, &cinfo->rho_sin_phi,
                        major, minor);
}

int get_mpc_code_info( mpc_code_t *cinfo, const char *buff)
{
   int i = 0, rval = -1;
   const size_t slen = strlen( buff);

   while( buff[i] > ' ' && buff[i] <= '~' && buff[i] != '!')
      i++;
   memset( cinfo, 0, sizeof( mpc_code_t));
   memcpy( cinfo->code, buff, 4);
   if( i >= 3 && i <= 4 && slen >= 30)
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
            const char *tptr = strchr( buff + 4, '@');

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
               {
               rval = atoi( tptr + 1);
               if( rval > 9000)        /* dynamical point */
                  cinfo->lat = cinfo->lon = cinfo->alt = 0.;
               }
            cinfo->planet = rval;
            if( cinfo->lat || cinfo->lon || cinfo->alt)   /* i.e.,  topocentric */
               _set_parallax_constants( cinfo);
            }
         }
      else if( buff[7] == '.' && strchr( "+- ", buff[21])
               && buff[14] == '.' && buff[23] == '.' && buff[3] == ' ')
         {                 /* 'standard' MPC format */
         cinfo->lon = limited_atof( buff + 4, 9);
         cinfo->rho_cos_phi = limited_atof( buff + 13, 8);
         cinfo->rho_sin_phi = limited_atof( buff + 21, 9);
         cinfo->name = buff + 30;
         cinfo->format = MPC_CODE_PARALLAXES;
         cinfo->lon *= PI / 180.;
         if( cinfo->rho_cos_phi || cinfo->rho_sin_phi)
            cinfo->lat = point_to_ellipse( 1., EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS,
                 cinfo->rho_cos_phi, cinfo->rho_sin_phi, &cinfo->alt);
         else
            cinfo->lat = cinfo->alt = 0.;
         cinfo->alt *= EARTH_MAJOR_AXIS;
         while( cinfo->prec1 < 5 && isdigit( buff[8 + cinfo->prec1]))
            cinfo->prec1++;      /* longitude precision,  in digits */
         while( cinfo->prec2 < 6 && isdigit( buff[15 + cinfo->prec2]))
            cinfo->prec2++;      /* parallax precision,  in digits */
         }
      else if( i == 30)
         {
         cinfo->name = buff + 30;
         cinfo->format = MPC_CODE_SATELLITE;
         rval = -2;
         }
      else
         rval = -1;
      }
   else if( slen >= 80 && buff[14] == 'v' && buff[80] < ' '
               && buff[48] == '.' && buff[37] == '.'
               && buff[77] == '2' && (buff[45] == '+' || buff[45] == '-'))
      {              /* roving observer line */
      double lat, lon, alt;

      if( 3 == sscanf( buff + 34, "%lf %lf %lf%n", &lon, &lat, &alt, &i))
         {
         cinfo->lon = lon * PI / 180.;
         cinfo->lat = lat * PI / 180.;
         cinfo->alt = alt;
         cinfo->name = "Roving observer";
         cinfo->planet = rval = 3;
         _set_parallax_constants( cinfo);
         memcpy( cinfo->code, buff + 77, 3);
         }
      }
   if( cinfo->code[3] == ' ')        /* standard 3-character code */
      cinfo->code[3] = '\0';
   else                       /* 'extended' 4-character code */
      cinfo->code[4] = '\0';
   if( rval != -1)
      {
      cinfo->planet = rval;
      if( cinfo->lon < 0.)
         cinfo->lon += PI + PI;
      }
   else        /* fail to 'reasonable' geocentric default */
      {
      cinfo->name = "Unkn";
      cinfo->planet = 3;
      }
   return( rval);
}

int extract_region_data_for_lat_lon( FILE *ifile, char *buff,
            const double lat_in_degrees, const double lon_in_degrees)
{
   char tbuff[90];
   size_t i = 0;

   *buff = '\0';
   while( !*buff && fgets( tbuff, sizeof( tbuff), ifile))
      if( *tbuff != '#')
         {
         double d_lon1 = quick_atof( tbuff)      - lon_in_degrees;
         double d_lon2 = quick_atof( tbuff + 20) - lon_in_degrees;

         while( d_lon1 > 180.)
            d_lon1 -= 360.;
         while( d_lon1 < -180.)
            d_lon1 += 360.;
         while( d_lon2 - d_lon1 > 180.)
            d_lon2 -= 360.;
         while( d_lon2 - d_lon1 < -180.)
            d_lon2 += 360.;
         if( d_lon1 * d_lon2 < 0.)
            {
            const double d_lat1 = quick_atof( tbuff + 10) - lat_in_degrees;
            const double d_lat2 = quick_atof( tbuff + 30) - lat_in_degrees;

            if( d_lat1 * d_lat2 < 0.)
               {
               strcpy( buff, tbuff + 40);
               while( buff[i] >= ' ')
                  i++;
               buff[i] = '\0';   /* remove trailing CR/LF */
               }
            }
         }
   return( *buff ? 0 : -2);
}

static double _get_number_for_angle( const char *buff, int *nbytes, int *is_integer)
{
   double rval;
   const char *tptr = buff, *endptr;

   while( *tptr == ' ')
      tptr++;
   if( !isdigit( *tptr))
      {
      *nbytes = 0;
      *is_integer = -1;
      return( 0.);
      }
   while( isdigit( *tptr))
      tptr++;
   rval = quick_strtod( buff, &endptr);
   *nbytes = (int)( endptr - buff);
   if( !*nbytes || rval < 0.)
      {
      *is_integer = -1;
      return( 0.);
      }
   *is_integer = (*tptr != '.');
   return( rval);
}

static double _get_angle( const char *buff, int *nbytes, int *n_fields)
{
   double rval;
   const char *tptr = buff;
   int is_integer, dummy_nbytes;

   assert( n_fields);
   if( !nbytes)
      nbytes = &dummy_nbytes;
   *n_fields = 0;
   rval = _get_number_for_angle( buff, nbytes, &is_integer);
   if( is_integer != -1)
      {
      tptr += *nbytes;
      *n_fields = 1;
      if( is_integer)
         {
         double fraction = _get_number_for_angle( tptr, nbytes, &is_integer);

         if( fraction > 60.)
            {
            *n_fields = 0;
            return( 0.);
            }
         if( is_integer != -1)
            {
            tptr += *nbytes;
            rval += fraction / 60.;
            *n_fields = 2;
            if( is_integer)
               {
               fraction = _get_number_for_angle( tptr, nbytes, &is_integer);
               if( fraction > 60.)
                  {
                  *n_fields = 0;
                  return( 0.);
                  }
               if( is_integer != -1)
                  {
                  tptr += *nbytes;
                  rval += fraction / 3600.;
                  *n_fields = 3;
                  }
               }
            }
         }
      }
   if( tptr != buff)
      while( *tptr == ' ')    /* skip trailing spaces */
         tptr++;
   *nbytes = (int)( tptr - buff);
   return( rval);
}

/* Return value is GOT_LAT if we got a latitude, GOT_LON if a longitude,
GOT_ALT if it was an altitude,  and GOT_NOTHING otherwise.     */

#define GOT_LON            0
#define GOT_LAT            1
#define GOT_ALT            2
#define GOT_NOTHING       -1

static int extract_lat_lon( const char *buff, size_t *bytes_read, double *value)
{
   const char *tptr = buff;
   int rval = GOT_NOTHING, nbytes, n_fields;
   bool is_negative = false;
   const char *compass = "nNeEsSwW";
   char compass_byte = 0;

   while( *tptr == ' ')
      tptr++;
   if( !*tptr)
      return( GOT_NOTHING);
   if( strchr( compass, *tptr))
      compass_byte = *tptr++;
   if( !strncasecmp( tptr, "alt", 3))
      {
      tptr += 3;
      if( *tptr == '.')
         tptr++;
      }
   *value = _get_angle( tptr, &nbytes, &n_fields);
   tptr += nbytes;
   if( !compass_byte && n_fields == 1)        /* possible height */
      {
      if( *tptr == ' ')
         tptr++;
      if( *tptr == 'm')
         {
         tptr++;
         rval = GOT_ALT;
         }
      else if( *tptr == 'f' && tptr[1] == 't')
         {
         const double intl_feet_to_meters = 0.3048;
/*       const double us_survey_feet_to_meters = 12. / 39.37;   */

         tptr += 2;               /* alt given in feet (shudder) */
         rval = GOT_ALT;
         *value *= intl_feet_to_meters;
         }
      }
/* printf( "From '%s',  got %d, value %f\n", buff, (int)( tptr - buff), *value); */
   if( n_fields && rval == GOT_NOTHING)
      {
      if( !compass_byte && strchr( compass, *tptr))
         {                 /* compass sign after the angle */
         compass_byte = *tptr++;
/*       printf( "Got byte '%c'; value %f\n", compass_byte, *value); */
         }
      compass_byte = tolower( compass_byte);
      if( compass_byte == 'n' || compass_byte == 's')
         {
         is_negative = (compass_byte == 's');
         rval = GOT_LAT;
         }
      else if( compass_byte == 'e' || compass_byte == 'w')
         {
         is_negative = (compass_byte == 'w');
         rval = GOT_LON;
         }
      }
   if( is_negative)
      *value = -*value;
   if( bytes_read)
      {
      if( *tptr == ',')
         tptr++;
      if( rval == GOT_NOTHING)
         *bytes_read = 0;
      else
         *bytes_read = tptr - buff;
      }
   return( rval);
}

/* Given text of the forms...

n44.01, W69.9
44 01 13.2 N 69 54 1.7W
E223.456,s56 20 23.3, Alt. 1700m
w 69.91, n 44.012, 1700ft

   etc.,  i.e.,  a lat/lon that would be readable by a human,
plus an optional altitude,  this function will puzzle it out.
Returns 0 on success,  -1 if it couldn't parse it.  As the above
shows,  the compass sign can be at the start or end of the text,
and one can use decimal degrees,  or degrees/decimal minutes,  or
degrees/minutes/decimal seconds.

   Feet are assumed to be 'international' feet,  1 foot = 30.48 cm,
not US survey feet.  (The difference is two parts per million.) */

int get_lat_lon_info( mpc_code_t *cinfo, const char *buff)
{
   double second_value, first_value, alt_value;
   size_t n_bytes;
   const int piece_type = extract_lat_lon( buff, &n_bytes, &first_value);
   int rval = -1;

   if( piece_type == GOT_LAT || piece_type == GOT_LON)
      {
      const char *tptr = buff + n_bytes;

      if( extract_lat_lon( tptr, &n_bytes, &second_value) ==
                                         GOT_LAT + GOT_LON - piece_type)
         {
         if( cinfo)
            {
            tptr += n_bytes;
            first_value *= PI / 180.;
            second_value *= PI / 180.;
            if( piece_type == GOT_LON)
               {
               cinfo->lon = first_value;
               cinfo->lat = second_value;
               }
            else
               {
               cinfo->lat = first_value;
               cinfo->lon = second_value;
               }
            if( cinfo->lon < 0.)
               cinfo->lon += PI + PI;
            if( cinfo->lat > PI / 2. || cinfo->lat < -PI / 2.
                     || cinfo->lon < 0. || cinfo->lon > PI + PI)
               return( -1);
            if( extract_lat_lon( tptr, &n_bytes, &alt_value) == GOT_ALT)
               cinfo->alt = alt_value;
            else              /* use a default 100m altitude if none specified */
               cinfo->alt = 100.;
            lat_alt_to_parallax( cinfo->lat, cinfo->alt,
                        &cinfo->rho_cos_phi, &cinfo->rho_sin_phi,
                        EARTH_MAJOR_AXIS, EARTH_MINOR_AXIS);
            cinfo->planet = 3;
            cinfo->name = buff;
            cinfo->format = MPC_CODE_LAT_LON_ALT;
            cinfo->prec1 = cinfo->prec2 = 0;
            strlcpy_error( cinfo->code, "Rov");
            }
         rval = 0;
         }
      }
   return( rval);
}

double get_ra_from_string( const char *buff, int *bytes_read)
{
   int n_fields;
   double rval = _get_angle( buff, bytes_read, &n_fields);

   if( n_fields > 2)     /* assume HH MM.mm or HH MM SS.sss */
      rval *= 15.;
   if( !n_fields || rval > 360.)
      *bytes_read = 0;
   return( rval * PI / 180.);
}

double get_dec_from_string( const char *buff, int *bytes_read)
{
   if( *buff != '+' && *buff != '-')
      {
      *bytes_read = 0;
      return( 0.);
      }
   else
      {
      int n_fields;
      double rval = _get_angle( buff + 1, bytes_read, &n_fields);

      if( !n_fields || rval > 90.)
         *bytes_read = 0;
      if( *buff == '-')
         rval = -rval;
      if( n_fields)
         (*bytes_read)++;
      return( rval * PI / 180.);
      }
}

/* https://www.minorplanetcenter.net/iau/info/Astrometry.html#HowObsCode
suggests that you start out using "observatory code" (XXX),  whilst
including a comment such as

COM Long. 239 18 45 E, Lat. 33 54 11 N, Alt. 100m, Google Earth

   Here,  the lat/lon are extracted using the above functions,  so the
formatting can be considerably more flexible (decimal degrees or
degrees/decimal minutes,  compass sign at start or end).  But don't
send such locations to MPC.               */

#define MALFORMED_XXX_LINE    (-2)

int get_xxx_location_info( mpc_code_t *cinfo, const char *buff)
{
   int rval = 0;

   if( strlen( buff) < 9 || memcmp( buff, "COM Long.", 9))
      rval = -1;
   else
      {
      const char *lat_ptr = strstr( buff, "Lat.");
      const char *alt_ptr = strstr( buff, "Alt.");

      if( !lat_ptr)
         rval = -4;
      if( !alt_ptr)
         rval -= 8;
      if( lat_ptr && alt_ptr
            && GOT_LON == extract_lat_lon( buff + 9, NULL, &cinfo->lon)
            && GOT_LAT == extract_lat_lon( lat_ptr + 4, NULL, &cinfo->lat))
         {
         if( cinfo->lon < 0.)
            cinfo->lon += 360.;
         cinfo->lon *= PI / 180.;
         cinfo->lat *= PI / 180.;
         cinfo->alt = atof( alt_ptr + 4);
         cinfo->format = MPC_CODE_LAT_LON_ALT;
         cinfo->prec1 = cinfo->prec2 = 0;
         cinfo->name = "Temporary MPC code";
         cinfo->planet = 3;
         strlcpy_error( cinfo->code, "XXX");
         _set_parallax_constants( cinfo);
         }
      else
         rval = MALFORMED_XXX_LINE;
      }
   return( rval);
}

#ifdef TEST_CODE

static int _text_search_and_replace( char *str, const char *oldstr,
                                     const char *newstr)
{
   size_t ilen = strlen( str), rval = 0;
   const size_t oldlen = strlen( oldstr);
   const size_t newlen = strlen( newstr);

   while( ilen >= oldlen)
      if( !memcmp( str, oldstr, oldlen))
         {
         memmove( str + newlen, str + oldlen, ilen - oldlen + 1);
         memcpy( str, newstr, newlen);
         str += newlen;
         ilen -= oldlen;
         rval++;
         }
      else
         {
         str++;
         ilen--;
         }
   return( (int)rval);
}

static FILE *fopen_ext( const char *exec_path, const char *filename, const char *permits)
{
   char obuff[255];
   size_t i;
   FILE *rval;

   strlcpy_error( obuff, exec_path);
   i = strlen( obuff);
   while( i && obuff[i - 1] != '/' && obuff[i - 1] != '\\')
      i--;
   obuff[i] = '\0';
   strlcat_err( obuff + i, filename, sizeof( obuff));
   rval = fopen( obuff, permits);
// if( !rval)
//    rval = fopen( obuff + i, permits);
   return( rval);
}

const char *kml_header_text =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
    "<Folder>\n"
    "      <name> MPC sites </name>\n"
    "      <description> MPC observatory sites </description>\n";

/* TODO : add in 'night earth' maps such as

https://www.nightearth.com/?@44.01,-69.9,11z&data=$bWVsMmQx    */

int main( const int argc, const char **argv)
{
   FILE *ifile;
   const char *ifilename = "ObsCodes.htm";
   char buff[200];
   mpc_code_t code;
   bool google_map_links = false, dump_comments = false;
   bool make_kml = false;
   int i, line_no = 0;
   size_t google_offset = 0;

   const char *header =
            "Pl Code Longitude  Latitude     Altitude rho_cos   rho_sin_phi  region";

   for( i = 1; i < argc; i++)
      if( argv[i][0] != '-')
         ifilename = argv[i];
      else
         switch( argv[i][1])
            {
            case 'v':
               dump_comments = true;
               break;
            case 'g':
               {
               const time_t t0 = time( NULL);
               FILE *hdr_ifile = fopen_ext( argv[0], "mpc_hdr.htm", "rb");
               bool showing_lines = true;

               assert( hdr_ifile);
               google_map_links = true;
               google_offset = 3;
               while( fgets( buff, sizeof( buff), hdr_ifile))
                  if( strstr( buff, "%.24s"))
                     printf( buff, asctime( gmtime( &t0)));
                  else if( *buff == '#')
                     {
                     if( strstr( ifilename, "rovers") &&
                                 !memcmp( buff, "# 'rovers' ", 11))
                        showing_lines = (buff[12] == 'n');
                     }
                  else if( showing_lines)
                     printf( "%s", buff);
               fclose( hdr_ifile);
               }
               break;
            case 'k':
               {
               make_kml = true;
               printf( "%s", kml_header_text);
               }
               break;
            default:
               printf( "Command line option '%s' unrecognized\n", argv[i]);
               return( -1);
            }
   ifile = fopen( ifilename, "rb");
   if( !ifile)
      ifile = fopen( "ObsCodes.html", "rb");
   if( !ifile)
      ifile = fopen_ext( argv[0], "ObsCodes.html", "rb");
   if( !ifile)
      {
      printf( "ObsCodes not opened\n");
      return( -1);
      }
   if( !make_kml)
      printf( "%s\n", header + google_offset);
   while( fgets( buff, sizeof( buff), ifile))
      if( get_mpc_code_info( &code, buff) != -1 || !get_xxx_location_info( &code, buff))
         {
         char region[100], obuff[200];
         bool show_link_for_this_line;

         code.lat *= 180. / PI;
         code.lon *= 180. / PI;
         *region = '\0';
         for( i = 0; buff[i] >= ' '; i++)
            ;
         buff[i] = '\0';
         if( code.planet == 3)
            {
            const char *geo_rect_filename = "geo_rect.txt";
            FILE *geo_rect_fp = fopen_ext( argv[0], geo_rect_filename, "rb");

            if( geo_rect_fp)
               {
               extract_region_data_for_lat_lon( geo_rect_fp, region,
                    code.lat, code.lon);
               fclose( geo_rect_fp);
               }
            }
         else                             /* not on Earth */
            snprintf_err( region, sizeof( region), "@%d", code.planet);
         snprintf_err( obuff, sizeof( obuff),
                 "%2d %-4s %10.6f %+10.6f %10.3f %9.7f %+10.7f %-15.15s ",
                  code.planet % 100,
                  code.code, code.lon, code.lat,
                  code.alt, code.rho_cos_phi, code.rho_sin_phi,
                  region);
         if( code.planet > 99)            /* Lagrange points */
            obuff[0] = obuff[1] = '-';
         if( 78 != strlen( obuff))
            fprintf( stderr, "'%s', len %d\n", obuff, (int)strlen( obuff));
         assert( 78 == strlen( obuff));
         if( code.prec1)         /* long. precision: blank unused digits */
            memset( obuff + 12 + code.prec1, ' ', 6 - code.prec1);
         if( code.prec2)         /* parallax data prec: blank unused */
            {
            memset( obuff + 43 + code.prec2, ' ', 7 - code.prec2);
            memset( obuff + 54 + code.prec2, ' ', 7 - code.prec2);
            }
         if( code.lon > 180.)
            code.lon -= 360.;
         if( google_map_links)        /* include HTML anchors */
            printf( "<a name=\"L%04d\"></a>", line_no++);
         show_link_for_this_line = (code.planet == 3 && google_map_links
                        && (code.lat || code.lon));
         if( show_link_for_this_line)
            {
            printf( "<a name=\"%s\"></a>", code.code);
            printf( "<a href=\"http://maps.google.com/maps?q=%f,%f\">",
                        code.lat, code.lon);
            }
         if( make_kml || show_link_for_this_line)
            _text_search_and_replace( buff, "&", "&amp;");
         if( make_kml && code.planet == 3 && (code.lat || code.lon))
            {
            i = 0;
            while( code.name[i] >= ' ')
               i++;
            printf( " <Placemark>\n");
            printf( "    <name>%.4s</name>\n", code.code);
            printf( "    <description>%.*s</description>\n", i, code.name);
            printf( "    <Point>\n");
            printf( "        <coordinates>%f,%f</coordinates>\n", code.lon, code.lat);
            printf( "    </Point>\n");
            printf( " </Placemark>\n\n");
            }
         if( !make_kml)
            {
            if( show_link_for_this_line)
               {
               char bing_link[100];

               snprintf( bing_link, sizeof( bing_link),
                     "<a href='https://www.bing.com/maps/?cp=%.7f~%.7f&amp;lvl=20&amp;style=a'>",
                              code.lat, code.lon);
               printf( "%.27s</a> %s%s %s</a>\n", obuff + google_offset,
                        bing_link, obuff + google_offset + 28, code.name);
               }
            else
               printf( "%s %s\n", obuff + google_offset, code.name);
            }
         }
      else if( dump_comments)    /* dump everything,  including */
         printf( "%s", buff);    /* comments from input file */
   fclose( ifile);
   if( !make_kml)
      printf( "%s\n", header + google_offset);
   else
      printf( "</Folder>\n</kml>\n");
   if( google_map_links)
      printf( "</pre></body></html>\n");
   return( 0);
}
#endif
