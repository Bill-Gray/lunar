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

#define N_EQUATORIAL_RADII 15

static const double equatorial_radii[N_EQUATORIAL_RADII] = {
      SUN_RADIUS, MERCURY_RADIUS, VENUS_RADIUS, EARTH_MAJOR_AXIS,
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

int get_mpc_code_info( mpc_code_t *cinfo, const char *buff)
{
   int i = 0, rval = -1;

   while( buff[i] > ' ' && buff[i] <= '~' && buff[i] != '!')
      i++;
   memset( cinfo, 0, sizeof( mpc_code_t));
   if( i >= 3 && i <= 4 && strlen( buff) >= 30)
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
      else if( buff[7] == '.' && strchr( "+- ", buff[21])
               && buff[14] == '.' && buff[23] == '.' && buff[3] == ' ')
         {                 /* 'standard' MPC format */
         cinfo->lon = atof( buff + 4);
         cinfo->rho_cos_phi = atof( buff + 13);
         cinfo->rho_sin_phi = atof( buff + 21);
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
   if( rval != -1)
      {
      cinfo->planet = rval;
      memcpy( cinfo->code, buff, 4);
      if( buff[3] == ' ')        /* standard 3-character code */
         cinfo->code[3] = '\0';
      else                       /* 'extended' 4-character code */
         cinfo->code[4] = '\0';
      if( cinfo->lon < 0.)
         cinfo->lon += PI + PI;
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
         double d_lon1 = atof( tbuff)      - lon_in_degrees;
         double d_lon2 = atof( tbuff + 20) - lon_in_degrees;
         const double d_lat1 = atof( tbuff + 10) - lat_in_degrees;
         const double d_lat2 = atof( tbuff + 30) - lat_in_degrees;

         while( d_lon1 > 180.)
            d_lon1 -= 360.;
         while( d_lon1 < -180.)
            d_lon1 += 360.;
         while( d_lon2 - d_lon1 > 180.)
            d_lon2 -= 360.;
         while( d_lon2 - d_lon1 < -180.)
            d_lon2 += 360.;
         if( d_lon1 * d_lon2 < 0. && d_lat1 * d_lat2 < 0.)
            {
            strcpy( buff, tbuff + 40);
            while( buff[i] >= ' ')
               i++;
            buff[i] = '\0';   /* remove trailing CR/LF */
            }
         }
   return( *buff ? 0 : -2);
}

#ifdef TEST_CODE

static int text_search_and_replace( char *str, const char *oldstr,
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

const char *kml_header_text =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
    "<Folder>\n"
    "      <name> MPC sites </name>\n"
    "      <description> MPC observatory sites </description>\n";

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( (argc < 2 ? "ObsCodes.htm" : argv[1]), "rb");
   char buff[200];
   mpc_code_t code;
   bool google_map_links = false, dump_comments = false;
   bool make_kml = false;
   int i;
   size_t google_offset = 0;

   const char *header =
            "Pl Code Longitude  Latitude     Altitude rho_cos   rho_sin_phi  region";

   if( !ifile)
      ifile = fopen( "ObsCodes.html", "rb");
   if( !ifile)
      {
      printf( "ObsCodes not opened\n");
      return( -1);
      }
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'v':
               dump_comments = true;
               break;
            case 'g':
               {
               const time_t t0 = time( NULL);
               FILE *hdr_ifile = fopen( "mpc_hdr.htm", "rb");

               assert( hdr_ifile);
               google_map_links = true;
               google_offset = 3;
               while( fgets( buff, sizeof( buff), hdr_ifile))
                  if( strstr( buff, "%s"))
                     printf( buff, ctime( &t0));
                  else if( *buff != '#')
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
   if( !make_kml)
      printf( "%s\n", header + google_offset);
   while( fgets( buff, sizeof( buff), ifile))
      if( get_mpc_code_info( &code, buff) != -1)
         {
         char region[100], obuff[200];
         bool show_link_for_this_line;

         code.lat *= 180. / PI;
         code.lon *= 180. / PI;
         *region = '\0';
         if( code.planet == 3)
            {
            FILE *ifile = fopen( "geo_rect.txt", "rb");

            if( ifile)
               {
               extract_region_data_for_lat_lon( ifile, region,
                    code.lat, code.lon);
               fclose( ifile);
               }
            }
#ifdef BITS_32
         sprintf( obuff,
#else
         snprintf( obuff, sizeof( obuff),
#endif
                 "%2d %-4s %10.6f %+10.6f %10.3f %9.7f %+10.7f %-15.15s ",
                  code.planet,
                  code.code, code.lon, code.lat,
                  code.alt, code.rho_cos_phi, code.rho_sin_phi,
                  region);
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
            printf( "<a name=\"L%04d\"></a>", i++);
         show_link_for_this_line = (code.planet == 3 && google_map_links
                        && (code.lat || code.lon));
         if( show_link_for_this_line)
            {
            printf( "<a name=\"%s\"></a>", code.code);
            printf( "<a href=\"http://maps.google.com/maps?q=%f,%f\">",
                        code.lat, code.lon);
            }
         if( make_kml || show_link_for_this_line)
            text_search_and_replace( buff, "&", "&amp;");
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
               printf( "%.27s</a> %s%s %s</a>", obuff + google_offset,
                        bing_link, obuff + google_offset + 28, code.name);
               }
            else
               printf( "%s %s", obuff + google_offset, code.name);
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
