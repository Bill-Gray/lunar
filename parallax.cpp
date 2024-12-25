#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "mpc_func.h"
#include "stringex.h"

/* Test program I run when I encounter a latitude/altitude and want to turn
it into parallax constants,  or vice versa.  See following 'error_exit'
function,  or run without command line arguments,  for usage.

   When compiled with CGI_VERSION #defined,  one gets a version suitable
for on-line use;  see https://www.projectpluto.com/parallax.htm for
an example of its usage.   */

#define EARTH_MAJOR_AXIS_IN_METERS    6378137.
     /* This code currently uses the GRS1980 value for the minor axis.
        That's about 0.105 mm less than the WGS1984 value.  I do not
        know of a case where the difference is actually measurable.  */
#define EARTH_MINOR_AXIS_IN_METERS    6356752.314140347
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

static char *show_angle( char *buff, double angle)
{
   const long microarcsec = (long)( fabs( angle) * 3600e+5);

   *buff = (angle < 0 ? '-' : '+');
   snprintf_err( buff + 1, 16, "%02ld %02ld %02ld.%05ld",
               (microarcsec / 360000000L),         /* degrees */
               (microarcsec / 6000000L) % 60L,     /* arcminutes */
               (microarcsec / 100000L) % 60L,     /* arcseconds */
                microarcsec % 100000L);           /* milliarcseconds */
   return( buff);
}

static int get_mpc_obscode_data( mpc_code_t *loc, const char *mpc_code, int pass)
{
   int rval = -1;
   char end_char = mpc_code[3];

   if( !end_char)
      end_char = ' ';
   while( rval && pass < 2)
      {
      const char *ifilename = (pass ? "ObsCodes.htm" : "rovers.txt");
      FILE *ifile;
      char buff[200];

#ifdef CGI_VERSION
      strlcpy_error( buff, "/home/projectp/public_html/cgi_bin/fo/");
#else
      strlcpy_error( buff, getenv( "HOME"));
      strlcat_error( buff, "/.find_orb/");
#endif
      strlcat_error( buff, ifilename);
      ifile = fopen( buff, "rb");
      if( !ifile)
         ifile = fopen( ifilename, "rb");
      if( !ifile)
         perror( buff);
      assert( ifile);
      while( rval && fgets( buff, sizeof( buff), ifile))
         if( !memcmp( buff, mpc_code, 3) && buff[3] == end_char)
            {
            const int err_code = get_mpc_code_info( loc, buff);
            char *tptr;

            if( 3 == err_code)
               {
               if( loc->lon > PI)
                  loc->lon -= PI + PI;
               printf( "%s !%+014.9f  %+013.9f %9.3f   %s%s", mpc_code,
                     loc->lon * 180. / PI,
                     loc->lat * 180. / PI,
                     loc->alt, loc->name, buff);
               }
            else if( -1 != err_code)
               printf( "%s\n", loc->name);
            tptr = strchr( buff, 10);
            assert( tptr);       /* remove trailing LF */
            *tptr = '\0';
#ifdef _WIN32                /* MS is different. */
            loc->name = _strdup( loc->name);
#else
            loc->name = strdup( loc->name);
#endif
            rval = 0;
            }
      fclose( ifile);
      pass++;
      }
   if( rval)
      {
      printf( "Couldn't find observatory code (%s)\n", mpc_code);
      exit( -1);
      }
   return( rval);
}

static void show_location( const mpc_code_t *loc)
{
   char buff[80];

   if( loc->lon)
      {
      printf( "Longitude %14.9f = %s\n", loc->lon, show_angle( buff, loc->lon));
      if( loc->lon < 0.)
         printf( "Longitude %14.9f = %s\n", loc->lon + 360., show_angle( buff, loc->lon + 360.));
      if( loc->lon > 180.)
         printf( "Longitude %14.9f = %s\n", loc->lon - 360., show_angle( buff, loc->lon - 360.));
      }
   printf( "Latitude  %11.9f = %s\n", loc->lat, show_angle( buff, loc->lat));
   printf( "Altitude %.5f meters above the WGS84 ellipsoid (_not_ above sea level/geoid)\n", loc->alt);
   printf( "Parallax constants %.11f %+.11f\n", loc->rho_cos_phi, loc->rho_sin_phi);
   printf( "In meters: %.5f %+.5f\n", loc->rho_cos_phi * EARTH_MAJOR_AXIS_IN_METERS,
                                      loc->rho_sin_phi * EARTH_MAJOR_AXIS_IN_METERS);
   if( loc->lon)
      {
      FILE *ifile = fopen( "geo_rect.txt", "rb");
      const double lon = loc->lon * PI / 180.;

#ifndef CGI_VERSION
      if( !ifile)
         {
         strlcpy_error( buff, getenv( "HOME"));
         strlcat_error( buff, "/.find_orb/geo_rect.txt");
         ifile = fopen( buff, "rb");
         }
#endif
      printf( "xyz in Earth radii %+.7f %+.7f %+.7f\n",
                                  cos( lon) * loc->rho_cos_phi,
                                  sin( lon) * loc->rho_cos_phi,
                                  loc->rho_sin_phi);
      printf( "xyz in meters      %+.5f %+.5f %+.5f\n",
                  cos( lon) * loc->rho_cos_phi * EARTH_MAJOR_AXIS_IN_METERS,
                  sin( lon) * loc->rho_cos_phi * EARTH_MAJOR_AXIS_IN_METERS,
                                   loc->rho_sin_phi * EARTH_MAJOR_AXIS_IN_METERS);
      if( ifile)
         {
         extract_region_data_for_lat_lon( ifile, buff, loc->lat, loc->lon);
         if( *buff)
            printf( "This point is somewhere in %s\n", buff + 2);
         fclose( ifile);
         }
      }
}

static double get_angle( const char *ibuff)
{
   bool is_negative = false;
   double rval;
   char field1[40], field2[40], field3[40];

   if( *ibuff == '-')
      {
      ibuff++;
      is_negative = true;
      }
   else if( *ibuff == '+')
      ibuff++;
   *field1 = *field2 = *field3 = '\0';
   sscanf( ibuff, "%20s %20s %20s", field1, field2, field3);
   rval = atof( field1) + atof( field2) / 60. + atof( field3) / 3600.;
   if( is_negative)
      rval = -rval;
   return( rval);
}

#ifndef CGI_VERSION

static void set_location_two_params( mpc_code_t *loc, double rho_cos_phi,
                                                 double rho_sin_phi)
{
   if( fabs( rho_cos_phi) > 2. || fabs( rho_sin_phi) > 2.)
      {             /* looks like parallax values in meters */
      rho_cos_phi /= EARTH_MAJOR_AXIS_IN_METERS;
      rho_sin_phi /= EARTH_MAJOR_AXIS_IN_METERS;
      }
   loc->rho_cos_phi = rho_cos_phi;
   loc->rho_sin_phi = rho_sin_phi;
   loc->lat = point_to_ellipse( 1., EARTH_MINOR_AXIS_IN_METERS / EARTH_MAJOR_AXIS_IN_METERS,
                 rho_cos_phi, rho_sin_phi, &loc->alt);
   loc->alt *= EARTH_MAJOR_AXIS_IN_METERS;
   loc->lon = 0.;
}

static void set_location_three_params( mpc_code_t *loc, double p1, double p2, double p3)
{
   if( fabs( p1) > 400 || fabs( p2) > 400)
      {             /* looks like parallax values in meters */
      p1 /= EARTH_MAJOR_AXIS_IN_METERS;
      p2 /= EARTH_MAJOR_AXIS_IN_METERS;
      p3 /= EARTH_MAJOR_AXIS_IN_METERS;
      }
   if( fabs( p1) < 2. && fabs( p2) < 2.)
      {
      const double rho_cos_phi = sqrt( p1 * p1 + p2 * p2);

      set_location_two_params( loc, rho_cos_phi, p3);
      loc->lon = atan2( p2, p1);
      }
   else        /* p1 = longitude, p2 = latitude, p3 = alt in meters */
      {
      loc->lon = p1 * PI / 180.;
      loc->lat = p2 * PI / 180.;
      loc->alt = p3;
      lat_alt_to_parallax( loc->lat, loc->alt, &loc->rho_cos_phi, &loc->rho_sin_phi,
               EARTH_MAJOR_AXIS_IN_METERS, EARTH_MINOR_AXIS_IN_METERS);
      }
}

static void error_exit( void)
{
   printf( "Run 'parallax' specifying a geographical location on the command line.\n"
           "This can be a comma-separated lat/lon/alt,  such as\n\n"
           "./parallax n44.02, w69.9, 381m\n"
           "./parallax 69 54W, n44 01 12, 1250 ft\n\n"
           "or similar (the above are all the same location),  or two arguments\n"
           "(rho_cos_phi, rho_sin_phi),  or just a three-characte MPC observatory code.\n"
           "The positions and parallax constants for that location will be shown.\n"
           "Run the program with two MPC codes as command-line arguments,  and the\n"
           "distance between them will be shown.\n");
   exit( -1);
}

int main( int argc, const char **argv)
{
   mpc_code_t loc;
   int i, j, n_commas = 0, use_only_obscodes_dot_html = 0;
   int show_mpc_format = 0;

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-' && !isdigit( argv[i][1]))
         {                    /* handle a command line option,  then remove it */
         switch( argv[i][1])
            {
            case 'f':
               use_only_obscodes_dot_html = 1;
               break;
            case 'm':
               show_mpc_format = 1;
               break;
            default:
               fprintf( stderr, "Unrecognized option '%s'\n", argv[i]);
               error_exit( );
            }
         argc--;
         memmove( argv + i, argv + i + 1, (argc - i) * sizeof( char *));
         i--;
         }
      else
         for( j = 0; argv[i][j]; j++)
            if( argv[i][j] == ',')
               n_commas++;

   if( 2 == n_commas)
      {
      char buff[100];

      strlcpy_error( buff, argv[1]);
      for( i = 2; i < argc; i++)
         {
         strlcat_error( buff, " ");
         strlcat_error( buff, argv[i]);
         }
      if( get_lat_lon_info( &loc, buff))
         error_exit( );
      }
   else if( argc < 2 || argc > 4)
      error_exit( );
   else if( argc == 2)          /* MPC code provided */
      get_mpc_obscode_data( &loc, argv[1], use_only_obscodes_dot_html);
   else if( argc == 3 && 3 == strlen( argv[1]) && 3 == strlen( argv[2]))
      {           /* two MPC codes provided */
      mpc_code_t loc2;
      double dist, posn_ang, p1[2], p2[2];

      get_mpc_obscode_data( &loc, argv[1], use_only_obscodes_dot_html);
      get_mpc_obscode_data( &loc2, argv[2], use_only_obscodes_dot_html);
      p1[0] = loc.lon;
      p1[1] = loc.lat;
      p2[0] = loc2.lon;
      p2[1] = loc2.lat;
      calc_dist_and_posn_ang( p1, p2, &dist, &posn_ang);
      printf( "(%s) is %.3f km from (%s),  at bearing %.2f (0=N, 90=E, 180=S, 270=W)\n",
               argv[2], dist * EARTH_MAJOR_AXIS_IN_METERS / 1000.,
               argv[1], 360. - posn_ang * 180. / PI);
      }
   else if( argc == 3)          /* parallax constants provided */
      set_location_two_params( &loc, atof( argv[1]), atof( argv[2]));
   else        /* argc == 4 */
      set_location_three_params( &loc, get_angle( argv[1]), get_angle( argv[2]), atof( argv[3]));
   loc.lat *= 180. / PI;
   loc.lon *= 180. / PI;
   show_location( &loc);
   if( show_mpc_format)             /* output in format desired when */
      {                             /* submitting corrections to MPC */
      printf( "contact_name: Bill Gray\n");
      printf( "email_adr: pluto@projectpluto.com\n");
      printf( "observatory_code: %s\n", loc.code);
      printf( "observatory_name: %s\n", loc.name);
      printf( "observatory_site: %s\n", loc.name);
      printf( "observatory_country: US\n");
      printf( "observatory_lat: %f\n", loc.lat);
      printf( "observatory_long: %f\n", loc.lon);
      printf( "observatory_alt: %f\n", loc.alt);
      printf( "telescope_height: 0\n");
      printf( "reference: WGS84\n");
      }
}
#else

#ifdef __has_include
   #if __has_include("cgi_func.h")
       #include "cgi_func.h"
   #else
       #error   \
         'cgi_func.h' not found.  This project depends on the 'lunar'\
         library.  See www.github.com/Bill-Gray/lunar .\
         Clone that repository,  'make'  and 'make install' it.
       #ifdef __GNUC__
         #include <stop_compiling_here>
            /* Above line suppresses cascading errors. */
       #endif
   #endif
#else
   #include "cgi_func.h"
#endif
#include <string.h>

int main( void)
{
   mpc_code_t loc;
   char field[30], buff[100];
   int rval;
   double xyz[3];

   printf( "Content-type: text/html\n\n");
   printf( "<pre>");
   avoid_runaway_process( 300);
   rval = initialize_cgi_reading( );
   if( rval <= 0)
      {
      printf( "<p> <b> CGI data reading failed : error %d </b>", rval);
      printf( "This isn't supposed to happen.</p>\n");
      return( 0);
      }
   memset( &loc, 0, sizeof( mpc_code_t));
   xyz[0] = xyz[1] = xyz[2] = 0.;
   while( !get_cgi_data( field, buff, NULL, sizeof( buff)))
      {
      if( !strcmp( field, "rho_sin_phi"))
         loc.rho_sin_phi = atof( buff);
      else if( !strcmp( field, "rho_cos_phi"))
         loc.rho_cos_phi = atof( buff);
      else if( !strcmp( field, "lat") && strchr( buff, '.'))
         loc.lat = get_angle( buff) * PI / 180.;
      else if( !strcmp( field, "lon") && strchr( buff, '.'))
         loc.lon = get_angle( buff) * PI / 180.;
      else if( !strcmp( field, "alt"))
         loc.alt = atof( buff) / EARTH_MAJOR_AXIS_IN_METERS;
      else if( !memcmp( field, "xyz", 3) && field[3] >= '0' && field[3] <= '2')
         xyz[field[3] - '0'] = atof( buff);
      else if( !strcmp( field, "mpc_code") && 3 == strlen( buff))
         get_mpc_obscode_data( &loc, buff, 0);
      }
   if( xyz[0] || xyz[1] || xyz[2])
      {
      loc.lon = atan2( xyz[1], xyz[0]);
      loc.rho_cos_phi = sqrt( xyz[0] * xyz[0] + xyz[1] * xyz[1]);
      loc.rho_sin_phi = xyz[2];
      }
   if( loc.rho_sin_phi && loc.rho_cos_phi)
      {
      if( fabs( loc.rho_sin_phi) > 100. || fabs( loc.rho_cos_phi) > 100.)
         {           /* assume meters were entered */
         loc.rho_sin_phi /= EARTH_MAJOR_AXIS_IN_METERS;
         loc.rho_cos_phi /= EARTH_MAJOR_AXIS_IN_METERS;
         }
      loc.lat = point_to_ellipse( 1., EARTH_MINOR_AXIS_IN_METERS / EARTH_MAJOR_AXIS_IN_METERS,
                 loc.rho_cos_phi, loc.rho_sin_phi, &loc.alt);
      }
   else if( loc.lat)
      lat_alt_to_parallax( loc.lat, loc.alt, &loc.rho_cos_phi, &loc.rho_sin_phi,
              1., EARTH_MINOR_AXIS_IN_METERS / EARTH_MAJOR_AXIS_IN_METERS);
   else
      {
      printf( "Must provide either rho sin(phi) and rho cos(phi),  in\n"
              "which case the latitude and altitude will be computed and\n"
              "shown;  <i>or</i> latitude and altitude,  in which case you'll\n"
              "get the parallax constants as output.  Or,  you can provide\n"
              "x, y, and z;  or an MPC code.  Hit the Back-arrow in your\n"
              "browser and review your options.\n");
      return( 0);
      }
   loc.lat *= 180. / PI;
   loc.lon *= 180. / PI;
   loc.alt *= EARTH_MAJOR_AXIS_IN_METERS;
   show_location( &loc);
   if( loc.lon)
      {
      if( loc.lon > 180.)
         loc.lon -= 360.;
      printf( "<a href='http://maps.google.com/maps?q=%.7f,%.7f'>",
                        loc.lat, loc.lon);
      printf( "Click here for a G__gle map of this location</a>\n");
      printf( "<a href='https://www.bing.com/maps/?cp=%.7f~%.7f&lvl=18&style=a'>",
                        loc.lat, loc.lon);
      printf( "Click here for a Bing map of this location</a>\n");
      }
   return( 0);
}
#endif
