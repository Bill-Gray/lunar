/* astcheck.cpp: an off-line 'MPChecker'

Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifdef _WIN32
   #define WIN32_LEAN_AND_MEAN
   #include <windows.h>
#else
   #include <sys/types.h>
   #include <unistd.h>
#endif
#include "watdefs.h"
#include "date.h"
#include "comets.h"
#include "afuncs.h"
#include "mpc_func.h"
#include "stringex.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define LOG_10 2.3025850929940456840179914546843642076011014886287729760333279009675726
const double radians_to_arcsec = 180. * 3600. / PI;

int get_earth_loc( const double t_millennia, double *results);

#if defined( __WATCOMC__) && !defined( _WIN32)
void usleep( const long microseconds)
{
   INTENTIONALLY_UNUSED_PARAMETER( microseconds);
}

int getpid( void)
{
   return( 1);
}
#endif

static int get_mpc_data( const char *buff, double *jd, double *ra, double *dec)
{
   int err_code;

   if( strlen( buff) < 80 || strlen( buff) > 83)
      err_code = -5;
   else
      err_code = get_ra_dec_from_mpc_report( buff, NULL, ra, NULL,
                                                NULL, dec, NULL);
   *jd = extract_date_from_mpc_report( buff, NULL);
   return( err_code);
}

static inline double law_of_cosines( const double a, const double b, const double c)
{
   return( .5 * (a * a + b * b - c * c) / (a * b));
}

static double calc_obs_magnitude( ELEMENTS *elem, const double obj_sun,
                      const double obj_earth, const double earth_sun)
{
   double magnitude;

   if( !elem->abs_mag)
      magnitude = 0.;
   else if( !elem->is_asteroid)
      magnitude = elem->slope_param * log10( obj_sun);
   else
      {
      const double cos_phase_ang =
                  law_of_cosines( obj_sun, obj_earth, earth_sun);

      magnitude = 5. * log10( obj_sun) + phase_angle_correction_to_magnitude(
                     acose( cos_phase_ang), elem->slope_param);
      }
   magnitude += 5. * log10( obj_earth) + elem->abs_mag;
   return( magnitude);
}

static int16_t integerize_angle( double angle)
{
   const int iangle = (int)floor( angle * 32768. / PI + .5);

   return( (int16_t)iangle);
}

#define AST_DATA struct AST_DATA

#pragma pack(2)
AST_DATA
   {
   int16_t ra, dec;                /* in 360/65536 degrees */
   };
#pragma pack( )

#define MAX_SOF_SIZE 200

static FILE *orbits_file;
static int n_asteroids, record_length;
int verbose = 0;
const char *data_path = NULL;
char sof_header[MAX_SOF_SIZE];
int32_t sof_checksum;

static FILE *get_file_from_path( const char *filename, const char *permits)
{
   FILE *fp = NULL;

   if( data_path && *data_path)
      {
      char buff[450];

      strlcpy_error( buff, data_path);
      if( buff[strlen( buff) - 1] != '/')
         strlcat_error( buff, "/");
      strlcat_error( buff, filename);
      fp = fopen( buff, permits);
      }
   if( !fp)
      fp = fopen( filename, permits);
   return( fp);
}

static FILE *get_sof_file( const char *filename)
{
   FILE *ifile = get_file_from_path( filename, "rb");

   if( ifile)
      {
      int filelen;
      size_t i, j;
      char buff[450];

      if( !fgets( buff, sizeof( buff), ifile))
         {
         fprintf( stderr, "Unable to read '%s'\n", filename);
         exit( -4);
         }
      record_length = (int)strlen( buff);
      assert( record_length < MAX_SOF_SIZE);
      strlcpy_error( sof_header, buff);
      fseek( ifile, 0L, SEEK_END);
      filelen = ftell( ifile);
      if( filelen % record_length)
         {
         printf( "'%s' appears to be corrupted.\n", filename);
         exit( -5);
         }
      n_asteroids = filelen / record_length - 1;      /* there's a header line */
      for( i = 0; i < 4; i++)
         {
         const int32_t big_prime = 1234567891;

         fseek( ifile, (long)( i * (filelen - sizeof( buff))) / 3L, SEEK_SET);
         if( fread( buff, sizeof( buff), 1, ifile))
            for( j = 0; j < sizeof( buff); j++)
               sof_checksum = sof_checksum * big_prime + (int32_t)buff[i];
         }
      fseek( ifile, 0L, SEEK_SET);
      if( verbose)
         printf( "'%s': %d objects; record size %d\n",
                      filename, n_asteroids, record_length);
      }
   return( ifile);
}

static double obj_sun_dist;

static double compute_asteroid_loc( const double *earth_loc,
            ELEMENTS *elem, const double jd, double *ra, double *dec)
{
   double r1 = 0., dist, asteroid_loc[4];
   int i, n_iterations = 0;

   do             /* light-time lag:  should converge _very_ fast       */
      {           /* it'll almost always require exactly two iterations */
      dist = r1;
      comet_posn( elem, jd - dist / AU_PER_DAY, asteroid_loc);
      for( i = 0; i < 3; i++)
         asteroid_loc[i] -= earth_loc[i];
      r1 = vector3_length( asteroid_loc);
      if( n_iterations >= 15)
         fprintf( stderr, "? e=%f, q=%f, jd=%f; diff %f/%f\n", elem->ecc, elem->q, jd, dist - r1, r1);
      assert( n_iterations++ < 20);
      }
      while( fabs( dist - r1) > .001);
   ecliptic_to_equatorial( asteroid_loc);
   *ra = atan2( asteroid_loc[1], asteroid_loc[0]);
   *dec = asin( asteroid_loc[2] / r1);
   obj_sun_dist = asteroid_loc[3];
   return( r1);
}

AST_DATA *compute_day_data( const long ijd)
{
   char tbuff[300];
   int i, counter = 0;
   AST_DATA *rval;
   const double jd = (double)ijd;
   double earth_loc[6];
   clock_t t0 = clock( );

   if( verbose)
      printf( "Computing data for %ld (%d asteroids)\n", ijd, n_asteroids);
   rval = (AST_DATA *)malloc( n_asteroids * sizeof( AST_DATA));
   if( !rval)
      {
      printf( "OUT OF MEMORY\n");
      return( NULL);
      }
   get_earth_loc( (jd      - 2451545.) / 365250., earth_loc);
   fseek( orbits_file, (long)record_length, SEEK_SET);
   tbuff[record_length] = '\0';
   for( i = 0; i < n_asteroids && fgets( tbuff, sizeof( tbuff), orbits_file); i++)
      {
      ELEMENTS class_elem;
      double ra, dec;

      if( extract_sof_data( &class_elem, tbuff, sof_header))
         {
         printf( "'mpcorb.sof' data doesn't parse correctly:\n%s\n", tbuff);
         exit( -1);
         }

      compute_asteroid_loc( earth_loc, &class_elem, jd, &ra, &dec);

      rval[i].ra = integerize_angle( ra);
      rval[i].dec = integerize_angle( dec);
      if( verbose && counter <= i * 80 / n_asteroids)
         {
         printf( "%d", counter % 10);
         counter++;
         }
      }
   if( verbose)
      printf( "\nTime: %.1f seconds\n",
                  (clock( ) - t0) / (double)CLOCKS_PER_SEC);
   return( rval);
}

static double centralize_angle( double ang)
{
   while( ang > PI)
      ang -= PI + PI;
   while( ang < -PI)
      ang += PI + PI;
   return( ang);
}

/* In loading up the 'day data' (basically a low-precision RA/dec for
the geocentric position of each asteroid as of a certain day),  we
attempt to open a file for that day of the form YYYYMMDD.chk.  If
the file is opened,  and the header indicates the correct version,
number of asteroids,  and checksum, we just load up the data
from the file and return.   Otherwise,  this function creates the
'day data' file,  saves it for future use,  and returns the data
it's just generated.

It used to be as simple as that,  but 'astcheck' now gets used in a
mode where many instances are run in parallel.  Run seventeen
instances for the same day, and all seventeen would start generating
.chk files.  Ideally, _one_ would do so and the other sixteen would
wait for the file to be ready.  So accessing and building of .chk
files now works as follows :

(1) We try to open the .chk day data file we need and read in the data.
The file may not exist,  or may contain outdated data (the 'sof_checksum'
or other header data may not match what we're looking for).  If that's
the case,  we drop through the code to the point where we generate the
day data and create the day data file.  Right away,  we write out the
header information and fflush() it.  That way,  the other sixteen
processes will see a header with no data (yet) and realize that they
should wait a bit for the data to become available.

(2) If we read in the header and it's what we're looking for (its
"magic number",  checksum,  and the number of asteroids match),
we try to read in the actual asteroid data.  That part may simply
succeed (because the .chk file in question has already been completed).
If we only get none or some of the data,  we assume we'll have to
wait for some other process to finish creating and writing the file.

In that case,  we 'sleep' for one second and try reading the data
again.  Ideally,  we eventually get all the data.  Just as a check,
if this takes more than sixty seconds, we assume something has gone
seriously wrong and bomb out;  this is to avoid having phantom
processes left hanging.  (Haven't seen that actually happen yet.) We
also close the input file and reopen it and fseek() to the start of
the remaining data.  That's just to force an update to the input
file size.         */

#define HEADER_SIZE 4

static AST_DATA *get_cached_day_data( const int ijd)
{
   AST_DATA *rval = NULL;
   char filename[20];
   FILE *ifile, *ofile;
   int32_t header[HEADER_SIZE];
   const int32_t magic_version_number = 1314159266;

                  /* Create a filename in 'YYYYMMDD.chk' form: */
   full_ctime( filename, (double)ijd, FULL_CTIME_YMD | FULL_CTIME_NO_SPACES
                     | FULL_CTIME_DATE_ONLY | FULL_CTIME_MONTHS_AS_DIGITS
                     | FULL_CTIME_LEADING_ZEROES);
   strlcat_error( filename, ".chk");
   if( verbose > 2)
      printf( "Creating '%s'\n", filename);
   ifile = get_file_from_path( filename, "rb");
   if( ifile)
      {
      if( !fread( header, HEADER_SIZE, sizeof( int), ifile))
         {
         printf( "Error reading header data in '%s'\n", filename);
         exit( -2);
         }
      if( header[0] == magic_version_number &&
                      header[1] == sof_checksum && header[2] == n_asteroids)
         {
         fseek( ifile, 0L, SEEK_END);
         if( ftell( ifile) == (n_asteroids + HEADER_SIZE) * 4L)
            rval = (AST_DATA *)malloc( n_asteroids * sizeof( AST_DATA));
         fseek( ifile, (long)sizeof( header), SEEK_SET);
         }
      if( rval)   /* appears to be legitimate cached data */
         {
         int n_read = 0, n_iterations = 0;

         while( n_read != n_asteroids)
            {
            n_read += (int)fread( rval + n_read, sizeof( AST_DATA),
                                        n_asteroids - n_read, ifile);
            if( n_read != n_asteroids)
               {
               if( n_iterations++ == 60)
                  {               /* Give up after one minute */
                  fprintf( stderr, "Overtime\n");
                  exit( -3);
                  }
#ifndef _WIN32
               if( verbose)
                  printf( "(%d) %d of %d read, iter %d\n", (int)getpid( ),
                        n_read, n_asteroids, n_iterations);
#endif
                                 /* Close/reopen file.  Gets around some */
                                 /* caching issues.   */
               fclose( ifile);
#ifndef _WIN32
               usleep( 1000000);    /* 1000000 microseconds = 1 second */
#else
               Sleep( 1000);        /* 1000 milliseconds = 1 second */
#endif
               ifile = get_file_from_path( filename, "rb");
               fseek( ifile, 16 + 4L * n_read, SEEK_SET);
               }
            }
         fclose( ifile);
         return( rval);
         }
      fclose( ifile);
      }
   header[0] = magic_version_number;
   header[1] = sof_checksum;
   header[2] = n_asteroids;
   header[3] = -1;         /* not currently used */
   ofile = get_file_from_path( filename, "wb");
   if( !ofile)
      {
      fprintf( stderr, "Couldn't open '%s'\n", filename);
      perror( "File open failure");
      exit( -1);
      }
   fwrite( header, HEADER_SIZE, sizeof( int), ofile);
   fflush( ofile);
   rval = compute_day_data( ijd);
   fwrite( rval, n_asteroids, sizeof( AST_DATA), ofile);
   fclose( ofile);

   return( rval);
}

         /* Figuring out if an RA 'x' is between two RAs 'bound1' */
         /* and 'bound2' is complicated by the discontinuity at   */
         /* +/-180 degrees (or +/-32768 of the angular units used */
         /* in this program.)                                     */
static int is_between( int bound1, int bound2, int x, int tolerance)
{
   int rval;

   bound1 -= x;
   bound2 -= x;
   if( bound1 - bound2 > 32768)     /* bound1 is more than 180 degrees greater than bound2 */
      bound1 -= 65536;
   if( bound2 - bound1 > 32768)    /* bound2 is more than 180 degrees greater than bound1 */
      bound2 -= 65536;
   while( bound1 + bound2 < -65536)
      {
      bound1 += 65536;
      bound2 += 65536;
      }
   while( bound1 + bound2 > 65536)
      {
      bound1 -= 65536;
      bound2 -= 65536;
      }
   if( bound1 < bound2)
      rval = (0 >= bound1 - tolerance && 0 <= bound2 + tolerance);
   else
      rval = (0 >= bound2 - tolerance && 0 <= bound1 + tolerance);
   return( rval);
}

int qsort_mpc_cmp( const void *elem1, const void *elem2)
{
   const char **buff1 = (const char **)elem1;
   const char **buff2 = (const char **)elem2;
   int compare = memcmp( *buff1, *buff2, 12);

   if( !compare)     /* same ID;  now compare times */
      compare = memcmp( (*buff1) + 15, (*buff2) + 15, 16);
   return( compare);
}

   /* To compute observed object motion,  this code grabs the RA/dec */
   /* of the first observation;  then looks ahead for observations   */
   /* of the same object from the same station.  The difference      */
   /* between the pair gives motion.  Except we want to use a pair   */
   /* far enough apart to show real motion,  but not so far apart    */
   /* that sky curvature or acceleration are factors.                */
   /*    Return value is the JD of the observation used to determine */
   /* the motion.  If there was only one observation,  we give the   */
   /* time of that observation plus epsilon,  to evade division by   */
   /* zero.                                                          */

static double compute_motion( const char **lines, const int n_lines,
                           double *ra_motion, double *dec_motion)
{
   double ra, dec, jd, rval;
   int i;

   *ra_motion = *dec_motion = 0.;
   get_mpc_data( lines[0], &jd, &ra, &dec);
   rval = jd + 1e-6;
   for( i = 1; i < n_lines && !memcmp( lines[0], lines[i], 12); i++)
      if( !memcmp( lines[0] + 77, lines[i] + 77, 3))
         {
         double ra1, dec1, jd1;
         const double five_degrees = PI / 36.;

         if( !get_mpc_data( lines[i], &jd1, &ra1, &dec1)
                && fabs( jd1 - jd) < 10.   /* within ten days... */
                && fabs( ra1 - ra) < five_degrees
                && fabs( dec1 - dec) < five_degrees)
            {
            *ra_motion = (ra1 - ra) / (jd1 - jd);
            *dec_motion = (dec1 - dec) / (jd1 - jd);
            *ra_motion *= cos( dec);
                     /* cvt radians/day to arcsec/hr: */
            *ra_motion *= radians_to_arcsec  / 24.;
            *dec_motion *= radians_to_arcsec / 24.;
            rval = jd1;
            }
         }
   return( rval);
}

void observer_cartesian_coords( const double jd, const double lon,
              const double rho_cos_phi, const double rho_sin_phi,
              double *vect)
{
   const double angle = lon + green_sidereal_time( jd);
   const double earth_major_axis = 6378140. / AU_IN_METERS;

   *vect++ = cos( angle) * rho_cos_phi * earth_major_axis;
   *vect++ = sin( angle) * rho_cos_phi * earth_major_axis;
   *vect++ = rho_sin_phi               * earth_major_axis;
}

static double get_topo_loc( const double jd, double *topo_loc, const double lon,
              const double rho_cos_phi, const double rho_sin_phi)
{
   double earth_loc[6];
   int i;

   get_earth_loc( (jd - 2451545.) / 365250., earth_loc);
   observer_cartesian_coords( jd, lon, rho_cos_phi, rho_sin_phi, topo_loc);
   equatorial_to_ecliptic( topo_loc);
   for( i = 0; i < 3; i++)
      topo_loc[i] += earth_loc[i];
   return( vector3_length( topo_loc));
}

static int get_mpcorb_dot_dat_line( const char *filename, const int line_no,
                                 char *buff)
{
   static long offset, line_len;
   FILE *ifile = get_file_from_path( filename, "rb");
   int rval = 0, n_iterations = 3;

   if( !ifile)
      return( -1);
   if( !line_len)
      {
      while( fgets( buff, 210, ifile) && memcmp( buff, "00001 ", 6))
         ;
      line_len = (long)strlen( buff);
      offset = ftell( ifile) - line_len;
      }
   *buff = '\0';
   if( fseek( ifile, offset + (long)line_no * line_len, SEEK_SET))
      rval = -1;
   else while( n_iterations && strlen( buff) != (size_t)line_len
                         && fgets( buff, 210, ifile))
      n_iterations--;
   if( !rval && strlen( buff) != (size_t)line_len)
      rval = -2;
   fclose( ifile);
   return( rval);
}

static char _dummy_filename[40];

static void make_fake_file( const char **argv)
{
   char buff[81];
   const double jd = get_time_from_string( 0., argv[2], 0, NULL);
   const double ra = atof( argv[3]);
   const double dec = atof( argv[4]);
   FILE *ofile = fopen( _dummy_filename, "wb");

   strlcpy_error( buff, "    dummy     C");
   snprintf_append( buff, 56, "%16.8f %011.7f %+011.7f", jd, ra, dec);
   strlcat_error( buff, "                 Synth");
   strlcat_error( buff, argv[5]);          /* MPC code */
   fprintf( ofile, "%s\n", buff);
   fclose( ofile);
   assert( strlen( buff) == 80);
}

static void err_message( void)
{
   printf( "\nastcheck needs the name of a file containing MPC-formatted (80-column)\n");
   printf( "astrometric data.  It will then attempt to match the records to known\n");
   printf( "objects in the 'mpcorb.sof' file.\n\n");
   printf( "Command-line options are:\n\n");
   printf( "   -r(dist)   Set search distance to 'dist' arcseconds. Default is 3600.\n");
   printf( "   -z(tol)    Set motion match tolerance to 'tol' arcsec/hr. Default is 10.\n");
   printf( "   -m(mag)    Set limiting mag to 'mag'.  Default is 22.\n");
   printf( "   -l         Show distance from line of variations. Experimental.\n");
   printf( "Alternatively,  one can get a list of asteroids/comets within a desired\n");
   printf( "area with\n\n");
   printf( "astcheck -c (date) (RA in degrees) (dec in degrees) (MPC code) (options)\n\n");
   printf( "For example, 'astcheck -c 2022apr3.1415 292.653 -7.653 E12 -r7200' would get\n");
   printf( "a list of asteroids within two degrees of that RA/dec as seen from (E12).\n");
}

static void show_astcheck_info( void)
{
   printf( "ASTCHECK version %s %s\n", __DATE__, __TIME__);
   printf( "%d objects\n", n_asteroids);
}

/* An oversimplified getopt(). */

static const char *get_arg( const int argc, const char **argv, const int idx)
{
   if( argv[idx][2] || idx == argc - 1)
      return( argv[idx] + 2);
   else
      return( argv[idx + 1]);
}

#define IS_POWER_OF_TWO( n)    (((n) & ((n)-1)) == 0)

#if defined(_MSC_VER) && _MSC_VER < 1900
                      /* For older MSVCs,  we have to supply our own  */
                      /* snprintf().  See snprintf.cpp for details.  */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

#ifdef CGI_VERSION
int astcheck_main( const int argc, const char **argv)
#else
int main( const int argc, const char **argv)
#endif
{
   double jd, ra, dec;
   FILE *ifile;
   const char *sof_filename = "mpcorb.sof";
   char buff[400];
   char **ilines = NULL;
   int show_lov = 0;
   int i, n_ilines = 0, n, max_results = 100;
   int n_lines_printed = 0;
   double tolerance_in_arcsec = 18000.;       /* = five degrees */
   double mag_limit = 22.;
   AST_DATA *day_data[2] = { NULL, NULL};
   long curr_loaded_day_data = 0;
   FILE *mpc_station_file;
   char curr_station[7];
   double rho_sin_phi = 0., rho_cos_phi = 0., longitude = 0.;
   double motion_tolerance = 10.;  /* require a match to within 10"/hr */
   int results_array_size = 5;
   char **results = (char **)calloc( results_array_size, sizeof( char *));
   bool is_list_file = false;
   bool show_header = true;
   const char *mpcorb_extracts = "";
   void *ades_context = init_ades2mpc( );

   if( argc < 2)
      {
      err_message( );
      return( -1);
      }
#ifdef _WIN32
   strlcpy_error( _dummy_filename, "astcheck.tmp");
#else
   snprintf( _dummy_filename, sizeof( _dummy_filename), "astcheck%d.tmp",
                        (int)getpid( ));
#endif
   if( !strcmp( argv[1], "-c"))
      {
      assert( argc > 5);
      make_fake_file( argv);
      is_list_file = true;
      max_results = 20000;
      }
   memset( curr_station, 0, sizeof( curr_station));
   for( i = (is_list_file ? 6 : 2); i < argc; i++)
      if( argv[i][0] == '-')
         {
         const char *arg = get_arg( argc, argv, i);

         assert( arg);
         switch( argv[i][1])
            {
            case 'e':
               mpcorb_extracts = arg;
               break;
            case 'h':
               show_header = false;
               break;
            case 'r':
               tolerance_in_arcsec = atof( arg);
               break;
            case 'v':
               setvbuf( stdout, NULL, _IONBF, 0);
               verbose = 1 + atoi( arg);
               break;
            case 'm':
               mag_limit = atof( arg);
               break;
            case 'p':
               data_path = arg;
               break;
            case 'l':
               show_lov = 1;
               break;
            case 'z':
               motion_tolerance = atof( arg);
               break;
            case 'M':
               max_results = atoi( arg);
               break;
#ifdef NOT_READY_QUITE_YET
            case 'e':
               show_uncertainty = 1;
               break;
#endif
            case 'f':
               sof_filename = arg;
               break;
            default:
               printf( "%s: unrecognized command-line option\n", argv[i]);
               break;
            }
         }
   mpc_station_file = get_file_from_path( "ObsCodes.html", "rb");
   if( !mpc_station_file)        /* perhaps stored with truncated extension? */
      mpc_station_file = get_file_from_path( "ObsCodes.htm", "rb");
   if( !mpc_station_file)
      {
      printf( "ObsCodes.html not found; parallax won't be included!\n");
      printf( "Astcheck can run without this file,  but will produce better\n");
      printf( "results if it has it :\n\n");
      printf( "https://www.minorplanetcenter.net/iau/lists/ObsCodes.html\n\n");
      printf( "Download this file and put it in the directory in which\n");
      printf( "astcheck is running.\n");
      }
   if( argc < 2)
      {
      printf( "No input file name specified\n");
      err_message( );
      return( -1);
      }

   orbits_file = get_sof_file( sof_filename);
   if( !orbits_file)
      {
      printf( "Couldn't open '%s'\n", sof_filename);
      printf( "Astcheck gets orbital elements from 'mpcorb.sof'.  See\n"
              "https://www.projectpluto.com/astcheck.htm#setup for details on\n"
              "how to create/maintain that file.\n");
      return( -2);
      }

   ifile = fopen( is_list_file ? _dummy_filename : argv[1], "rb");
   if( !ifile)
      {
      printf( "%s not opened\n", argv[1]);
      err_message( );
      return( -3);
      }

               /* Read astrometry lines and allocate memory for them : */
   ilines = (char **)malloc( (n_ilines + 1) * sizeof( char *));
   n_ilines = 0;
   while( fgets_with_ades_xlation( buff, sizeof( buff), ades_context, ifile))
      if( !get_mpc_data( buff, &jd, &ra, &dec))
         {
         n_ilines++;
         if( IS_POWER_OF_TWO( n_ilines))
            ilines = (char **)realloc( ilines, n_ilines * 2 * sizeof( char *));
         ilines[n_ilines - 1] = (char *)malloc( strlen( buff) + 1);
         strcpy( ilines[n_ilines - 1], buff);
         }
   fclose( ifile);
   free_ades2mpc_context( ades_context);
   if( is_list_file)
#ifdef _WIN32                /* MS is different. */
      _unlink( _dummy_filename);
#else
      unlink( _dummy_filename);
#endif
   if( !n_ilines)
      {
      printf( "No astrometry found in '%s'\n", argv[1]);
      err_message( );
      return( -1);
      }
   qsort( ilines, n_ilines, sizeof( char **), qsort_mpc_cmp);
   for( n = 0; n < n_ilines; n++)
      if( strlen( ilines[n]) >= 80
                     && (!n || memcmp( ilines[n], ilines[n - 1], 12))
                     && !get_mpc_data( ilines[n], &jd, &ra, &dec))
         {
         double earth_loc[6], earth_loc2[6];
         double ra_motion = 0., dec_motion = 0., earth_sun_dist;
         const double cos_dec = cos( dec);
         const int16_t int_dec = (int16_t)integerize_angle( dec);
         const int16_t int_ra  = (int16_t)integerize_angle( ra);
         const double delta_t = td_minus_ut( jd) / seconds_per_day;
         double jd2;
         const int16_t tolerance = (int16_t)
                          ( tolerance_in_arcsec * 65536 / (360. * 3600.));
         char tbuff[300];
         int n_results = 0;
         int n_checked = 0;
         bool singleton_observation;

         jd += delta_t;
         if( mpc_station_file && memcmp( ilines[n] + 77, curr_station, 3))
            {
            int got_station_data = 0;

            strlcpy_error( curr_station, ilines[n] + 77);
            curr_station[3] = '\0';
            fseek( mpc_station_file, 0L, SEEK_SET);
            rho_sin_phi = rho_cos_phi = longitude = 0.;
            while( !got_station_data &&
                           fgets( tbuff, sizeof( tbuff), mpc_station_file))
               got_station_data = !memcmp( tbuff, curr_station, 3);
            if( !got_station_data)
               {
               FILE *rovers_file = get_file_from_path( "rovers.txt", "rb");

               if( !rovers_file)
                  {
                  fprintf( stderr, "Couldn't open 'rovers.txt'\n");
                  return( -1);
                  }
               while( !got_station_data &&
                           fgets( tbuff, sizeof( tbuff), rovers_file))
                  got_station_data = !memcmp( tbuff, curr_station, 3);
               fclose( rovers_file);
               }
            if( got_station_data)
               {
               mpc_code_t code_info;
               const int err_code = get_mpc_code_info( &code_info, tbuff);

               if( err_code < 0)
                  {
                  fprintf( stderr, "Code '%s' not found; error %d\n", tbuff, err_code);
                  got_station_data = 0;
                  }
               else
                  {
                  longitude = code_info.lon * 180. / PI;
                  rho_cos_phi= code_info.rho_cos_phi;
                  rho_sin_phi= code_info.rho_sin_phi;
                  }
               }
            if( !got_station_data)
               printf( "FAILED to find MPC code %s\n", curr_station);
            longitude *= PI / 180.;
            }

         if( curr_loaded_day_data != (int)jd || !day_data[0] || !day_data[1])
            {
            if( day_data[0])
               free( day_data[0]);
            if( day_data[1])
               free( day_data[1]);
            day_data[0] = get_cached_day_data( (int)jd);
            day_data[1] = get_cached_day_data( (int)jd + 1);
            curr_loaded_day_data = (int)jd;
            }
         if( !n && show_header)     /* on our very first object: */
            {
            show_astcheck_info( );
            printf( "An explanation of these data is given at the bottom of the list.\n");
            if( is_list_file)
               {
               printf( "                           RA  (J2000)  dec       mag ");
               printf( "  dRA/dt    dDec/dt\n");
               }
            else
               printf( "                             d_ra   d_dec    dist    mag  motion \n");
            }
         if( verbose)
            printf( "JD %f, RA %f, dec %f\n",
                     jd, ra * 180. / PI, dec * 180. / PI);
         jd2 = compute_motion( (const char **)ilines + n, n_ilines - n,
                                 &ra_motion, &dec_motion);
         jd2 += delta_t;
         earth_sun_dist =
               get_topo_loc( jd, earth_loc, longitude, rho_cos_phi, rho_sin_phi);
         get_topo_loc( jd2, earth_loc2, longitude, rho_cos_phi, rho_sin_phi);
         memcpy( buff, ilines[n], 12);
         buff[12] = '\0';
         singleton_observation = ( !ra_motion && !dec_motion);
         if( singleton_observation)
            {
            if( !is_list_file)
               printf( "\n%s: only one observation\n", buff);
            }
         else
#ifdef CGI_VERSION
            printf( "\n<b>%s: %.0f\"/hr in RA, %.0f\"/hr in dec (%.2f hours)</b>\n",
#else
            printf( "\n%s: %.0f\"/hr in RA, %.0f\"/hr in dec (%.2f hours)\n",
#endif
                        buff, ra_motion, dec_motion, (jd2 - jd) * 24.);
         n_lines_printed++;
         for( i = 0; i < n_asteroids; i++)
            {
            const int16_t tolerance2 = tolerance;

            assert( day_data[0]);
            assert( day_data[1]);
            if( is_between( day_data[0][i].ra, day_data[1][i].ra, int_ra, tolerance2 + 5))
               if( is_between( day_data[0][i].dec, day_data[1][i].dec, int_dec, tolerance2 + 5))
                  {
                  ELEMENTS class_elem;
                  double ra1, dec1, mag;
                  double d_ra, d_dec;
                  double earth_obj_dist, dist;
                  int sof_rval = -999;

                  n_checked++;
                  fseek( orbits_file, (i + 1) * record_length, SEEK_SET);
                  if( fgets( tbuff, sizeof( tbuff), orbits_file))
                     sof_rval = extract_sof_data( &class_elem, tbuff, sof_header);
                  if( sof_rval)
                     {
                     fprintf( stderr, "Couldn't read .sof elements: ast %d, rval %d\n",
                                    i, sof_rval);
                     if( sof_rval != -999)
                        fprintf( stderr, "%s", tbuff);
                     exit( -1);
                     }
                  class_elem.is_asteroid = 1;
                  if( tbuff[1] == '/' || strchr( "APXCD", tbuff[3]))
                     class_elem.is_asteroid = 0;   /* it's a comet */
                  earth_obj_dist = compute_asteroid_loc( earth_loc, &class_elem, jd,
                           &ra1, &dec1);
                  mag = calc_obs_magnitude( &class_elem, obj_sun_dist,
                              earth_obj_dist, earth_sun_dist);
                  d_ra = centralize_angle( ra1 - ra) * cos_dec;
                  d_dec = dec1 - dec;
                  dist = sqrt( d_ra * d_ra + d_dec * d_dec);
                  dist *= radians_to_arcsec;
                  if( mag < mag_limit && dist < tolerance_in_arcsec)
                     {
                     double computed_ra_motion, computed_dec_motion;
                     double dt_in_hours = (jd2 - jd) * 24.;

                         /* Compute asteroid posn at second time for motion: */
                     compute_asteroid_loc( earth_loc2, &class_elem, jd2,
                              &computed_ra_motion, &computed_dec_motion);
                     computed_ra_motion =
                           centralize_angle( computed_ra_motion - ra1) * cos_dec;
                     computed_dec_motion -= dec1;
                                 /* cvt motions from radians/day to "/hour: */
                     computed_ra_motion *=  radians_to_arcsec / dt_in_hours;
                     computed_dec_motion *= radians_to_arcsec / dt_in_hours;
                     if( (fabs( computed_dec_motion - dec_motion) < motion_tolerance &&
                           fabs( computed_ra_motion - ra_motion) < motion_tolerance)
                                    || singleton_observation)
                        {
                        char mpcorb_info[240];
                        double ra2, dec2, lov_len, dist_from_lov;
                        int j;

                             /* Compute asteroid posn .1 days later, but same */
                             /* earth loc, for LOV computation:               */
                        compute_asteroid_loc( earth_loc, &class_elem, jd + .1,
                              &ra2, &dec2);
                        ra2 = centralize_angle( ra2 - ra) * cos_dec;
                        dec2 -= dec;
                        ra2 -= d_ra;       /* (ra2, dec2) is now a vector pointing */
                        dec2 -= d_dec;     /* along the LOV                        */
                        lov_len = sqrt( ra2 * ra2 + dec2 * dec2);
                        dist_from_lov = (d_ra * dec2 - ra2 * d_dec) / lov_len;

                        if( is_list_file)
                           {
                           if( ra1 < 0.)
                              ra1 += PI + PI;
                           snprintf( tbuff + 26, sizeof( tbuff) - 26,
                                 "%010.6f %+010.6f  %5.2f  %8.4f %8.4f",
                                          ra1 * 180. / PI, dec1 * 180. / PI, mag,
                                          computed_ra_motion /60., computed_dec_motion / 60.);
                           }
                        else
                           snprintf( tbuff + 26, sizeof( tbuff) - 26,
                              "%6.0f %6.0f  %6.0f  %4.1f %5.0f%5.0f",
                              -d_ra * radians_to_arcsec,
                              -d_dec * radians_to_arcsec, dist,
                              mag, computed_ra_motion, computed_dec_motion);
                        if( !class_elem.abs_mag)
                           memset( tbuff + 49, '-', 5);
                        memset( tbuff + 12, ' ', 14);
//                      snprintf( tbuff + strlen( tbuff), sizeof( tbuff) - strlen( tbuff),
//                                            "  %.4f", earth_obj_dist);
                        if( show_lov)
                           snprintf( tbuff + strlen( tbuff),
                                           sizeof( tbuff) - strlen( tbuff),
                                           "  %6.0f",
                                           dist_from_lov * radians_to_arcsec);
                        if( mpcorb_extracts &&
                                   (!get_mpcorb_dot_dat_line( "mpcorb.dat", i, mpcorb_info)
                                 || !get_mpcorb_dot_dat_line( "MPCORB.DAT", i, mpcorb_info)))
                           {
                           const char *tptr = mpcorb_extracts;

                           while( *tptr)
                              {
                              int start, count;
                              char *endptr = tbuff + strlen( tbuff);

                              if( sscanf( tptr, "%d,%d", &start, &count) != 2)
                                 {
                                 fprintf( stderr, "Error parsing mpcorb extracts at '%s'\n", tptr);
                                 exit( -1);
                                 }
                              *endptr++ = ' ';
                              memcpy( endptr, mpcorb_info + start - 1, count);
                              endptr[count] = '\0';
                              while( *tptr > ' ' && *tptr != ';')
                                 tptr++;
                              while( *tptr == ' ' || *tptr == ';')
                                 tptr++;
                              }
                           }
                        if( is_list_file)
                           j = n_results;
                        else
                           for( j = 0; j < n_results
                                     && atof( results[j] + 39) < dist; j++)
                              ;
                        if( n_results > results_array_size - 2)
                           {
                           results_array_size <<= 1;
                           results = (char **)realloc( results,
                                       results_array_size * sizeof( char *));
                           }
                        memmove( results + j + 1, results + j,
                                           (n_results - j) * sizeof( char *));
                        results[j] = (char *)malloc( strlen( tbuff) + 1);
                        strcpy( results[j], tbuff);
                        n_results++;
                        }
                     }
                  }
            }
         for( i = 0; i < n_results; i++)
            {
            if( i < max_results)
               {
               printf( "%s\n", results[i]);
               n_lines_printed++;
               }
            free( results[i]);
            }
         if( verbose)
            printf( "%d objects had to be checked\n", n_checked);
         }
   for( i= 0; i < n_ilines; i++)
      free( ilines[i]);
   free( ilines);
   free( results);
   if( day_data[0])
      free( day_data[0]);
   if( day_data[1])
      free( day_data[1]);
   fclose( orbits_file);
   if( show_header)
      printf( "The apparent motion and arc length for each object are shown,  followed\n"
           "by a list of possible matches,  in order of increasing distance.  For\n"
           "each match,  the separation is shown,  both in RA and dec,  and then\n"
           "the 'total' separation,  all in arcseconds.  Next,  the magnitude and\n"
           "apparent motion of the possible match are shown.  All motions are in\n"
           "arcseconds per hour.\n");
   if( !mpc_station_file)
      printf( "ObsCodes.html not found; parallax wasn't included!\n");
   else
      fclose( mpc_station_file);
   if( show_header)
      printf( "\nRun time: %.1f seconds\n",
                  (double)clock( ) / (double)CLOCKS_PER_SEC);
                     /* If the output was quite long,  re-display */
                     /* the explanation of the output :           */
   if( n_lines_printed > 40 && show_header)
      show_astcheck_info( );
   return( 0);
}
