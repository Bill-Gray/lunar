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
#include "watdefs.h"
#include "date.h"
#include "comets.h"
#include "afuncs.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define LOG_10 2.302585
const double radians_to_arcsec = 180. * 3600. / PI;

int get_earth_loc( const double t_millennia, double *results);
long extract_mpcorb_dat( ELEMENTS *elem, const char *buff);   /* mpcorb.c */
long extract_astorb_dat( ELEMENTS *elem, const char *buff);   /* mpcorb.c */

static int get_mpc_data( const char *buff, double *jd, double *ra, double *dec)
{
   int i1, i2, year, month;
   double tval, day;

   if( strlen( buff) < 80 || strlen( buff) > 83)
      return( -1);
   if( sscanf( buff + 32, "%d %d %lf", &i1, &i2, &tval) != 3)
      return( -2);
   *ra = ((double)i1 + (double)i2 / 60. + tval / 3600.) * (PI / 12.);

   if( sscanf( buff + 45, "%d %d %lf", &i1, &i2, &tval) != 3)
      return( -3);
   *dec = ((double)i1 + (double)i2 / 60. + tval / 3600.) * (PI / 180.);
   if( buff[44] == '-')
      *dec = -*dec;
   else if( buff[44] != '+')
      return( -4);

               /* Read in the day/month/year from the record... */
   if( sscanf( buff + 15, "%d %d %lf", &year, &month, &day) != 3)
      return( -5);
               /* ...and convert to a JD: */
   *jd = day + (double)dmy_to_day( 0, month, year, 0) - .5;
   return( 0);
}

static inline double law_of_cosines( const double a, const double b, const double c)
{
   return( .5 * (a * a + b * b - c * c) / (a * b));
}

static double calc_obs_magnitude( ELEMENTS *elem, const double obj_sun,
                      const double obj_earth, const double earth_sun)
{
   double magnitude;

   if( !elem->is_asteroid)
      magnitude = elem->slope_param * log( obj_sun);
   else
      {
      const double cos_phase_ang =
                  law_of_cosines( obj_sun, obj_earth, earth_sun);
      const double half_phase_ang = acose( cos_phase_ang) / 2.;
      const double log_tan_half_phase = log( tan( half_phase_ang));
      const double phi1 = exp( -3.33 * exp( log_tan_half_phase * 0.63));
      const double phi2 = exp( -1.87 * exp( log_tan_half_phase * 1.22));

      if( cos_phase_ang > 1. || cos_phase_ang < -1.)
         printf( "???? Triangle error: %f %f %f\n", obj_sun, obj_earth, earth_sun);
      magnitude = 5. * log( obj_sun)
                           -2.5 * log( (1. - elem->slope_param) * phi1
                                            + elem->slope_param * phi2);
      }
   magnitude += 5. * log( obj_earth);
   magnitude /= LOG_10;      /* cvt from natural logs to common (base 10) */
   magnitude += elem->abs_mag;
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

static FILE *orbits_file;
static int n_asteroids, astorb_epoch, record_length, verbose = 0;
static int n_numbered, using_mpcorb = 0;
static size_t mpcorb_header_length;
const char *data_path = NULL;

FILE *get_mpcorb_file( void)
{
   const char *filename = "mpcorb.dat";
   FILE *ifile = fopen( filename, "rb");
   char buff[250];

   if( !ifile && data_path)
      {
      strcpy( buff, data_path);
      strcat( buff, filename);
      ifile = fopen( buff, "rb");
      }

   if( ifile)
      {
      int lines_read = 0, i;

      if( !fgets( buff, sizeof( buff), ifile))
         {
         printf( "Unable to read mpcorb.dat\n");
         exit( -4);
         }
      for( i = 0; buff[i]; i++)
         if( buff[i] == '.')
            buff[i] = ' ';
      if( i > 70)       /* put date into YYYYMMDD form */
         {
         const double jd = get_time_from_string( 0., buff + 60, 0, NULL);

         full_ctime( buff, jd, FULL_CTIME_YMD | FULL_CTIME_NO_SPACES
                     | FULL_CTIME_DATE_ONLY | FULL_CTIME_MONTHS_AS_DIGITS
                     | FULL_CTIME_LEADING_ZEROES);
         astorb_epoch = atol( buff);
         }
      while( lines_read < 50 && !mpcorb_header_length &&
                                      fgets( buff, sizeof( buff), ifile))
         {
         lines_read++;
         if( *buff == '-')    /* we've read the entire header */
            mpcorb_header_length = ftell( ifile);
         }
      record_length = 203;
      fseek( ifile, 0L, SEEK_END);
      n_asteroids = (int)(ftell( ifile) - mpcorb_header_length) / record_length;
      fseek( ifile, 0L, SEEK_SET);
      if( verbose)
         printf( "MPCORB epoch: %d; %d objects; record size %d\n",
               astorb_epoch, n_asteroids, record_length);
      }
   return( ifile);
}

FILE *get_astorb_file( void)
{
   const char *filename = "astorb.dat";
   FILE *ifile = fopen( filename, "rb");
   char tbuff[300];

   if( !ifile && data_path)
      {
      strcpy( tbuff, data_path);
      strcat( tbuff, filename);
      ifile = fopen( tbuff, "rb");
      }

   if( ifile)
      {
      long filesize;

      if( !fgets( tbuff, sizeof( tbuff), ifile))
         {
         printf( "Unable to read astorb.dat\n");
         exit( -5);
         }
      record_length = (int)strlen( tbuff);
      fseek( ifile, 0L, SEEK_END);
      filesize = ftell( ifile);
      n_asteroids = filesize / record_length;
      fseek( ifile, 0L, SEEK_SET);
      astorb_epoch = atoi( tbuff + 207);
                  /* File size should be an even multiple of the */
                  /* record size... if not,  something's off:    */
      if( filesize != n_asteroids * record_length)
         {
         printf( "ASTORB file size appears to be wrong: %ld\n",
                  filesize);
         verbose++;
         }
      if( verbose)
         printf( "ASTORB epoch: %d; %d objects; record size %d\n",
               astorb_epoch, n_asteroids, record_length);
      }
   return( ifile);
}

static double obj_sun_dist;

static double compute_asteroid_loc( const double *earth_loc,
            ELEMENTS *elem, const double jd, double *ra, double *dec)
{
   double r1 = 0., dist, asteroid_loc[4];
   int i;

   do             /* light-time lag:  should converge _very_ fast       */
      {           /* it'll almost always require exactly two iterations */
      dist = r1;
      comet_posn( elem, jd - dist / AU_PER_DAY, asteroid_loc);
      for( i = 0; i < 3; i++)
         asteroid_loc[i] -= earth_loc[i];
      r1 = vector3_length( asteroid_loc);
      }
      while( fabs( dist - r1) > .001);
   ecliptic_to_equatorial( asteroid_loc);
   *ra = atan2( asteroid_loc[1], asteroid_loc[0]);
   *dec = asin( asteroid_loc[2] / r1);
   obj_sun_dist = asteroid_loc[3];
   return( r1);
}


/* Six records in 'astorb.dat' contained ASCII zeroes in columns where the */
/* taxonomic class ought to be.  So for 'astorb',  we zero everything out. */

/* MPCORB has an intentional design flaw:  a line feeds is inserted after  */
/* the numbered objects,  and again after the multi-opp objects.  Hence,   */
/* we read in two bytes and check for line feeds,  which would indicate    */
/* we're actually a bit ahead of the record.                               */

static int get_orbit_rec( char *buff)
{
   int rval;

   if( using_mpcorb)
      {
      int bytes_read;

      if( !fread( buff, 2, 1, orbits_file))
         return( 0);
      if( buff[1] == 10)         /* single-opp object */
         bytes_read = 0;
      else if( buff[0] == 10)    /* multi-opp object */
         {
         buff[0] = buff[1];
         bytes_read = 1;
         }
      else        /* numbered object */
         bytes_read = 2;
      rval = (int)fread( buff + bytes_read, record_length - bytes_read, 1, orbits_file);
      }
   else
      {
      int i;

      rval = (int)fread( buff, record_length, 1, orbits_file);
      if( rval)
         for( i = 0; i < record_length; i++)
            if( !buff[i])
               buff[i] = ' ';
      }
   return( rval);
}

int16_t *ephem_uncertainties = NULL;

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
   fseek( orbits_file, (long)mpcorb_header_length, SEEK_SET);
   tbuff[record_length] = '\0';
   for( i = 0; i < n_asteroids && get_orbit_rec( tbuff); i++)
      {
      ELEMENTS class_elem;
      const double curr_ephem_unc = atof( tbuff + 190);
      double ra, dec;

      if( using_mpcorb)
         {
         if( !extract_mpcorb_dat( &class_elem, tbuff))
            {
            printf( "'mpcorb.dat' data doesn't parse correctly:\n%s\n", tbuff);
            exit( -1);
            }
         }
      else
         if( !extract_astorb_dat( &class_elem, tbuff))
            {
            printf( "'astorb.dat' data doesn't parse correctly:\n%s\n", tbuff);
            exit( -1);
            }

      compute_asteroid_loc( earth_loc, &class_elem, jd, &ra, &dec);

      rval[i].ra = integerize_angle( ra);
      rval[i].dec = integerize_angle( dec);
      if( curr_ephem_unc < 32000.)
         ephem_uncertainties[i] = (int16_t)curr_ephem_unc;
      else        /* store large values in arcminutes */
         ephem_uncertainties[i] = (int16_t)(-curr_ephem_unc / 60.);
      if( verbose && counter <= i * 80 / n_asteroids)
         {
         printf( "%d", counter % 10);
         counter++;
         }
      if( using_mpcorb && tbuff[172] == ')')
         n_numbered = i + 1;
      if( !using_mpcorb && tbuff[5] != ' ')
         n_numbered = i + 1;
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

/* In loading up the 'day data' (basically a low-precision RA/dec for */
/* the geocentric position of each asteroid as of a certain day),  we */
/* attempt to open a file for that day of the form YYYYMMDD.chk.  If  */
/* the file is opened,  and the header indicates the correct version, */
/* number of asteroids,  and astorb epoch,  we just load up the data  */
/* from the file and return.                                          */
/*   If we don't find the file,  or the header doesn't match,  then   */
/* we call the above 'compute_day_data',  and write that data out to  */
/* a .chk file so we don't have to recompute it all the next time.    */
/*   The first time this function is called,  we also load current    */
/* ephemeris uncertainty data,  either from the 'curr_unc' file (if   */
/* the header matched) or from astorb.dat itself (if it didn't.)      */
/*    NOTE that if we're using mpcorb.dat instead of astorb.dat,  the */
/* extension is changed to .chl.  Thus,  one can maintain .chl files  */
/* for mpcorb,  and .chk ones for astorb,  without collisions.        */

#define HEADER_SIZE 4

static AST_DATA *get_cached_day_data( const int ijd)
{
   char filename[20];
   FILE *ifile, *ofile;
   int32_t header[HEADER_SIZE];
   const int32_t magic_version_number = 1314159266;
   AST_DATA *rval;
   int uncertainties_loaded = (ephem_uncertainties != NULL);

   if( !ephem_uncertainties)
      ephem_uncertainties = (int16_t *)malloc( n_asteroids * sizeof( int16_t));
   if( !ephem_uncertainties)
      return( NULL);


                  /* Create a filename in 'YYYYMMDD.chk' form: */
   full_ctime( filename, (double)ijd, FULL_CTIME_YMD | FULL_CTIME_NO_SPACES
                     | FULL_CTIME_DATE_ONLY | FULL_CTIME_MONTHS_AS_DIGITS
                     | FULL_CTIME_LEADING_ZEROES);
   strcat( filename, (using_mpcorb ? ".chl" : ".chk"));
   if( verbose > 2)
      printf( "Creating '%s'\n", filename);
   ifile = fopen( filename, "rb");
   if( ifile)
      {
      if( !fread( header, HEADER_SIZE, sizeof( int), ifile))
         {
         printf( "Error reading header data in '%s'\n", filename);
         exit( -2);
         }
      if( header[0] != magic_version_number
                   || header[1] != astorb_epoch || header[2] != n_asteroids)
         fclose( ifile);
      else   /* appears to be legitimate cached data */
         {
         rval = (AST_DATA *)malloc( n_asteroids * sizeof( AST_DATA));
         if( !rval)
            {
            printf( "Ran out of memory\n");
            exit( -4);
            }
         n_numbered = header[3];
         if( !fread( rval, n_asteroids, sizeof( AST_DATA), ifile))
            {
            printf( "Error in asteroid data in '%s'\n", filename);
            exit( -3);
            }
         fclose( ifile);
         if( !uncertainties_loaded)
            {
            ifile = fopen( "curr_unc", "rb");
            if( !fread( ephem_uncertainties, n_asteroids, sizeof( int16_t), ifile))
               {
               printf( "Error in asteroid uncertainty data\n");
               exit( -4);
               }
            fclose( ifile);
            }
         return( rval);
         }
      }
   rval = compute_day_data( ijd);
   header[0] = magic_version_number;
   header[1] = astorb_epoch;
   header[2] = n_asteroids;
   header[3] = n_numbered;
   ofile = fopen( filename, "wb");
   fwrite( header, HEADER_SIZE, sizeof( int), ofile);
   fwrite( rval, n_asteroids, sizeof( AST_DATA), ofile);
   fclose( ofile);

   ofile = fopen( "curr_unc", "wb");
   fwrite( ephem_uncertainties, n_asteroids, sizeof( int16_t), ofile);
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

int qsort_stricmp( const void *elem1, const void *elem2)
{
   const char **buff1 = (const char **)elem1;
   const char **buff2 = (const char **)elem2;

   return( strcmp( *buff1, *buff2));
}

   /* To compute observed object motion,  this code grabs the RA/dec */
   /* of the first observation;  then looks ahead for observations   */
   /* of the same object from the same station.  The difference      */
   /* between the pair gives motion.  Except we want to use a pair   */
   /* far enough apart to show real motion,  but not so far apart    */
   /* that sky curvature or acceleration are factors.                */
   /*    Return value is the JD of the observation used to determine */
   /* the motion.                                                    */

static double compute_motion( const char **lines, const int n_lines,
                           double *ra_motion, double *dec_motion)
{
   double ra, dec, jd, rval = 0.;
   int i;

   *ra_motion = *dec_motion = 0.;
   get_mpc_data( lines[0], &jd, &ra, &dec);
   for( i = 1; i < n_lines && !memcmp( lines[0], lines[i], 12); i++)
      if( !memcmp( lines[0] + 77, lines[i] + 77, 3))
         {
         double ra1, dec1, jd1;
         const double five_degrees = PI / 36.;

         get_mpc_data( lines[i], &jd1, &ra1, &dec1);
         if( fabs( jd1 - jd) < 10.)   /* within ten days... */
            if( fabs( ra1 - ra) < five_degrees)
               if( fabs( dec1 - dec) < five_degrees)
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

static void err_message( void)
{
   printf( "\nastcheck needs the name of a file containing MPC-formatted (80-column)\n");
   printf( "astrometric data.  It will then attempt to match the records to known\n");
   printf( "objects in the 'astorb.dat' file.\n\n");
   printf( "Command-line options are:\n\n");
   printf( "   -r(dist)   Set search distance to 'dist' arcseconds. Default is 3600.\n");
   printf( "   -z(tol)    Set motion match tolerance to 'tol' arcsec/hr. Default is 10.\n");
   printf( "   -m(mag)    Set limiting mag to 'mag'.  Default is 22.\n");
   printf( "   -M         Use MPCORB,  not ASTORB.\n");
   printf( "   -l         Show distance from line of variations. Experimental.\n");
}

#ifdef USE_CEU_DATA
static int16_t cvt_tolerance( int16_t arcseconds)
{
// double tval = (double)arcseconds;
//
// if( tval < 0.)        /* "huge" value,  in arcminutes: */
//    tval = -tval / 60.;
// return( (int16_t)( tval * 65536. / (360. * 3600.)));
   if( arcseconds >= 0)
      return( arcseconds / 3);
   else          /* more than 32000 arcsec */
      return( 32000 / 3);
}
#endif

static void show_astorb_info( void)
{
   printf( "ASTCHECK version %s %s\n", __DATE__, __TIME__);
   printf( "astorb.dat version %d %d %d,  with %d objects (%d numbered)\n",
            astorb_epoch / 10000, (astorb_epoch / 100) % 100,
            astorb_epoch % 100, n_asteroids, n_numbered);
}

#ifdef _MSC_VER
     /* Microsoft Visual C/C++ has no snprintf.  Yes,  you read that      */
     /* correctly.  MSVC has an _snprintf which doesn't add a '\0' at the */
     /* end if max_len bytes are written.  You can't pass a NULL string   */
     /* to determine the necessary buffer size.  The following, however,  */
     /* is a "good enough" replacement:  for a non-NULL string, the       */
     /* output will be "correct" (the '\0' will always be present and in  */
     /* the right place).  The only deviation from "proper" sprintf() is  */
     /* that the maximum return value is max_len;  you can never know how */
     /* many bytes "should" have been written.                            */
     /*   Somewhat ancient MSVCs don't even have vsnprintf;  you have to  */
     /* use vsprintf and work things out so you aren't overwriting the    */
     /* end of the buffer.                                                */
#include <stdarg.h>

int snprintf( char *string, const size_t max_len, const char *format, ...)
{
   va_list argptr;
   int rval;

   va_start( argptr, format);
#if _MSC_VER <= 1100
   rval = vsprintf( string, format, argptr);
#else
   rval = vsnprintf( string, max_len, format, argptr);
#endif
   string[max_len - 1] = '\0';
   va_end( argptr);
   return( rval);
}
#endif    /* #ifdef _MSC_VER */

#define MAX_RESULTS 500

int main( const int argc, const char **argv)
{
   double jd, ra, dec;
   FILE *ifile;
   char buff[90];
   char **ilines;
   int show_lov = 0, show_uncertainty = 0;
   int i, n_ilines = 0, n, unnumbered_only = 0;
   int n_lines_printed = 0;
   double tolerance_in_arcsec = 18000.;       /* = five degrees */
   double mag_limit = 22.;
   AST_DATA *day_data[2] = { NULL, NULL};
   long curr_loaded_day_data = 0;
   FILE *mpc_station_file = fopen( "ObsCodes.html", "rb");
   char curr_station[7];
   double rho_sin_phi = 0., rho_cos_phi = 0., longitude = 0.;
   double motion_tolerance = 10.;  /* require a match to within 10"/hr */

   if( !mpc_station_file)        /* perhaps stored with truncated extension? */
      mpc_station_file = fopen( "ObsCodes.htm", "rb");
   if( !mpc_station_file)
      printf( "ObsCodes.html not found; parallax won't be included!\n");
   curr_station[0] = '\0';

   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'r':
               tolerance_in_arcsec = atof( argv[i] + 2);
               break;
            case 'v':
               setvbuf( stdout, NULL, _IONBF, 0);
               verbose = 1 + atoi( argv[i] + 2);
               break;
            case 'm':
               mag_limit = atof( argv[i] + 2);
               break;
            case 'M':
               using_mpcorb = 1;
               break;
            case 'p':
               data_path = argv[i] + 2;
               break;
            case 'l':
               show_lov = 1;
               break;
            case 'z':
               motion_tolerance = atof( argv[i] + 2);
               break;
            case 'u':
               unnumbered_only = 1;
               break;
            case 'e':
               show_uncertainty = 1;
               break;
            default:
               printf( "%s: unrecognized command-line option\n", argv[i]);
               break;
            }
   if( argc < 2)
      {
      printf( "No input file name specified\n");
      err_message( );
      return( -1);
      }
   ifile = fopen( argv[1], "rb");
   if( !ifile)
      printf( "%s not opened\n", argv[1]);
   if( using_mpcorb)
      orbits_file = get_mpcorb_file( );
   else
      {
      orbits_file = get_astorb_file( );
      if( !orbits_file)     /* astorb failed;  let's try mpcorb.dat */
         {
         orbits_file = get_mpcorb_file( );
         using_mpcorb = 1;
         }
      }
   if( !orbits_file)
      {
      printf( "Couldn't find astorb.dat");
      err_message( );
      return( -2);
      }
   if( !ifile)
      {
      err_message( );
      return( -3);
      }

               /* Run through input file and count lines of astrometry: */
   while( fgets( buff, sizeof( buff), ifile))
      if( !get_mpc_data( buff, &jd, &ra, &dec))
         n_ilines++;

   if( !n_ilines)
      {
      printf( "No astrometry found in '%s'\n", argv[1]);
      err_message( );
      return( -1);
      }
               /* Allocate memory for astrometry lines, then read 'em: */
   ilines = (char **)malloc( n_ilines * sizeof( char *));
   n_ilines = 0;
   fseek( ifile, 0L, SEEK_SET);
   while( fgets( buff, sizeof( buff), ifile))
      if( !get_mpc_data( buff, &jd, &ra, &dec))
         {
         ilines[n_ilines] = (char *)malloc( strlen( buff) + 1);
         strcpy( ilines[n_ilines], buff);
         n_ilines++;
         }
   qsort( ilines, n_ilines, sizeof( char **), qsort_stricmp);
   for( n = 0; n < n_ilines; n++)
      if( (!n || memcmp( ilines[n], ilines[n - 1], 12)) &&
                       !get_mpc_data( ilines[n], &jd, &ra, &dec))
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
         char *results[MAX_RESULTS + 1];
         char tbuff[300];
         int n_results = 0;
         int n_checked = 0;

         jd += delta_t;
         if( mpc_station_file && memcmp( ilines[n] + 77, curr_station, 3))
            {
            int got_station_data = 0;

            strcpy( curr_station, ilines[n] + 77);
            curr_station[3] = '\0';
            fseek( mpc_station_file, 0L, SEEK_SET);
            rho_sin_phi = rho_cos_phi = longitude = 0.;
            while( !got_station_data &&
                           fgets( tbuff, sizeof( tbuff), mpc_station_file))
               got_station_data = !memcmp( tbuff, curr_station, 3);
            if( got_station_data)
               sscanf( tbuff + 3, "%lf%lf%lf",
                                     &longitude, &rho_cos_phi, &rho_sin_phi);
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
         if( !n)     /* on our very first object: */
            {
            show_astorb_info( );
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
         if( ra_motion || dec_motion)
            printf( "\n%s: %.0f\"/hr in RA, %.0f\"/hr in dec (%.2f hours)\n",
                        buff, ra_motion, dec_motion, (jd2 - jd) * 24.);
         else
            printf( "\n%s: only one observation\n", buff);
         n_lines_printed++;
         for( i = (unnumbered_only ? n_numbered : 0); i < n_asteroids; i++)
            {
#ifdef USE_CEU_DATA
            const int16_t tolerance2 = tolerance + cvt_tolerance( day_data[i].uncertainty);
#else
            const int16_t tolerance2 = tolerance;
#endif

            if( is_between( day_data[0][i].ra, day_data[1][i].ra, int_ra, tolerance2 + 5))
               if( is_between( day_data[0][i].dec, day_data[1][i].dec, int_dec, tolerance2 + 5))
                  {
                  ELEMENTS class_elem;
                  double ra1, dec1, mag;
                  double earth_obj_dist, dist;

                  n_checked++;
                  fseek( orbits_file, i * record_length + (long)mpcorb_header_length, SEEK_SET);
                  get_orbit_rec( tbuff);
                  if( using_mpcorb)
                     {
                     extract_mpcorb_dat( &class_elem, tbuff);
                     memcpy( tbuff, tbuff + 166, 26);
                     }
                  else
                     extract_astorb_dat( &class_elem, tbuff);
                  tbuff[26] = '\0';
                  earth_obj_dist = compute_asteroid_loc( earth_loc, &class_elem, jd,
                           &ra1, &dec1);
                  mag = calc_obs_magnitude( &class_elem, obj_sun_dist,
                              earth_obj_dist, earth_sun_dist);
                           /* Cvt ra1, dec1 to be relative to observation point: */
                  ra1 = centralize_angle( ra1 - ra) * cos_dec;
                  dec1 -= dec;
                  dist = sqrt( ra1 * ra1 + dec1 * dec1);
                  dist *= radians_to_arcsec;
                  if( mag < mag_limit && dist < tolerance_in_arcsec)
                     {
                     double computed_ra_motion, computed_dec_motion;
                     double dt_in_hours = (jd2 - jd) * 24.;

                         /* Compute asteroid posn at second time for motion: */
                     compute_asteroid_loc( earth_loc2, &class_elem, jd2,
                              &computed_ra_motion, &computed_dec_motion);
                     computed_ra_motion =
                           centralize_angle( computed_ra_motion - ra) * cos_dec;
                     computed_ra_motion -= ra1;
                     computed_dec_motion -= dec + dec1;
                                 /* cvt motions from radians/day to "/hour: */
                     computed_ra_motion *=  radians_to_arcsec / dt_in_hours;
                     computed_dec_motion *= radians_to_arcsec / dt_in_hours;
                     if( fabs( computed_dec_motion - dec_motion) < motion_tolerance &&
                           fabs( computed_ra_motion - ra_motion) < motion_tolerance)
                        {
                        double ra2, dec2, lov_len, dist_from_lov;
                        int j;

                             /* Compute asteroid posn .1 days later, but same */
                             /* earth loc, for LOV computation:               */
                        compute_asteroid_loc( earth_loc, &class_elem, jd + .1,
                              &ra2, &dec2);
                        ra2 = centralize_angle( ra2 - ra) * cos_dec;
                        dec2 -= dec;
                        ra2 -= ra1;       /* (ra2, dec2) is now a vector pointing */
                        dec2 -= dec1;     /* along the LOV                        */
                        lov_len = sqrt( ra2 * ra2 + dec2 * dec2);
                        dist_from_lov = (ra1 * dec2 - ra2 * dec1) / lov_len;

                        snprintf( tbuff + 26, sizeof( tbuff) - 26,
                              "%6.0f %6.0f  %6.0f  %4.1f %5.0f%5.0f",
                              -ra1 * radians_to_arcsec,
                              -dec1 * radians_to_arcsec, dist,
                              mag, computed_ra_motion, computed_dec_motion);
//                      snprintf( tbuff + strlen( tbuff), sizeof( tbuff) - strlen( tbuff),
//                                            "  %.4f", earth_obj_dist);
                        if( show_lov)
                           snprintf( tbuff + strlen( tbuff),
                                           sizeof( tbuff) - strlen( tbuff),
                                           "  %6.0f",
                                           dist_from_lov * radians_to_arcsec);
                        if( show_uncertainty)
                           {
                           int uncertainty = (int)ephem_uncertainties[i];

                           if( uncertainty < 0)
                              uncertainty *= -60;
                           snprintf( tbuff + strlen( tbuff),
                                       sizeof( tbuff) - strlen( tbuff),
                                       "  %d\"", uncertainty);
                           }
                        for( j = 0; j < n_results
                                     && atof( results[j] + 39) < dist; j++)
                           ;
                        memmove( results + j + 1, results + j,
                                           (n_results - j) * sizeof( char *));
                        results[j] = (char *)malloc( strlen( tbuff) + 1);
                        strcpy( results[j], tbuff);
                        if( n_results < MAX_RESULTS)
                           n_results++;
                        }
                     }
                  }
            }
         for( i = 0; i < n_results; i++)
            {
            char number[10];

            number[0] = ' ';
            number[8] = '\0';
            memcpy( number + 1, results[i], 7);
            if( number[6] != ' ')         /* yes,  it's a numbered object */
               {                          /* put it in parentheses        */
               int j;

               number[7] = ')';
               for( j = 6; number[j] != ' '; j--)
                  ;
               number[j] = '(';
               }
            printf( "%s %s\n", number, results[i] + 7);
            free( results[i]);
            }
         n_lines_printed += n_results;
         if( verbose)
            printf( "%d objects had to be checked\n", n_checked);
         }
   for( i= 0; i < n_ilines; i++)
      free( ilines[i]);
   free( ilines);
   if( day_data[0])
      free( day_data[0]);
   if( day_data[1])
      free( day_data[1]);
   if( !mpc_station_file)
      printf( "ObsCodes.html not found; parallax wasn't included!\n");
   else
      fclose( mpc_station_file);
   printf( "\nRun time: %.1f seconds\n",
                  (double)clock( ) / (double)CLOCKS_PER_SEC);
                     /* If the output was quite long,  re-display */
                     /* the 'astorb' information:                 */
   if( n_lines_printed > 40)
      show_astorb_info( );
   return( 0);
}
