/* Code to add satellite offset data to 80-column MPC formatted
astrometry,  using coordinates downloaded from JPL Horizons.  Can be
compiled both as a standalone utility and as the code behind an
on-line one;  see https://www.projectpluto.com/add_off.htm .

   Note that SOHO data,  at least,  is also available at

https://sohowww.nascom.nasa.gov/data/ancillary/orbit/

   Haven't checked that source out carefully,  since Horizons
has worked well thus far,  but the above URL could come in handy. */

#include <errno.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "mpc_func.h"

typedef struct
{
   double jd, xyz[3];
   char mpc_code[4];
} offset_t;

const double tolerance = 1e-5;
int verbose = 0;
int n_positions_set = 0, n_positions_failed = 0;

/* If the observation is from a spacecraft,  return the JDE of the
observation.  (Horizons expects times for vector ephems in JDE,  not
UTC JDs.) */

static double get_sat_obs_jd( const char *buff)
{
   double jd = 0.;

   if( strlen( buff) > 80 && (buff[14] == 'S' || buff[14] == 's'))
      jd = extract_date_from_mpc_report( buff, NULL);
   if( jd < 2450000.)
      jd = 0.;
   else
      jd += td_minus_utc( jd) / seconds_per_day;
   return( jd);
}

/* The following conversion table is going to need occasional fixes.
Cas = Cassini,  SoO = Solar Orbiter,  and PSP = Parker Solar Probe
are _not_ official MPC codes.  */

typedef struct
{
   char mpc_code[4];
   int jpl_desig;
} xref_t;

static int get_horizons_idx( const char *mpc_code)
{
   static const xref_t xrefs[] = {
          {"245",     -79 },   /* Spitzer            */
          {"249",     -21 },   /* SOHO               */
          {"250",     -48 },   /* Hubble             */
          {"258", -139479 },   /* Gaia               */
          {"Cas",     -82 },   /* Cassini            */
          {"C49",    -234 },   /* STEREO-A           */
          {"C50",    -235 },   /* STEREO-B           */
          {"C51",    -163 },   /* WISE               */
          {"C52", -128485 },   /* Swift              */
          {"C53", -139089 },   /* NEOSSAT            */
          {"C54",     -98 },   /* New Horizons       */
          {"C55",    -227 },   /* Kepler             */
          {"C56", -141043 },   /* LISA Pathfinder    */
          {"C57",     -95 },   /* TESS               */
          {"PSP",     -96 },   /* Parker Solar Probe */
          {"JWT",    -170 },   /* James Webb (Space) Telescope */
          {"SoO",    -144 }};  /* Solar Orbiter      */
   size_t i;
   const size_t n_xrefs = sizeof( xrefs) / sizeof( xrefs[0]);

   for( i = 0; i < n_xrefs; i++)
      if( !memcmp( mpc_code, xrefs[i].mpc_code, 3))
         return( xrefs[i].jpl_desig);
   return( 0);
}

/* The following modifies an 'S' (satellite RA/dec line) into an 's' line
(satellite offset from the center of the earth).  If the offset is large,
it is stored in AU;  otherwise,  km.  MPC appears to be somewhat more
forgiving of where the decimal place goes than their documentation states.
We basically pick it to provide maximum precision.  */

static int set_mpc_style_offsets( char *buff, const double *xyz)
{
   double maxval = 0;
   int i;

   for( i = 0; i < 3; i++)
      if( maxval < fabs( xyz[i]))
         maxval = fabs( xyz[i]);
   memset( buff + 33, ' ', 39);
   buff[33] = ' ';
   if( maxval > 9999999.0)     /* show in AU */
      {
      buff[32] = '2';         /* mark as units = AU */
      for( i = 0; i < 3; i++)
         sprintf( buff + 35 + i * 12,
                (maxval > 9.9 * AU_IN_KM ? "%11.8f" : "%11.9f"),
                     fabs( xyz[i] / AU_IN_KM));
      }
   else
      {
      const char *format;

      buff[32] = '1';         /* mark as units = km */
      if( maxval > 999999.0)   /* Gaia,  perhaps SOHO eventually  */
         format = "%11.3f";
      else if( maxval > 99999.0)   /* TESS,  e.g.  */
         format = "%11.4f";
      else
         format = "%11.5f";
      for( i = 0; i < 3; i++)
         sprintf( buff + 35 + i * 12, format, fabs( xyz[i]));
      }
   for( i = 0; i < 3; i++)
      buff[34 + i * 12] = (xyz[i] > 0. ? '+' : '-');
   buff[14] = 's';
   buff[70] = ' ';
   return( 0);
}

const char *cmd_start = "wget -O /tmp/locs %s"
    "\"https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='%d'"
    "&REF_PLANE='FRAME'"
    "&OBJ_DATA='NO'&TABLE_TYPE='V'&TLIST=";

const char *cmd_end = "&VEC_TABLE='2'&VEC_LABELS='N'\"";

/* To set the time offsets,  we send a query to JPL Horizons that will
look something like the following (split here over four lines for ease
of explanation) :

https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='-163'
&REF_PLANE='FRAME'&OBJ_DATA='NO'&TABLE_TYPE='V'&TLIST=
'2458843.421181','2458843.486631','2458843.551951','2458843.616891',
&VEC_TABLE='2'&VEC_LABELS='N'

This requests positions on the four JDEs given on the third line for
object -163 (which is Horizons' index for (C51) WISE.)  REF_PLANE='FRAME'
specifies J2000 equatorial coordinates.  TABLE_TYPE='V' specifies
vectors.  VEC_TABLE='2' specifies positions and velocities... we actually
could skip the latter,  but I have some possible use cases for them
(adjusting for timing errors,  for example).

   Each time adds 17 bytes to our URL.  I've yet to determine what
JPL's limit is on URL length.  A buffer size of 1700 bytes,  after
allowing for the header and trailer data,  accommodates 89 times.
So if we encounter an unset offset,  we look for up to 88 other
instances where that particular obscode was used,  form a query to
ask for all of them,  and then set up to 89 offsets at a go.  */

static int set_offsets( offset_t *offsets, const int n_offsets)
{
   char buff[1700];     /* supports 89 times at a go */
   int i;
   const int horizons_idx = get_horizons_idx( offsets->mpc_code);

   if( !horizons_idx)
      {
      printf( "ERROR! MPC code '%s' wasn't found.\n", offsets->mpc_code);
      printf( "Either it's not an MPC code,  or it's not one of the spacecraft\n");
      printf( "that this software knows about.  Check the 'add_loc.c' source\n");
      printf( "code,  and/or contact the author.\n");
      for( i = 0; i < n_offsets; i++)
         if( !strcmp( offsets[i].mpc_code, offsets[0].mpc_code))
            offsets[i].xyz[0] = -0.1;     /* mark as "don't try again" */
      return( 0);
      }
   snprintf( buff, sizeof( buff), cmd_start,
            (verbose ? "" : "-q "), horizons_idx);
   for( i = 0; i < n_offsets; i++)
      if( !strcmp( offsets[i].mpc_code, offsets[0].mpc_code))
         {
         if( i)
            strcat( buff, ",");
         sprintf( buff + strlen( buff), "'%.6f'", offsets[i].jd);
         if( strlen( buff) + 60 > sizeof( buff))
            break;
         }
   strcat( buff, cmd_end);
   if( verbose)
      printf( "%s\n", buff);
   if( !system( buff))
      {
      FILE *ifile = fopen( "/tmp/locs", "rb");

      assert( ifile);
      while( fgets( buff, sizeof( buff), ifile))
         if( strstr( buff, " = A.D. ") && strstr( buff, " TDB"))
            {
            const double jd = atof( buff);
            int j;

            if( verbose)
               printf( "Found locations\n%s", buff);
            if( fgets( buff, sizeof( buff), ifile))
               for( j = 0; j < i; j++)
                  if( !strcmp( offsets[j].mpc_code, offsets[0].mpc_code)
                           && fabs( offsets[j].jd - jd) < tolerance)
                     {
                     const int n_found = sscanf( buff, "%lf %lf %lf",
                           &offsets[j].xyz[0], &offsets[j].xyz[1], &offsets[j].xyz[2]);

                     if( verbose)
                        printf( "Set for offset %d\n", j);
                     assert( n_found == 3);
                     n_positions_set++;
                     }
            }
         else if( verbose > 1 || !memcmp( buff, "No ephemeris", 12))
            printf( "%s", buff);
      fclose( ifile);
      }
   else
      {
      printf( "Error with system() : '%s'\n", strerror( errno));
      printf( "'%s'\n", buff);
      }
         /* If some or all obs weren't set,  zero their MPC codes.  That */
         /* will keep us from making repeated failed requests for them.  */
   for( i--; i > 0; i--)
      if( !strcmp( offsets[i].mpc_code, offsets[0].mpc_code))
         if( !offsets[i].xyz[0])
            {
            n_positions_failed++;
            offsets[i].mpc_code[0] = '\0';
            }
   return( 0);
}

/* The following reads the input file and looks for 80-column obs from
spacecraft.  On a second pass,  it removes any existing 's' (spacecraft
position) records and replaces them with 's' records created from
the Horizons ephems.   */

int process_file( const char *filename)
{
   FILE *ifile = fopen( filename, "rb");
   char buff[300];
   double jd;
   offset_t *offsets = NULL;
   int i, n_offsets = 0;

   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      if( buff[14] == 'S' && (jd = get_sat_obs_jd( buff)) != 0.)
         {
         if( verbose)
            printf( "Sat obs: %.5f\n%s", jd, buff);
         n_offsets++;
         offsets = (offset_t *)realloc( offsets,
                              n_offsets * sizeof( offset_t));
         memset( offsets + n_offsets - 1, 0, sizeof( offset_t));
         offsets[n_offsets - 1].jd = jd;
         memcpy( offsets[n_offsets - 1].mpc_code, buff + 77, 3);
         }
   for( i = 0; i < n_offsets; i++)
      {
      if( verbose)
         printf( "%d: JD %.5f; code '%s'\n", i, offsets[i].jd, offsets[i].mpc_code);
      if( !offsets[i].xyz[0] && offsets[i].mpc_code[0])
         set_offsets( offsets + i, n_offsets - i);
      }
   fseek( ifile, 0, SEEK_SET);
   while( fgets( buff, sizeof( buff), ifile))
      if( (jd = get_sat_obs_jd( buff)) <= 0.)
         printf( "%s", buff);
      else if( buff[14] != 's')
         {
         printf( "%s", buff);
         for( i = 0; i < n_offsets; i++)
            if( !memcmp( buff + 77, offsets[i].mpc_code, 3)
                       && fabs( jd - offsets[i].jd) < tolerance)
               set_mpc_style_offsets( buff, offsets[i].xyz);
         printf( "%s", buff);
         }
   fclose( ifile);
   free( offsets);
   printf( "COM %d positions set by add_off; %d failed in %.2f seconds\n",
         n_positions_set, n_positions_failed,
         (double)clock( ) / (double)CLOCKS_PER_SEC);
   printf( "COM add_off ver " __DATE__ " " __TIME__ "\n");
   return( 0);
}

#ifdef ON_LINE_VERSION
int dummy_main( const int argc, const char **argv)
#else
int main( const int argc, const char **argv)
#endif
{
   int i;

   if( argc < 2)
      {
      fprintf( stderr, "'add_loc' takes the name of an input file of astrometry\n"
                       "as a command-line argument.\n");
      return( -1);
      }
   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'v':
               verbose = 1;
               if( atoi( argv[i] + 2))
                  verbose = atoi( argv[i] + 2);
               printf( "Verbose = %d\n", verbose);
               break;
            default:
               printf( "Option '%s' unrecognized\n", argv[i]);
               break;
            }
   return( process_file( argv[1]));
}

#ifdef ON_LINE_VERSION
#include "cgi_func.h"

int main( void)
{
   const char *argv[20];
   const size_t max_buff_size = 40000;
   char *buff = (char *)malloc( max_buff_size);
   char field[30];
   FILE *lock_file = fopen( "lock.txt", "w");
   extern char **environ;
   int cgi_status;
   const char *temp_filename = "/tmp/add_off.txt";
   size_t bytes_written = 0;

   avoid_runaway_process( 15);
   printf( "Content-type: text/html\n\n");
   printf( "<html> <body> <pre>\n");
   if( !lock_file)
      {
      printf( "<p> Server is busy.  Try again in a minute or two. </p>");
      printf( "<p> Your TLEs are very important to us! </p>");
      return( 0);
      }
   setbuf( lock_file, NULL);
   fprintf( lock_file, "'add_off' : We're in\n");
   for( size_t i = 0; environ[i]; i++)
      fprintf( lock_file, "%s\n", environ[i]);
   cgi_status = initialize_cgi_reading( );
   fprintf( lock_file, "CGI status %d\n", cgi_status);
   if( cgi_status <= 0)
      {
      printf( "<p> <b> CGI data reading failed : error %d </b>", cgi_status);
      printf( "This isn't supposed to happen.</p>\n");
      return( 0);
      }
   while( !get_cgi_data( field, buff, NULL, max_buff_size))
      {
      if( !strcmp( field, "TextArea") || !strcmp( field, "upfile"))
         {
         if( strlen( buff) > 70)
            {
            char *tptr = strstr( buff, "COM verbo");
            int i;
            FILE *ofile = fopen( temp_filename,
                               (bytes_written ? "ab" : "wb"));

            assert( ofile);
            for( i = 1; buff[i + 1]; i++)     /* cvt Mac-style CR endings */
               if( buff[i] == 13 && buff[i + 1] != 10)   /* to LF endings */
                  buff[i] = 10;
            bytes_written += fwrite( buff, 1, strlen( buff), ofile);
            fprintf( ofile, "\n");
            fclose( ofile);
            if( tptr)
               verbose = 1;
            }
         }
      }
   fprintf( lock_file, "%d bytes written\n", (int)bytes_written);
   free( buff);
   argv[0] = "add_off";
   argv[1] = temp_filename;
   argv[2] = NULL;
   dummy_main( 2, argv);
   fprintf( lock_file, "Done\n");
   printf( "</pre> </body> </html>");
   fclose( lock_file);
   return( 0);
}
#endif
