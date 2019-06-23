/* Code to read in files of the form

https://lasp.colorado.edu/mms/sdc/public/data/ancillary/mms1/predeph/MMS1_PREDEPH_2019155_2019171.V00

and spit out ephemeris files of the sort readable by 'eph2tle' (q.v.,
in the Find_Orb repository).  The files give times in days since
1958 Jan 1 = JD 2436204.5 TAI.  TAI differs from TDT (the ephemeris
time scale we actually want) by a solidly fixed 32.184 seconds.

   This is somewhat analogous to 'themis.cpp' (q.v.),  written to
produce accurate TLEs for an object where Space-Track is not entirely
effectual.  However,  while the TLEs generated in this manner are
really good and would enable us to determine which of the four MMS
satellites we're looking at,  the Space-Track TLEs appear to be good
enough for the general problem of at least being able to say : "you
got one of the four MMS sats;  how much do we really care which one
you got?"

   If we decide that _really_ matters,  I'll set this up to generate
TLEs for all four MMS sats,  as 'themis' currently does for
THEMIS-A, D,  and E.  But that's not a priority at present. */

#define TDT_MINUS_TAI      32.184
#define JAN_1_1958         (2436204.5 + TDT_MINUS_TAI / 86400.)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "watdefs.h"
#include "date.h"

static void error_exit( const char *message)
{
   if( message)
   fprintf( stderr, "%s\n", message);
   fprintf( stderr, "'mms' requires the name of an MMS ephemeris file as\n"
                    "a command-line argument.  Output defaults to 'mms.eph'.\n");
   exit( -1);
}

/* I really should use getopt() or a portable variant.  However,  this has
been sufficiently effective thus far... */

static const char *get_arg( const int argc, const char **argv, const int idx)
{
   if( argv[idx][2] || idx == argc - 1)
      return( argv[idx] + 2);
   else
      return( argv[idx + 1]);
}

static double mms_ephem_jd( const char *buff)
{
   if( strlen( buff) > 150 && buff[4] == '-' && buff[8] == '/'
               && buff[11] == ':')
      return( atof( buff + 23) + JAN_1_1958);
   else
      return( 0.);
}

int main( const int argc, const char **argv)
{
   double jd_start = 0.;
   char buff[300];
   FILE *ifile, *ofile = stdout;
   int n_lines = 0, out_freq = 72, i;
   int mms_number = 0;

   if( argc < 2)
      error_exit( NULL);
   ifile = fopen( argv[1], "rb");
   if( !ifile)
      {
      fprintf( stderr, "Couldn't open '%s' ", argv[1]);
      perror( NULL);
      error_exit( NULL);
      }
   if( argc > 2 && argv[2][0] != '-')
      {
      ofile = fopen( argv[2], "wb");
      if( !ifile)
         {
         fprintf( stderr, "Couldn't open '%s' ", argv[2]);
         perror( NULL);
         error_exit( NULL);
         }
      }
   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         {
         const char *arg = get_arg( argc, argv, i);

         assert( arg);
         switch( argv[i][1])
            {
            case 'f':
               out_freq = atoi( arg);
               if( !out_freq || (720 % out_freq))
                  error_exit( "The output frequency must be a factor of 720.\n");
               break;
            default:
               fprintf( stderr,  "Unrecognized parameter '%s ", argv[i]);
               error_exit( NULL);
            }
         }
   while( fgets( buff, sizeof( buff), ifile))
      {
      const double jd = mms_ephem_jd( buff);

      if( jd)
         {
         n_lines++;
         if( !jd_start)
            jd_start = jd;
         }
      if( !memcmp( buff, "Spacecraft = MMS", 16))
         mms_number = buff[16];
      }
   assert( mms_number >= '1' && mms_number <= '4');
   fseek( ifile, 0L, SEEK_SET);
                     /* '0,149597870.7,86400' = coords are J2000 equatorial,
                        in km and km/s. */
   fprintf( ofile, "%.6f %.6f %d 0,149597870.7,86400 (500) Geocentric: MMS%c = 2015-011%c\n",
            jd_start, (double)out_freq / 720., n_lines / out_freq,
            mms_number, mms_number + 'A' - '1');

   n_lines = 0;
   while( fgets( buff, sizeof( buff), ifile))
      {
      const double jd = mms_ephem_jd( buff);

      if( jd)
         {
         buff[184] = '\0';
         if( n_lines % out_freq == 0)
            fprintf( ofile, "%.6f %s\n", jd, buff + 44);
         n_lines++;
         }
      }
   fclose( ifile);
   fprintf( ofile, "\nCreated from 'mms.cpp' (q.v.)\n");
   full_ctime( buff, jd_start, CALENDAR_JULIAN_GREGORIAN);
   fprintf( ofile, "Ephemeris start: %s\n", buff);

   full_ctime( buff, jd_start + (n_lines / out_freq) * out_freq / 720., CALENDAR_JULIAN_GREGORIAN);
   fprintf( ofile, "Ephemeris end:   %s\n", buff);
   if( ofile != stdout)
      fclose( ofile);
   return( 0);
}
