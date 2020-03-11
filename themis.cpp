/* Code to convert THEMIS ephemerides such as those at

http://soleil.ssl.berkeley.edu/ground_systems/products/THEMIS/THEMIS_D/2019_041/THEMIS_D.2019_041.OEM_CARA

   into the sort of state vectors that Find_Orb would emit,  in preparation
for making TLEs using eph2tle (see my Find_Orb repo on GitHub).

   The conversion is fairly straightforward.  The input files give an 'object
ID' equal to the NORAD catalogue number,  from which we can get the international
(YYYY-NNNA) designation.  We get a start and stop time in Gregorian calendar form.
Positions are given at one-minute intervals;   we take every 'out_freq' line and
output it.

   The resulting TLEs tend to cover about 18 days;  'them_cat' can then be
used to concatenate Themis TLE sets.            */

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
   fprintf( stderr, "'themis' requires the name of a THEMIS-A, D,  or E ephemeris file as\n"
                    "a command-line argument.  Output defaults to stdout.\n");
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

         /* JD 2454101.5 = 2007 Jan 1.  The THEMIS satellites were launched */
         /* early that year.  We shouldn't see ephems before this date.    */
#define MIN_JD 2454101.5

int main( const int argc, const char **argv)
{
   double jd_start = 0., jd_end = 0.;
   char buff[200];
   FILE *ifile, *ofile = stdout;
   int output_line = 0, norad_id = 0, out_freq = 144, i;

   if( argc < 2)
      error_exit( NULL);
   ifile = fopen( argv[1], "rb");
   if( !ifile)
      {
      fprintf( stderr, "Couldn't open '%s' ", argv[1]);
      perror( NULL);
      error_exit( NULL);
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
               if( !out_freq || (1440 % out_freq))
                  error_exit( "The output frequency must be a factor of 1440.\n");
               break;
            case 'o':
               ofile = fopen( arg, "wb");
               if( !ofile)
                  {
                  fprintf( stderr,  "Couldn't open '%s ", arg);
                  perror( NULL);
                  error_exit( NULL);
                  }
               break;
            default:
               fprintf( stderr,  "Unrecognized parameter '%s ", argv[i]);
               error_exit( NULL);
            }
         }
   while( fgets( buff, sizeof( buff), ifile))
      {
      if( !memcmp( buff, "OBJECT_ID", 9))
         norad_id = atoi( buff + 23);
      if( !memcmp( buff, "START_TIM", 9))
         jd_start = get_time_from_string( 0., buff + 23, FULL_CTIME_YMD, NULL);
      if( !memcmp( buff, "STOP_TIME", 9))
         jd_end   = get_time_from_string( 0., buff + 23, FULL_CTIME_YMD, NULL);
      if( buff[4] == '-' && buff[7] == '-' && buff[10] == 'T' && buff[13] == ':')
         {
         if( !output_line)
            {
            int intl_letter = 0;

            assert( jd_start > MIN_JD);
            assert( jd_end > MIN_JD);
            assert( norad_id);
            if( norad_id >= 30580 && norad_id <= 30582)
               intl_letter = 'A' + norad_id - 30580;
            if( norad_id == 30797 || norad_id == 30798)
               intl_letter = 'D' + norad_id - 30797;
            assert( intl_letter);
                     /* '-1,149597870.7,86400' = coords are mean
                        equatorial of date, in km and km/s. */

            fprintf( ofile, "%f %f %d -1,149597870.7,86400 (500) Geocentric: ",
                           jd_start, (double)out_freq / 1440.,
                           (int)( jd_end - jd_start) * 1440 / out_freq);
            fprintf( ofile, "NORAD %d = 2007-004%c = THEMIS %c\n",
                           norad_id, intl_letter, intl_letter);
            }
         if( !(output_line % out_freq))
            fprintf( ofile, "%f %s", jd_start + (double)output_line / 1440., buff + 24);
         output_line++;
         }
      }
   fclose( ifile);
   assert( jd_start > MIN_JD);
   assert( jd_end > MIN_JD);
   fprintf( ofile, "\nCreated from 'themis.cpp' (q.v.)\n");
   full_ctime( buff, jd_start, CALENDAR_JULIAN_GREGORIAN);
   fprintf( ofile, "Ephemeris start: %s\n", buff);

   full_ctime( buff, jd_end, CALENDAR_JULIAN_GREGORIAN);
   fprintf( ofile, "Ephemeris end:   %s\n", buff);
   if( ofile != stdout)
      fclose( ofile);
   return( 0);
}
