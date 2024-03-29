/* Code to process Gaia asteroid observation data into ADES format.

   Gaia has provided asteroid astrometry in three releases,  each time
in slightly different CSV formats that are non-standard and a little
oddly organized.  The most recent (and generally preferred) data are
in the FPR (Focussed Product Release) :

https://cdn.gea.esac.esa.int/Gaia/gfpr/Solar_system/sso_observation/

   These do lack photometry,  but the preceding Gaia-DR3 files,

http://cdn.gea.esac.esa.int/Gaia/gdr3/Solar_system/sso_observation/

   provide that.  This program no longer processes the eldest Gaia-DR2
data from

http://cdn.gea.esac.esa.int/Gaia/gdr2/sso_observation/csv/

   To use this program,  download and un-GZip the necessary CSV file(s)
from the above sites.  Run this program as,  for example,

./csv2ades 4179 > output.ades
./csv2ades 4179,101955 > output.ades

   The first example would output all the observations for (4179)
Toutatis to 'output.ades'. The second example would output observations
for all objects from (4179) Toutatis to (101955) Bennu,  inclusive.
Note that the output files can get to be pretty big.  You may be
better off just running an object or a few at a time.

   There are 20 GZIPped files in the above directory for FPR,  each about
1.3 GBytes,  27 GBytes for all twenty.  Each file contains astrometry for
about 1/20 of the numbered minor planets observed by Gaia. Numbered
object (N) is placed in the file

SsoObservation_xx.csv

   where xx = (N + 11) mod 20.  Thus,  for example,  the data for
(118186) would be in the file SsoObservation_17.csv.  The reason
for the offset by 11 eludes me.  It does mean that if you're looking
for a specific object,  you can figure out which file will have data
for it and then download only that file.  (Or you can run the above
command and get an error message telling you which file wasn't found.) */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"

static char *get_csv( char *obuff, const char *ibuff,
            const char *header, const char *tag, const size_t max_len)
{
   size_t i;
   unsigned nval = 0;
   const char *found = strstr( header, tag);

   if( !found)
      fprintf( stderr, "Couldn't find '%s'\n", tag);
   assert( found);
   assert( found >= header);
   for( i = 0; i < (size_t)(found - header); i++)
      if( header[i] == ',')
         nval++;
   while( nval-- && ibuff)
      {
      ibuff = strchr( ibuff, ',');
      if( ibuff)
         ibuff++;
      }
   assert( ibuff);
   for( i = 0; i < max_len && ibuff[i] && ibuff[i] != ','; i++)
      obuff[i] = ibuff[i];
   obuff[i] = '\0';
   return( obuff);
}

static void csv_to_ades( const char *ibuff, const char *header,
                          const int low_num, const int high_num)
{
   char tbuff[90];
   int num;

   get_csv( tbuff, ibuff, header, "number_mp,", 9);
   num = atoi( tbuff);
   if( num >= low_num && num <= high_num)
      {
      const double jd_epoch = 2455197.5;    /* 2010 Jan 1.0 */
      double jd_obs;
      int i;

      printf( "   <optical>\n");
      printf( "    <permID>%s</permID>\n", tbuff);
      printf( "    <mode>TDI</mode>\n");
      printf( "    <stn>258</stn>\n");
      printf( "    <sys>ICRF_KM</sys>\n");
      printf( "    <ctr>399</ctr>\n");
      for( i = 0; i < 6; i++)
         {
         char tag[30];
         double value;
         const double LB = 1.550519768e-8;   /* conversion from TCB to TDB;  see */
                  /* https://github.com/IAU-ADES/ADES-Master/issues/52#issuecomment-2024035881 */

         snprintf( tag, sizeof( tag), "%s%c_gaia_geocentric",
                                          (i < 3 ? "" : "v"), 'x' + i % 3);
         value = atof( get_csv( tbuff, ibuff, header, tag, sizeof( tag) - 1));
         value *= AU_IN_KM;
         if( i < 3)
            printf( "    <pos%d>%.4f</pos%d>\n", i + 1, value * (1 - LB), i + 1);
         else
            printf( "    <vel%d>%.10f</vel%d>\n", i - 2, value / seconds_per_day, i - 2);
         }
      jd_obs = atof( get_csv( tbuff, ibuff, header, "epoch_utc,", 99)) + jd_epoch;
      full_ctime( tbuff, jd_obs, FULL_CTIME_MILLISECS | FULL_CTIME_YMD
                           | FULL_CTIME_LEADING_ZEROES | FULL_CTIME_MONTHS_AS_DIGITS);
      tbuff[4] = tbuff[7] = '-';
      tbuff[10] = 'T';
      printf( "    <obsTime>%sZ</obsTime>\n", tbuff);
      printf( "    <ra>%.9f</ra>\n",           atof( get_csv( tbuff, ibuff, header, "ra,", 99)));
      printf( "    <dec>%.9f</dec>\n",         atof( get_csv( tbuff, ibuff, header, "dec,", 99)));
      printf( "    <rmsRA>%.5f</rmsRA>\n",     atof( get_csv( tbuff, ibuff, header, "ra_error_random,", 99)) / 1000.);
      printf( "    <rmsDec>%.5f</rmsDec>\n",   atof( get_csv( tbuff, ibuff, header, "dec_error_random,", 99)) / 1000.);
      printf( "    <rmsCorr>%s</rmsCorr>\n", get_csv( tbuff, ibuff, header, "ra_dec_correlation_random,", 10));
      printf( "    <astCat>Gaia3</astCat>\n");
      if( strstr( header, "g_mag,"))   /* FPR doesn't have magnitude data */
         {
         get_csv( tbuff, ibuff, header, "g_mag,", 9);
         if( *tbuff >= '0' && *tbuff <= '9')
            {
            double g_flux, g_flux_sigma;

            printf( "    <mag>%s</mag>\n", tbuff);
            g_flux       = atof( get_csv( tbuff, ibuff, header, "g_flux", 99));
            g_flux_sigma = atof( get_csv( tbuff, ibuff, header, "g_flux_error", 99));
            printf( "    <rmsMag>%.3f</rmsMag>\n",
                              (g_flux_sigma / g_flux) / log( 2.512));
            printf( "    <band>G</band>\n");
            printf( "    <photCat>Gaia2</photCat>\n");
            }
         }
      printf( "   </optical>\n");
      }
}

static void error_exit( void)
{
   fprintf( stderr, "Run as 'csv2ades (asteroid number)' for ADES data on a single\n"
            "object,  or 'csv2ades (lownum,highnum)' for ADES data for a range\n"
            "of numbered objects.  Output is to stdout.  See 'csv2ades.c' for\n"
            "information about where to get the Gaia data in CSV files,  needed\n"
            "as input to this program.\n");
   exit( -1);
}

/* See above comments about the layout of numbered observations in the
various SsoObservation_XX.csv files.  Unless you ask for more than 20
asteroids,  we don't need to look at all twenty files.      */

static bool file_needed( const int file_no, const int low_num, const int high_num)
{
   int i;

   for( i = low_num; i <= high_num && i < low_num + 20; i++)
      if( ((i + 11) % 20) == file_no)
         return true;
   return( false);
}

int main( const int argc, const char **argv)
{
   FILE *hdr_file = fopen( "gaia.hdr", "rb");
   char buff[800];
   int low_num, high_num = -1, file_num;

   if( !hdr_file)
      {
      perror( "Can't open 'gaia.hdr'");
      return( -1);
      }
   if( argc != 2 || sscanf( argv[1], "%d,%d", &low_num, &high_num) < 1)
      {
      fprintf( stderr, "No asteroids specified on the command line\n");
      error_exit( );
      }
   if( high_num < 0)    /* just specified one asteroid */
      high_num = low_num;
   while( fgets( buff, sizeof( buff), hdr_file))
      if( *buff != '*')
         printf( "%s", buff);
      else for( file_num = 0; file_num < 20; file_num++)
         if( file_needed( file_num, low_num, high_num))
            {
            FILE *ifile;
            char header[800];

            snprintf( buff, sizeof( buff), "SsoObservation_%02d.csv", file_num);
            ifile = fopen( buff, "rb");
            if( !ifile)
               {
               fprintf( stderr, "File '%s' not found\n"
                     "'csv2ades.c' contains directions on where to get Gaia data\n", buff);
               exit( -2);
               }
            *header = '\0';
            while( fgets( buff, sizeof( buff), ifile))
               {
               if( !memcmp( buff, "solution_id", 11))
                  strcpy( header, buff);
               else if( *buff != '#')
                  csv_to_ades( buff, header, low_num, high_num);
               }
            fclose( ifile);
            }
   fclose( hdr_file);
   return( 0);
}

/*
   <optical>
    <permID>3200</permID>        number_mp
    <mode>TDI</mode>
    <stn>258</stn>
    <sys>ICRF_AU</sys>
    <ctr>0</ctr>
    <pos1>-0.9241930322</pos1>    x_gaia
    <pos2>+0.3423918499</pos2>    y_gaia
    <pos3>+0.1490296836</pos3>    z_gaia
    <obsTime>2016-02-27T15:45:22.593Z</obsTime>     ! epoch_utc (days since 2016 Feb 27 = JD 2457445.5)
    <ra>29.28052137</ra>                              ra
    <dec>+21.81708923</dec>                           dec
    <rmsRA>0.30945110</rmsRA>                       ! ra_error_random (IN milliarcsec)
    <rmsDec>0.52852813</rmsDec>                     ! dec_error_random
    <rmsCorr>0.9999607648</rmsCorr>                   ra_dec_correlation_random
    <astCat>Gaia2</astCat>
    <exp>4.4</exp>
    <remarks>TCB=2016-02-27T15:46:49.938Z, transitId=47101903417545381, observationId=471019034175453811, rmsRaSyst=0.00848502 as, rmsDecSyst=0.00394375 as, rmsCorrSyst=0.88192202, positionScanAngle=300.348 deg</remarks>
   </optical>

solution_id,source_id,observation_id,number_mp,         0-3
epoch,epoch_err,epoch_utc,ra,                           4-7
dec,ra_error_systematic,dec_error_systematic,ra_dec_correlation_systematic,   8-11
ra_error_random,dec_error_random,ra_dec_correlation_random,g_mag,           12-15
g_flux,g_flux_error,x_gaia,y_gaia,                                          16-19
z_gaia,vx_gaia,vy_gaia,vz_gaia,                                             20-23
position_angle_scan,level_of_confidence                                     24-26
*/
