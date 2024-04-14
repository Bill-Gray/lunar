#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "watdefs.h"
#include "date.h"
#include "jpl_xref.h"

/* You can access vector (and other) ephemerides from JPL Horizons
with a URL.  This program generates the URL required for a particular
object covering a particular year at 0.1 day intervals;  the resulting
ephemerides can be fed through 'jpl2mpc',  then 'eph2tle' to generate
TLEs for that year.  The URL looks like (broken into multiple lines
for clarity/explanations;  this is for Gaia in 2020):

https://ssd.jpl.nasa.gov/api/horizons_api?format=text
&COMMAND=%27-139479%27            object -139479 = JPL ID for Gaia
&OBJ_DATA=%27YES%27               give us the header/boilerplate info
&TABLE_TYPE=%27V%27               V = vector ephems
&START_TIME=%272020-01-01%27      starts 2020 Jan 1
&START_TIME=%27JD2458849.5%27     (alternative start format)
&STOP_TIME=%272021-01-01%27       ends the following year
&STOP_TIME=%27JD2459215.5%27      (similarly alternative end format)
&STEP_SIZE=%273660%27             3660 steps (2020 was a leap year)
&VEC_TABLE=%272%27                table type 2 = positions & velocities
&VEC_LABELS=%27NO%27              skip the XYZ VXVYVZ labelling

Compiles with gcc -Wall -Wextra -pedantic -o jpl_url jpl_url.c liblunar.a -lm */

static void error_exit( void)
{
   const size_t n_xrefs = sizeof( jpl_xrefs) / sizeof( jpl_xrefs[0]);
   size_t i;

   fprintf( stderr, "Usage : jpl_url <MPC code or YYYY-NNNA desig> <start date> <end date>\n"
           "See comments in 'jpl_url.c' for details.\n");
   for( i = 0; i < n_xrefs; i++)
      fprintf( stderr, "%s %-11s %s\n", jpl_xrefs[i].mpc_code,
                                        jpl_xrefs[i].intl_desig,
                                        jpl_xrefs[i].name);
   exit( -1);
}

int main( const int argc, const char **argv)
{
   size_t i = 0;
   const size_t n_xrefs = sizeof( jpl_xrefs) / sizeof( jpl_xrefs[0]);
   double jd1, jd2;
   const double min_jd = 2447892.5;    /* 1990 Jan 1 */

   if( argc < 4)
      error_exit( );
   while( i < n_xrefs && strcmp( argv[1], jpl_xrefs[i].mpc_code)
                      && strcmp( argv[1], jpl_xrefs[i].intl_desig)
                      && atoi( argv[1]) != jpl_xrefs[i].jpl_desig)
      i++;
   if( i == n_xrefs)
      {
      fprintf( stderr, "'%s' is not a spacecraft MPC code,  YYYY-NNNA"
                        " designation,  or Horizons object number\n", argv[1]);
      fprintf( stderr, "Valid designations are :\n");
      error_exit( );
      }
   jd1 = get_time_from_string( 0., argv[2], 0, NULL);
   jd2 = get_time_from_string( jd1, argv[3], 0, NULL);
   if( jd1 < min_jd || jd2 < min_jd)
      {
      fprintf( stderr, "Can't go before JD %f\n", min_jd);
      error_exit( );
      }
   printf( "https://ssd.jpl.nasa.gov/api/horizons.api?format=text");
   printf( "&COMMAND='%d'", jpl_xrefs[i].jpl_desig);
   printf( "&OBJ_DATA='YES'&TABLE_TYPE='V'");
   printf( "&START_TIME='JD%f'", jd1);
   printf( "&STOP_TIME='JD%f'", jd2);
   printf( "&STEP_SIZE='%d'", (int)( (jd2 - jd1) * 10. + 0.5));
   printf( "&VEC_TABLE='2'&VEC_LABELS='NO'\n");
   return( 0);
}
