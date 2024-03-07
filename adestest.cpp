#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpc_func.h"
#include "watdefs.h"
#include "date.h"

/* Given the name of a file containing XML or PSV ADES data as a command
line argument,  'adestest' will read it and output 80-column MPC-like
astrometry.  This can be used to make XML ADES slightly more readable
(and PSV ADES data slightly _less_ readable),  but its main purpose
was to let me test my code for parsing ADES data.  It also does serve
as an example of how to use the 'ades2mpc.cpp' ADES-parsing functions.

Note that the 80-column output contains various extensions to the MPC
format.  RA/decs are stored in decimal degrees,  both to match what
we get from ADES and to allow greater precision.  Uncertainties (which
the 80-column format knows nothing about) are stored in COM (comment)
lines.  Dates/times are stored in a compacted form to allow millisecond
precision (the usual MPC format allows only 10^-6 day = 86.4 ms
precision).  The resulting "80-column data" will work with all of my
tools,  but probably not with anyone else's. */

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( argv[1], "rb");
   char buff[400];
   void *ades_context = init_ades2mpc( );
   int i, rval, show_data = 0;
   bool make_mpc80 = false, comments = true;

   assert( ifile);
   assert( ades_context);
   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'd':
               show_data = 1;
               break;
            case 'm':
               make_mpc80 = true;
               break;
            case 'c':
               comments = false;
               break;
            default:
               fprintf( stderr, "'%s' not recognized\n", argv[i]);
               break;
            }
   while( fgets_with_ades_xlation( buff, sizeof( buff), ades_context, ifile))
      {
      if( show_data || make_mpc80)
         {
         unsigned time_format;
         const double jd = extract_date_from_mpc_report( buff, &time_format);
         double ra, dec;

         if( jd)
            {
            int ra_format, dec_format;
            double ra_prec, dec_prec;
            const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923;

            get_ra_dec_from_mpc_report( buff, &ra_format, &ra, &ra_prec,
                                             &dec_format, &dec, &dec_prec);

            if( show_data)
               printf( "MJD %f  RA %f  dec %f\n", jd - 2400000.5, ra, dec);
            if( make_mpc80)
               {                 /* switch date, RA, dec to MPC80 format */
               char tbuff[40];
               int format = FULL_CTIME_YMD | FULL_CTIME_LEADING_ZEROES
                                 | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_FORMAT_DAY;

               format |= FULL_CTIME_6_PLACES;
               full_ctime( tbuff, jd, format);
               memcpy( buff + 15, tbuff, 17);
               output_angle_to_buff( tbuff, ra * 12. / PI, 3);
               memcpy( buff + 32, tbuff, 12);
               output_signed_angle_to_buff( tbuff, dec * 180. / PI, 2);
               memcpy( buff + 44, tbuff, 12);
               }
            }
         }
      if( comments || memcmp( buff, "COM ", 4))
         printf( "%s\n", buff);
      }
   rval = free_ades2mpc_context( ades_context);
   fclose( ifile);
   printf( "rval = %d\n", rval);
   return( 0);
}
