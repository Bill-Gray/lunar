#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpc_func.h"

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

   assert( ifile);
   assert( ades_context);
   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'd':
               show_data = 1;
               break;
            default:
               fprintf( stderr, "'%s' not recognized\n", argv[i]);
               break;
            }
   while( fgets_with_ades_xlation( buff, sizeof( buff), ades_context, ifile))
      {
      printf( "%s\n", buff);
      if( show_data)
         {
         const double jd = extract_date_from_mpc_report( buff, NULL);
         double ra, dec;

         if( jd)
            {
            get_ra_dec_from_mpc_report( buff, NULL, &ra, NULL,
                                              NULL, &dec, NULL);
            printf( "MJD %f  RA %f  dec %f\n", jd - 2400000.5, ra, dec);
            }
         }
      }
   rval = free_ades2mpc_context( ades_context);
   fclose( ifile);
   printf( "rval = %d\n", rval);
   return( 0);
}
