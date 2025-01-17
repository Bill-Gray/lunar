#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpc_func.h"
#include "watdefs.h"
#include "date.h"

static void error_exit( )
{
   fprintf( stderr,
      "\n"
      "Given the name of a file containing XML or PSV ADES data as a command\n"
      "line argument,  'adestest' will read it and output 80-column MPC-like\n"
      "astrometry.  This can be used to make XML ADES slightly more readable\n"
      "(and PSV ADES data slightly _less_ readable),  but its main purpose\n"
      "was to let me test my code for parsing ADES data.  It also serves\n"
      "as an example of how to use the 'ades2mpc.cpp' ADES-parsing functions.\n"
      "\n"
      "Note that the 80-column output contains various extensions to the MPC\n"
      "format.  RA/decs are stored in decimal degrees,  both to match what\n"
      "we get from ADES and to allow greater precision.  Uncertainties (which\n"
      "the 80-column format knows nothing about) are stored in COM (comment)\n"
      "lines.  Dates/times are stored in a compacted form to allow millisecond\n"
      "precision (the usual MPC format allows only 10^-6 day = 86.4 ms\n"
      "precision).  The resulting '80-column data' will work with all of my\n"
      "tools,  but probably not with anyone else's.\n"
      "\n"
      "    Add a '-m' command line switch to get 'real' 80-column data,  with\n"
      "times in decimal days and RA/decs in base-60 form.\n");
   exit( -1);
}

void ades_artsat_desigs( void *ades_context, const bool ignore_artsat_desigs);

int main( const int argc, const char **argv)
{
   FILE *ifile;
   char buff[400];
   void *ades_context = init_ades2mpc( );
   int i, rval, show_data = 0;
   bool make_mpc80 = false, comments = true;

   if( argc < 2)
      {
      fprintf( stderr, "No file of ADES astrometry supplied on command line\n");
      error_exit( );
      }
   ifile = fopen( argv[1], "rb");
   if( !ifile)
      {
      fprintf( stderr, "Couldn't open '%s' : ", argv[1]);
      perror( NULL);
      error_exit( );
      }
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
            case 'a':
               ades_artsat_desigs( ades_context, true);
               break;
            default:
               fprintf( stderr, "'%s' not recognized\n", argv[i]);
               error_exit( );
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
               if( buff[14] != 'v')
                  {
                  output_angle_to_buff( tbuff, ra * 12. / PI, 3);
                  memcpy( buff + 32, tbuff, 12);
                  output_signed_angle_to_buff( tbuff, dec * 180. / PI, 2);
                  memcpy( buff + 44, tbuff, 12);
                  }
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
