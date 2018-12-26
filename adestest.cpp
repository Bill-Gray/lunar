#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpc_func.h"

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
   printf( "rval = %d\n", rval);
   return( 0);
}
