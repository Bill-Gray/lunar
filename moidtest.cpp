#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "comets.h"
#include "afuncs.h"
#include "mpc_func.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define J2000 2451545.0

static bool verbose = false;

static void show_elements( const ELEMENTS *elem)
{
#if 0
   printf( "q=%13.10f e=%12.10f i=%10.6f asc_node=%10.6f arg_per=%10.6f\n",
#endif
   printf( "%13.10f %12.10f %10.6f %10.6f %10.6f\n",
               elem->q, elem->ecc,
               elem->incl * 180. / PI,
               elem->asc_node * 180. / PI,
               elem->arg_per * 180. / PI);
}

                        /* `Oumuamua : MOID = 0.0958212 */
                        /* (433) Eros : MOID = 0.149341 */
                        /* (4179) Toutatis : MOID = 0.00609937 */
                        /* (588) Achilles : MOID = 3.47306 */
                        /* (163693) Atira : MOID = 0.206823 */
                        /* C/1995 O1 (Hale-Bopp) : MOID = 0.116025 */
                        /* 1999 AN10 = (137108) : MOID = 0.000128404  */

static int run_intraobject_check( const char *header, FILE *ifile)
{
   const long loc = ftell( ifile);
   int i, n_objects = 0;
   char buff[300];
   ELEMENTS *elem, telem;
   double max_diff = 0.;

   while( fgets( buff, sizeof( buff), ifile))
      if( !extract_sof_data( &telem, buff, header))
         n_objects++;
   printf( "%d objects\n", n_objects);
   elem = (ELEMENTS *)calloc( n_objects, sizeof( ELEMENTS));
   fseek( ifile, loc, SEEK_SET);
   n_objects = 0;
   while( fgets( buff, sizeof( buff), ifile))
      if( !extract_sof_data( &elem[n_objects], buff, header))
         {
         if( verbose)
            printf( "%d %.30s\n", n_objects, buff);
         derive_quantities( &elem[n_objects], SOLAR_GM);
         for( i = 0; i < n_objects; i++)
            {
            moid_data_t mdata;
#ifdef TEMP_REMOVE
            double moid1 = find_moid_full( &elem[i], &elem[n_objects], &mdata);
            double moid2 = find_moid_full( &elem[n_objects], &elem[i], &mdata);
            double diff = fabs( moid1 - moid2);

            if( verbose)
               printf( "%d %d  %.13lf %.13lf %.13lf\n", i, n_objects, moid1,
                        moid2, diff);
            if( max_diff < diff)
               max_diff = diff;
#endif
            double moid;

            if( elem[n_objects].ecc < elem[i].ecc)
               moid = find_moid_full( &elem[i], &elem[n_objects], &mdata);
            else
               moid = find_moid_full( &elem[n_objects], &elem[i], &mdata);
            printf( "MOID = %.12f\n", moid);
            }
         n_objects++;
         }
   printf( "%d objects\n", n_objects);
   free( elem);
   fclose( ifile);
   printf( "Max diff %.13f\n", max_diff);
   return( 0);
}

int main( const int argc, const char **argv)
{
   ELEMENTS elem, earth_elem;
   double t_cen, moid, reference_moid = 0.;
   const char *ifilename = "moidtest.txt";
   FILE *ifile;
   char header[300], buff[300], obj_name[60];
   int i, planet_number = 3;
   bool elems_found = false, reversing = false, intraobject_check = false;
   moid_data_t mdata;

   *obj_name = '\0';
   for( i = 1; i < argc; i++)
      if( argv[i][0] != '-')
         {
         if( *obj_name)
            strcat( obj_name, " ");
         strcat( obj_name, argv[i]);
         }
      else
         {
         const char *arg = (i == argc - 1 || argv[i][2]) ? argv[i] + 2 :
                                    argv[i + 1];

         switch( argv[i][1])
            {
            case 'i':
               intraobject_check = true;
               break;
            case 'f':
               ifilename = arg;
               break;
            case 'n':
               planet_number = atoi( arg);
               break;
            case 'r':
               reversing = true;
               break;
            case 'v':
               verbose = true;
               break;
            default:
               fprintf( stderr, "Command-line option '%s' not recognized\n",
                           argv[i]);
               break;
            }
         }
   ifile = fopen( ifilename, "rb");
   assert( ifile);
   if( !fgets( header, sizeof( header), ifile))
      {
      fprintf( stderr, "Didn't get test file header line\n");
      return( -1);
      }
   memset( &elem, 0, sizeof( ELEMENTS));
   if( intraobject_check)
      return( run_intraobject_check( header, ifile));
   while( fgets( buff, sizeof( buff), ifile))
      if( strstr( buff, obj_name))
         {
         if( !memcmp( buff, "MOID", 4))
            reference_moid = atof( buff + 26);
         else if( !extract_sof_data( &elem, buff, header))
            elems_found = true;
         }
   fclose( ifile);
   if( !elems_found)
      {
      fprintf( stderr, "No elements read for '%s'\n", obj_name);
      return( -1);
      }
   t_cen = (elem.epoch - J2000) / 36525.;
   setup_planet_elem( &earth_elem, planet_number, t_cen);
   derive_quantities( &elem, SOLAR_GM);
   if( reversing)
      {
      ELEMENTS telem = elem;

      elem = earth_elem;
      earth_elem = telem;
      }
   moid = find_moid_full( &earth_elem, &elem, &mdata);
   printf( "Barbee speed %f km/s\n", mdata.barbee_speed * AU_IN_KM / seconds_per_day);
   printf( "JD closest approach %f\n", mdata.jd1);
   printf( "MOID = %.13f\n", moid);
   if( reference_moid)
      printf( "Reference MOID = %.13f;  diff %.13f\n",
            reference_moid, reference_moid - moid);
   show_elements( &elem);
   show_elements( &earth_elem);
   return( 0);
}
