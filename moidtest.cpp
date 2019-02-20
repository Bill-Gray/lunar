#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "brentmin.h"
#include "comets.h"
#include "afuncs.h"
#include "mpc_func.h"

int extract_sof_data( ELEMENTS *elem, const char *buff, const char *header);

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define J2000 2451545.0

static void show_elements( const ELEMENTS *elem)
{
   printf( "q=%8.5f e=%8.6f i=%8.4f asc_node=%8.4f arg_per=%8.4f\n",
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

int main( const int argc, const char **argv)
{
   ELEMENTS elem, earth_elem;
   double t_cen, moid, reference_moid = 0.;
   FILE *ifile = fopen( "moidtest.txt", "rb");
   char header[300], buff[300], obj_name[60];
   int i;
   bool elems_found = false;
   moid_data_t mdata;

   assert( ifile);
   assert( argc > 1);
   if( !fgets( header, sizeof( header), ifile))
      {
      fprintf( stderr, "Didn't get test file header line\n");
      return( -1);
      }
   memset( &elem, 0, sizeof( ELEMENTS));
   strcpy( obj_name, argv[1]);
   for( i = 2; i < argc; i++)
      if( argv[i][0] != '-')
         {
         strcat( obj_name, " ");
         strcat( obj_name, argv[i]);
         }
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
   setup_planet_elem( &earth_elem, 3, t_cen);
   derive_quantities( &elem, SOLAR_GM);
   for( i = 2; i < argc; i++)
      if( argv[i][0] != '-')
         switch( argv[i][1])
            {
            case 'r':
               {
               ELEMENTS telem = elem;

               elem = earth_elem;
               earth_elem = telem;
               }
               break;
            }
   moid = find_moid_full( &earth_elem, &elem, &mdata);
   printf( "Barbee speed %f km/s\n", mdata.barbee_speed * AU_IN_KM / seconds_per_day);
   printf( "JD closest approach %f\n", mdata.jd1);
   printf( "MOID = %.13f\n", moid);
   if( reference_moid)
      printf( "Reference MOID = %.13f;  diff %.13f\n",
            reference_moid, reference_moid - moid);
   show_elements( &elem);
   return( 0);
}
