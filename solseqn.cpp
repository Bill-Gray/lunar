/* solseqn.cpp: computes date/times of solstices and equinoxes

Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

/* Some test code I wrote a while back to compute dates/times for
solstices and equinoxes.  I don't actually use it anywhere.  You
can disregard it freely.   */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"
#include "date.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

double get_solstice_equinox_date( double jd);         /* solseqn.cpp */

static double daily_variation( const double t_cen)
{
   double rval = 3548.193;
   int i;
   static const double constants[34] = {
           7.311, 333.4515, 359993.7286,  /* tau */
            .305, 330.9814, 719987.4571,  /* tau */
            .010, 328.5170, 1079981.1857, /* tau */
         118.568,  87.5287, 359993.7286,
           2.476,  85.0561, 719987.4571,
           1.376,  27.8502, 4452671.1152,
            .119,  73.1375, 450368.8564,
            .114, 337.2264, 329644.6718,
            .086, 222.5400, 659289.3436,
            .078, 162.8136, 9224659.7915,
            .054,  82.5823, 1079981.1857,
            0. };
   for( i = 0; constants[i] != 0.; i += 3)
      {
      const double arg = constants[i + 1] + constants[i + 2] * t_cen / 10.;
      double term = sin( arg * PI / 180.);

      if( i < 9)
         term *= t_cen / 10.;
      rval += constants[i] * term;
      }
   return( rval);
}

double get_solstice_equinox_date( double jd)
{
   double err = 1.;

   while( fabs( err) > .000001)
      {
      const double t_c = (jd - 2451545.) / 36525.;
      double lon, ovals[3], d_lon;

      if( calc_big_vsop_loc( NULL, 3, ovals, t_c, 0.))
         return( 0.);         /* failed to find the 'big_vsop.bin' file */
      lon = ovals[0] - .09033 * (PI / 180.) / 3600.;
      nutation( t_c, &d_lon, NULL);
      lon += d_lon * (PI / 180.) / 3600.;
                  /* now include aberration: */
      lon -= .005775518 * ovals[2] *
                        daily_variation( t_c) * (PI / 180.) / 3600.;
                  /* now convert to degrees: */
      lon *= 180. / PI;
      while( lon > 45.)
         lon -= 90.;
      while( lon < -45.)
         lon += 90.;
      jd -= lon * 365.25 / 360.;
      err = lon;
      }
   return( jd);
}

#ifdef TEST_CODE
int main( const int argc, const char **argv)
{
   double year = atof( argv[1]), jd = 1721139.29189 + 365.242137 * year;
   double max_jd = 1721139.29189 + 365.242137 * 10000.;
               /* jd is approx March equinox */
   int season = 0, i, verbose = 0, output_mask = 15;
   int use_ut = 0;
   FILE *output_file = NULL;

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'u':
               use_ut = 1;
               break;
            case 'm':
               max_jd = 1721139.29189 + 365.242137 * atof( argv[i] + 2);
               break;
            case 'v':
               verbose = 1 + atoi( argv[i] + 2);
               break;
            case 'o':
               output_mask = atoi( argv[i] + 2);
               break;
            case 'f':
               output_file = fopen( argv[i] + 2, "w");
               break;
            default:
               printf( "'%s' is not a recognized option\n", argv[i]);
               break;
            }

   setvbuf( stdout, NULL, _IONBF, 0);
   if( verbose > 1)
      printf( "Starting up...\n");

   while( jd < max_jd)
      {
      char time_buff[80];
      const char *season_text[4] =
                     { "Mar equ", "Jun sol", "Sep equ", "Dec sol" };

      if( verbose > 1)
         printf( "Iterating...\n");
      jd = get_solstice_equinox_date( jd);
      if( !jd)
         {
         printf( "Couldn't find 'big_vsop.bin'\n");
         exit( -1);
         }
      if( verbose > 2)
         printf( "   full_ctime\n");
      if( use_ut)
         jd -= td_minus_ut( jd) / 86400.;
      full_ctime( time_buff, jd, 0);
      if( (output_mask >> season) & 1)
         {
         if( output_file)
            fprintf( output_file, "%s: %s %.5lf\n", season_text[season], time_buff, jd);
         printf( "%s: %s %.5lf\n", season_text[season], time_buff, jd);
         }
      season = (season + 1) % 4;
      jd += 91.3;        /* advance to "mean" next season */
      }
   if( output_file)
      fclose( output_file);
   return( 0);
}
#endif
