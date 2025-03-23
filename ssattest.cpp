/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "watdefs.h"
#include "lunar.h"
#include "date.h"
#include "afuncs.h"

#define JUPITER_R (71492. / AU_IN_KM)
#define J2000     2451545.

int main( const int argc, const char **argv)
{
   int i;
   double loc[12], jd, precess_matrix[9], t_years, obliquity;
   char buff[80];

   if( argc != 2)
      {
      printf( "'ssattest' takes a JD on the command line,  and outputs\n"
              "J2000 ecliptic Cartesian coordinates for the eight main\n"
              "satellites of Saturn and the four Galileans.\n");
      return( -1);
      }
   jd = get_time_from_string( 0., argv[1], 0, NULL);
   full_ctime( buff, jd, FULL_CTIME_YMD);
   t_years = (jd - J2000) / 365.25;
   obliquity = mean_obliquity( t_years / 100.);
   printf( "Coordinates are in AU,  ecliptic J2000\n");
   printf( "Date: %.5f = %s\n", jd, buff);
   for( i = 0; i < 8; i++)
      {
      const char *names[8] = {  "Mimas", "Enceladus", "Tethys",
                  "Dione", "Rhea", "Titan", "Hyperion", "Iapetus" };

      calc_ssat_loc( jd, loc, i, 0L);
      printf( "%d: %10.7f %10.7f %10.7f %s\n", i, loc[0], loc[1], loc[2],
                     names[i]);
      }
   printf( "\n");
   setup_precession( precess_matrix, 2000. + t_years, 2000.);
   calc_jsat_loc( jd, loc, 15, 0L);
   for( i = 0; i < 4; i++)
      {
      const char *names[4] = {  "Io", "Europa", "Ganymede", "Callisto" };
      double tloc[3];

                        /* turn ecliptic of date to equatorial: */
      rotate_vector( loc + i * 3, obliquity, 0);
                        /* then to equatorial J2000: */
      precess_vector( precess_matrix, loc + i * 3, tloc);
                        /* then to ecliptic J2000: */
      equatorial_to_ecliptic( tloc);
      printf( "%d: %10.7f %10.7f %10.7f %s\n", i, tloc[0] * JUPITER_R,
               tloc[1] * JUPITER_R, tloc[2] * JUPITER_R, names[i]);
      }
   return( 0);
}
