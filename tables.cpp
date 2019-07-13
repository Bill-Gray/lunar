/* tables.cpp: prints rise/set "almanac"-type tables

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

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "watdefs.h"
#include "lunar.h"
#include "date.h"
#include "afuncs.h"
#include "riseset3.h"

const static double pi =
     3.1415926535897932384626433832795028841971693993751058209749445923078;

static void get_rise_set_times( double *rise_set, const int planet_no,
                  double jd,
                  const double observer_lat, const double observer_lon,
                  const char *vsop_data)
{
   int i;

                                    /* Mark both the rise and set times     */
                                    /* as -1,  to indicate that they've     */
                                    /* not been found.  Of course,  it may  */
                                    /* turn out that one,  or both,  do     */
                                    /* not occur during the given 24 hours. */
   rise_set[0] = rise_set[1] = -1;
                                    /* Compute the altitude for each hour:  */
   for( i = 0; i < 24; i++)
      {
      int idx;
      double jd_riseset;

      jd_riseset = look_for_rise_set( planet_no, jd, jd + 1. / 24.,
                  observer_lat, observer_lon, vsop_data, &idx);

      if( idx != -1)
         rise_set[idx] = jd_riseset;
      jd += 1. / 24.;
      }
}

#if defined(_MSC_VER) && _MSC_VER < 1900
                      /* For older MSVCs,  we have to supply our own  */
                      /* snprintf().  See snprintf.cpp for details.  */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif


   /* The 'quadrant' function helps in figuring out dates of lunar phases
and solstices/equinoxes.  If the solar longitude is in one quadrant at
the start of a day,  but in a different quadrant at the end of a day,
then we know that there must have been a solstice or equinox during that
day.  Also,  if (lunar longitude - solar longitude) changes quadrants
from the start of a day to the end of a day,  we know there must have
been a lunar phase change during that day.

   In this code,  I don't bother finding the exact instant of these
events.  The code just checks for a quadrant change and reports the
corresponding event. */

static int quadrant( double angle)
{
   angle = fmod( angle, 2. * pi);
   if( angle < 0.)
      angle += 2. * pi;
   return( (int)( angle * 2. / pi));
}

int main( int argc, char **argv)
{
   char *vsop_data = load_file_into_memory( "vsop.bin", NULL);
   int i, year = atoi( argv[1]);
   int month_start = 1, month_end = 12, month;
   const double observer_lon = -69.90 * pi / 180.;
   const double observer_lat = 44.01 * pi / 180.;

   if( !vsop_data)
      {
      printf( "VSOP.BIN wasn't loaded.\n");
      return( -1);
      }

   if( argc > 2)        /* month specified,  rather than "entire year" */
      month_start = month_end = atoi( argv[2]);

   printf( "       Sun          Moon\n");
   printf( "Day  Rise Set     Rise Set\n");
   for( month = month_start; month <= month_end; month++)
      {
      long jd_start, jd_end;
      const int time_zone = -5;

      jd_start = dmy_to_day( 1, month, year, 0);
      if( month == 12)
         jd_end = dmy_to_day( 1, 1, year + 1, 0);
      else
         jd_end = dmy_to_day( 1, month + 1, year, 0);


      for( i = 0; i < (int)( jd_end - jd_start); i++)
         {
         double rise_set[4];
         double lunar_lon[2], solar_lon[2];
         double jd = (double)( jd_start + i) - .5 - (double)time_zone / 24.;
         char buff[80];
         int j, quad0, quad1;

         memset( buff, 0, 40);
         assert( i <= 30);
         get_rise_set_times( rise_set, 3,  jd, observer_lat, observer_lon,
                                                                vsop_data);
         get_rise_set_times( rise_set + 2, 10, jd, observer_lat, observer_lon,
                                                                vsop_data);
         if( (jd_start + i) % 7 == 6)        /* Sunday */
            strcpy( buff, "Su");
         else
            snprintf( buff, 3, "%2d", i + 1);
         for( j = 0; j < 4; j++)
            {
            static const int offsets[5] = { 4, 10, 17, 23, 29 };

            if( rise_set[j] < 0.)
               strcpy( buff + offsets[j], "--:--");
            else
               {
               unsigned minutes;
               double fraction;

               fraction = rise_set[j] + .5 + (double)time_zone / 24.;
               minutes = (unsigned)( (fraction - floor( fraction)) * 1440.0);

               snprintf( buff + offsets[j], 6, "%02u:%02u",
                                   (minutes / 60) % 24, minutes % 60);
               }
            }

         for( j = 0; j < 2; j++)
            {
            PLANET_DATA pdata;

            fill_planet_data( &pdata, 3, jd + (double)j,
                         observer_lat, observer_lon, vsop_data);
            solar_lon[j] = pdata.ecliptic_lon;
            fill_planet_data( &pdata, 10, jd + (double)j,
                         observer_lat, observer_lon, vsop_data);
            lunar_lon[j] = pdata.ecliptic_lon;
            }

         quad1 = quadrant( lunar_lon[1] - solar_lon[1]);
         quad0 = quadrant( lunar_lon[0] - solar_lon[0]);
         if( quad1 != quad0)
            {
            const char *phase_names = "1Q FM 3Q NM";

            memcpy( buff + 29, phase_names + quad0 * 3, 2);
            }
         quad1 = quadrant( solar_lon[1]);
         quad0 = quadrant( solar_lon[0]);
         if( quad1 != quad0)
            {
            static const char *strings[4] =
                            { "Summ Sol", "Autu Eq", "Wint Sol", "Vern Eq" };

            strcpy( buff + 29, strings[quad0]);
            }

         for( j = 0; j < 39; j++)
            if( !buff[j])
               buff[j] = ' ';
         printf( "%s\n", buff);
         }
      }
   free( vsop_data);
   return( 0);
}
