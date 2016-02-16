/* conjunct.cpp: finds dates/times of planetary conjunctions

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <math.h>
#include <time.h>
#include "watdefs.h"
#include "date.h"
#include "lunar.h"
#include "jpleph.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define AU_PER_CENTURY (AU_PER_DAY * 36525.)
#define N_OBJECTS 11

FILE *log_file = NULL;

static char is_used[N_OBJECTS];
static int home_planet = 3, occultations_only = 0;

static void get_rectang( void *jpleph, double *loc, double t_c, int planet_no)
{
   if( !planet_no)
      loc[0] = loc[1] = loc[2] = 0.;
   else
      {
      double tloc[6];
//    static int count = 100;

//    if( count)
//       {
//       count--;
//       printf( "t: %.5lf pl %d: ", t_c, planet_no);
//       }
      jpl_pleph( jpleph, t_c * 36525. + 2451545., planet_no, 11, tloc, 0);
//    if( count)
//       {
//       count--;
//       printf( "got it\n");
//       }
      memcpy( loc, tloc, 3 * sizeof( double));
      }
}

static void find_posns( void *jpleph, double *loc, double t, int n_passes,
                              const char *used)
{
   double t_c = (t - 2000.) / 100., earth_loc[3];
   int i, j, pass;
   static double first_dist[N_OBJECTS];

   get_rectang( jpleph, earth_loc, t_c, home_planet);
   for( i = 0; i < N_OBJECTS; i++)
      if( used[i])
         {
         double xyz[3];
         double r, *tptr = loc + i * 3;

         r = first_dist[i];
         for( pass = 0; pass < n_passes; pass++)
            {
            if( i)
               get_rectang( jpleph, xyz, t_c - (r / AU_PER_CENTURY), i);
            else
               xyz[0] = xyz[1] = xyz[2] = 0.;
            r = 0.;
            for( j = 0; j < 3; j++)
               {
               xyz[j] -= earth_loc[j];
               r += xyz[j] * xyz[j];
               }
            r = sqrt( r);
            }
         tptr[0] = atan2( xyz[1], xyz[0]);
         tptr[1] = asin( xyz[2] / r);
         tptr[2] = r;
         first_dist[i] = r;
         }
}

#ifdef UNUSED
static double centralize( double ang, double center)
{
   while( ang - center > PI)
      ang -= PI + PI;
   while( ang - center < -PI)
      ang += PI + PI;
   return( ang);
}
#endif

static double rough_dist2( const double *p1, const double *p2)
{
   double dlon = p1[0] - p2[0], dlat = p1[1] - p2[1];

   while( dlon > PI)
      dlon -= PI + PI;
   while( dlon < -PI)
      dlon += PI + PI;
   dlon *= cos( (p1[1] + p2[1]) * .5);
   return( dlon * dlon + dlat * dlat);
}

static void search_conjunctions( void *jpleph,
              const double t0, const double step_size0,
              const double *loc1, const double *loc2, const double *loc3)
{
   int i, j;

   for( i = 0; i < N_OBJECTS; i++)
      for( j = 0; j < i; j++)
         if( is_used[i] && is_used[j] && (j != 3 || i != 10))
            {
            double dsquared[3], occult_margin, inner_margin;
            int pass;
            double t = t0, step_size = step_size0, jd;
            char obuff[80];

            dsquared[0] = rough_dist2( loc1 + i * 3, loc1 + j * 3);
            dsquared[1] = rough_dist2( loc2 + i * 3, loc2 + j * 3);
            dsquared[2] = rough_dist2( loc3 + i * 3, loc3 + j * 3);
            if( dsquared[1] < dsquared[0] && dsquared[1] < dsquared[2])
               {
               char buff[80];
               double loc[3 * N_OBJECTS], idist, jdist;
               static double planet_radii[N_OBJECTS] = { 696265., 2439.,
                        6051., 6378., 3393., 71398., 60000.,
                        25400., 25295., 1137., 1737.53 };
               static double planet_inner[N_OBJECTS] = { 696265., 2439.,
                        6051., 6378. * (1. - .00335281),
                               3393. * (1. - .0051865),
                              71398. * (1. - .0648088),
                              60000. * (1. - .1076209),
                              25400. * (1. - .03),
                              25295. * (1. - .022), 1137., 1737.53 };
               char t_used[N_OBJECTS];

               memset( t_used, 0, N_OBJECTS);
               t_used[home_planet] = t_used[i] = t_used[j] = 1;
               for( pass = 5; pass; pass--)
                  {
                  double minval = (dsquared[0] - dsquared[2]) /
                        (2. * (dsquared[2] + dsquared[0] - 2. * dsquared[1]));

                  t += (minval - 1.) * step_size;
                  step_size *= .1;
                  if( pass > 1)        /* still gotta refine estimate */
                     {
                     int k;

                     for( k = 0; k < 3; k++)
                        {
                        find_posns( jpleph, loc,
                                  t + (double)(k - 2) * step_size, 2, t_used);
                        dsquared[k] = rough_dist2( loc + i * 3, loc + j * 3);
                        }
                     }
//                if( i == 3 && j == 0)
//                   printf( "%d: %lf (%lf)\n", pass, t,
//                                sqrt( dsquared[1]) * 180. / PI);
                  }
               idist = loc[i * 3 + 2];
               jdist = loc[j * 3 + 2];
               jd = 2451545.0 + (t - 2000.) * 365.25;
               if( jdist > idist)
                  {
                  occult_margin = planet_radii[home_planet] + idist *
                             (planet_radii[j] - planet_radii[home_planet]) / jdist;
                  occult_margin = (occult_margin + planet_radii[i]) / idist;
                  inner_margin = planet_inner[home_planet] + idist *
                             (planet_inner[j] - planet_inner[home_planet]) / jdist;
                  inner_margin = (inner_margin + planet_inner[i]) / idist;
                  }
               else
                  {
                  occult_margin = planet_radii[home_planet] + jdist *
                             (planet_radii[i] - planet_radii[home_planet]) / idist;
                  occult_margin = (occult_margin + planet_radii[j]) / jdist;
                  inner_margin = planet_inner[home_planet] + jdist *
                             (planet_inner[i] - planet_inner[home_planet]) / idist;
                  inner_margin = (inner_margin + planet_inner[j]) / jdist;
                  }
               occult_margin /= AU_IN_KM;
               inner_margin /= AU_IN_KM;
               full_ctime( buff, jd, 0);
               dsquared[1] = sqrt( dsquared[1]);
               sprintf( obuff, "%2d %2d%7.3lf  %s", i, j,
                              (180. / PI) * dsquared[1], buff);
               if( occult_margin > dsquared[1])
                  {
                  if( j)      /* planet-planet event */
                     obuff[5] = 'M';
                  else if( idist < jdist)
                     obuff[5] = 'T';
                  if( inner_margin < dsquared[1])
                     obuff[6] = '?';         /* indicate 'possible' */
                  }
               if( !occultations_only || obuff[5] == 'M' || obuff[5] == 'T')
                  {
                  printf( "%s\n", obuff);
                  if( log_file)
                     {
                     obuff[13] = '^';
                     fprintf( log_file, "%s//tdj%12.4lf;gp%02d;gp%02d^\n", obuff,
                                          jd, i, j);
                     }
                  }
               }
            }
}

int main( int argc, char **argv)
{
   int i, loop_count = 0;
   double t = atof( argv[1]), delta_t = atof( argv[2]);
   double max_t = 99999.;
   double val1[3 * N_OBJECTS], val2[3 * N_OBJECTS], val3[3 * N_OBJECTS];
   const char *ephemeris_filename = "e:\\jpl_eph\\unix.406";
   void *jpleph;
   long seconds_elapsed = 0;

   memset( is_used, 1, N_OBJECTS);
   for( i = 3; i < argc; i++)
      switch( argv[i][0])
         {
         case 'e':
            ephemeris_filename = argv[i] + 1;
            break;
         case 'm':
            max_t = atof( argv[i] + 1);
            break;
         case 'o':
            occultations_only = 1;
            break;
         case 'u':
            {
            int j;

            for( j = 0; j < N_OBJECTS; j++)
               is_used[j] = (char)( argv[i][j + 1] == '1');
            }
            break;
         case 'v':
            home_planet = atoi( argv[i] + 1);
            break;
         }
   is_used[home_planet] = 0;

   jpleph = jpl_init_ephemeris( ephemeris_filename, NULL, NULL);
   if (!jpleph)
      {
      printf( "'%s' not opened\n", ephemeris_filename);
      return( -1);
      }

   log_file = fopen( "d:\\z2", "wb");

   while( !kbhit( ) && t < max_t)
      {
      find_posns( jpleph, val3, t, 1, is_used);
      if( loop_count > 2)
         search_conjunctions( jpleph, t, delta_t, val1, val2, val3);
      memcpy( val1, val2, 3 * N_OBJECTS * sizeof( double));
      memcpy( val2, val3, 3 * N_OBJECTS * sizeof( double));
      t += delta_t;
      loop_count++;
      if( (long)( clock( ) / CLOCKS_PER_SEC) != seconds_elapsed)
         {
         seconds_elapsed++;
         printf( "%ld seconds elapsed; t = %.1lf\r", seconds_elapsed, t);
         }
      }
   return( 0);
}
