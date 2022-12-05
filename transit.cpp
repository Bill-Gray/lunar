/* Copyright (C) 2019, Project Pluto

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
02110-1301, USA.

   This is a small utility actually written about 25 years ago
to compute times of transits of Mercury and Venus as seen from
Earth (or for any planet as seen from one outside its orbit),
or to compute times of greatest elongation.  It works,  but was
written when I was a younger and more foolish programmer.  I
wrote it to generate tables which were then shown in my Guide
desktop planetarium software,  so making it user-friendly was
not a goal;  it only had to work once.  Hence the lack of
comments and strange 'logic' in parts of the code.     */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <math.h>
#include "watdefs.h"
#include "date.h"
#include "lunar.h"

#define GREGORIAN 0
#define VSOP_SIZE 60874U
#define PI 3.14159265358979323846
#define SOLAR_RADIUS (959.63 * (PI / 180.) / 3600.)

static double planet_radius[10] = { 696265., 2439., 6051., 6378.,
                      3393., 71398., 60000., 25400., 25295., 1500. };

static double radius[10] = { 0., .39, .72, 1., 1.52, 5.2, 9.54, 19., 30., 31. };

static double t0[10] = { 0., 88., 224.7, 365.25, 687., 4332.1, 10825.9,
                        30676.1, 59911.1, 90824.2 };
#ifdef UNUSED
static double r0[10] = { 0., .25, .72, 1., 1.52, 5.2, 9.54, 19.2, 30.0, 40.0 };
#endif

void lon_lat_diff( char *buff, double t, int planet1, int planet2,
                  double *d_lon, double *d_lat, double *elong)
{
   double lat1, lat2, lon1, lon2;

   const double t_c = (t - 2451545.) / 36525.;     /* cvt to julian centuries */
   const double r1 = calc_vsop_loc( buff, planet1, 2, t_c, 0.);
   const double r2 = calc_vsop_loc( buff, planet2, 2, t_c, 0.);
   const double t_c2 = t_c - (r1 - r2) / (AU_PER_DAY * 36525.);

   assert( d_lat || d_lon);
   if( d_lat)
      {
      lat1 = calc_vsop_loc( buff, planet1, 1, t_c, 0.);
      lat2 = calc_vsop_loc( buff, planet2, 1, t_c2, 0.);
      *d_lat = lat1 - lat2;
      }
   if( d_lon)
      {
      lon1 = calc_vsop_loc( buff, planet1, 0, t_c, 0.);
      lon2 = calc_vsop_loc( buff, planet2, 0, t_c2, 0.);
      *d_lon = lon1 - lon2;
      if( *d_lon < -PI) *d_lon += 2. * PI;
      if( *d_lon >  PI) *d_lon -= 2. * PI;
      }
   if( elong && d_lat && d_lon)
      {
      const double dx = r1 * cos( lat1) * cos( lon1) - r2 * cos( lat2) * cos( lon2);
      const double dy = r1 * cos( lat1) * sin( lon1) - r2 * cos( lat2) * sin( lon2);
      const double dz = r1 * sin( lat1)              - r2 * sin( lat2);
      const double dist = sqrt( dx * dx + dy * dy + dz * dz);
      const double cos_elong = ( r1 * r1 + dist * dist - r2 * r2) / (2. * r1 * dist);

      *elong = acos( cos_elong);
      }
}

int main( const int argc, const char **argv)
{
   int planet1 = atoi( argv[1]), planet2 = atoi( argv[2]);
   int day = 0, month = 1, i, verbose = 0, getting_elong = 0;
   int for_etb = 0;
   long year = atol( argv[3]), loop_count = 0;
   double max_t = 3e+6;
   double synodic, t, d_lon1, d_lon2, d_lon, d_lat1, check_range;
   char *buff = (char *)malloc( VSOP_SIZE);
   FILE *ifile = fopen( "vsop.bin", "rb");
   FILE *logfile = NULL;
   char str[80];

   assert( ifile);
   assert( buff);
         /* assume planet1 = outer planet,  planet2 = inner */
   if( planet1 < planet2)
      {
      int temp = planet1;
      planet1 = planet2;
      planet2 = temp;
      }

   for( i = 0; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'v': case 'V':
               verbose = 1;
               break;
            case 'l': case 'L':
               logfile = fopen( argv[i] + 2, "wb");
               break;
            case 'e': case 'E':
               getting_elong = 1;
               break;
            case 't': case 'T':
               for_etb = 1;
               break;
            case 'm': case 'M':
               max_t = (atof( argv[i] + 2) - 2000.) * 365.25 + 2451545;
               break;
            default:
               break;
            }

   synodic = 1. / (1. / t0[planet2] - 1. / t0[planet1]);
   if( verbose)
      printf( "Synodic period is %lf days\n", synodic);

   if( fread( buff, 1, VSOP_SIZE, ifile) != VSOP_SIZE)
      {
      fprintf( stderr, "vsop.bin not properly read\n");
      return( -1);
      }
   fclose( ifile);

   for( i = 0; i < 10; i++)
      planet_radius[i] /= AU_IN_KM;    /* cvt planet radii to AU */

   check_range = (planet_radius[0] + planet_radius[planet2]) / radius[planet2]
              - (planet_radius[0] - planet_radius[planet1]) / radius[planet1];
   printf( "Default check range is %lf degrees\n", check_range * 180. / PI);
   check_range *= 1.5;      /* safety margin */

            /* calc longitude difference for first approx */
   t = dmy_to_day( day, month, year, GREGORIAN);
   for( i = 0; i < 3; i++)
      {
      lon_lat_diff( buff, t, planet1, planet2, &d_lon, NULL, NULL);
                  /* make sure we BACK UP first.... */
      if( !i && d_lon > 0.)
         d_lon -= 2. * PI;
      t += synodic * d_lon / (2. * PI);
      full_ctime( str, t, 0);
      if( verbose)
         printf( "approx %d: %s\n", i + 1, str);
      }

   while( t < max_t)
      {
      double corr1, corr2, dt;

      if( !getting_elong || !(loop_count & 1))
         {
         lon_lat_diff( buff, t, planet1, planet2, &corr1, NULL, NULL);
         if( getting_elong && (loop_count & 2))    /* superior conj */
            {
            corr1 += PI;
            while( corr1 > PI)
               corr1 -= PI + PI;
            }
         dt = corr1 * synodic / (2. * PI);
         t += dt;
         lon_lat_diff( buff, t, planet1, planet2, &corr2, NULL, NULL);
         if( getting_elong && (loop_count & 2))    /* superior conj */
            {
            corr2 += PI;
            while( corr2 > PI)
               corr2 -= PI + PI;
            }
         t += corr2 * dt / (corr1 - corr2);
         lon_lat_diff( buff, t, planet1, planet2, NULL, &d_lat1, NULL);
         }

      if( getting_elong)
         {
         double stepsize, elong[3], delta = 9999.;
         const char *text[4] = { "Inf Conj", "West Elo",
                                 "Sup Conj", "East Elo" };
         char obuff[120];

         while( delta > .001 || delta < -.001)
            {
            stepsize = delta;
            if( fabs( stepsize) > 1.)
               stepsize = 1.;
            for( i = 0; i < 3; i++)
               lon_lat_diff( buff, t + (double)(i-1) * stepsize, planet1,
                              planet2, &corr1, &corr2, elong + i);
            delta = stepsize * .5 * (elong[0] - elong[2])
                                / (elong[0] + elong[2] - 2. * elong[1]);
            if( verbose)
               printf( "Shift: %lf\n", delta);
            t += delta;
            }
         full_ctime( str, t, 0);
         assert( strlen( str) < 40);
         assert( elong[1] > 0. && elong[1] < PI);
         snprintf( obuff, sizeof( obuff), "%s %s %5.2lf", text[loop_count & 3], str,
                                              elong[1] * 180. / PI);
         printf( "%s\n", obuff);
         if( logfile && for_etb)
            fprintf( logfile, "^%s//tdj%.4lf;gp%d^\n", obuff, t, planet2);
         }

      if( !getting_elong && fabs( d_lat1) < check_range)
         {                /* close enough to merit closer look */
         double d_lat2, t_c, time_diff;
         double a, b, c;         /* quadratic coeffs */
         double dx, dy, discrim;
         double r1, r2, real_range;

         lon_lat_diff( buff, t, planet1, planet2, &d_lon1, NULL, NULL);
         time_diff = synodic * d_lon1 / (2. * PI);
         lon_lat_diff( buff, t + time_diff, planet1, planet2, &d_lon2, &d_lat2, NULL);
               /* now calculate closest approach time and dist... */

         dx = d_lon2 - d_lon1;
         dy = d_lat2 - d_lat1;
         a = dx * dx + dy * dy;              /* quadratic for intercept */
         b = 2. * (dx * d_lon1 + dy * d_lat1);
         t -= b * time_diff / (2. * a);

         t_c = (t - 2451545.) / 36525.;     /* cvt to julian centuries */
         r1 = calc_vsop_loc( buff, planet1, 2, t_c, 0.);
         r2 = calc_vsop_loc( buff, planet2, 2, t_c, 0.);
         real_range = (planet_radius[0] + planet_radius[planet2]) / r2
                 - (planet_radius[0] - planet_radius[planet1]) / r1;
         full_ctime( str, t, 0);
         if( verbose)
            printf( "Conjunction time: %s; %.3lf\r", str, real_range * 180. / PI);
         c = d_lon1 * d_lon1 + d_lat1 * d_lat1 - real_range * real_range;
         discrim = b * b - 4. * a * c;
         if( discrim > 0.)
            {
            time_diff = fabs( time_diff * sqrt( discrim) * .5 / a);
            full_ctime( str, t - time_diff, 0);
            printf( "TRANSIT starts at %s        \n", str);
            if( logfile && !for_etb)
               fprintf( logfile, "TRANSIT starts at %s        \n", str);
            full_ctime( str, t + time_diff, 0);
            printf( "TRANSIT ends at   %s        \n\n", str);
            if( logfile && !for_etb)
               fprintf( logfile, "TRANSIT ends at   %s        \n\n", str);
            if( logfile && for_etb)
               {
               double lat_min = d_lat1 - dy * b / (2. * a);
               double lon_min = d_lon1 - dx * b / (2. * a);
               double sep_min = sqrt( lat_min * lat_min + lon_min * lon_min);

               if( lat_min > 0.)
                  sep_min *= -1.;
               sep_min *= r2 / (r1 - r2);
               full_ctime( str, t, 0);
               fprintf( logfile, "^%s//tdj%12.4lf;gp%d^ %4.1lf %5.2lf\n", str,
                   t, planet2, 2. * time_diff * 24., sep_min * 180. / PI);
               }
            }
         }
      if( !getting_elong)
         t += synodic;
      else
         {
         double conj_to_elong = (planet2 == 1 ? +21.6 : 70.8);
         switch( loop_count & 3)
            {
            case 0:        /* just computed inferior conjunction */
            case 3:        /* just computed evening elong        */
               t += conj_to_elong;
               break;
            case 2:        /* just computed superior conjunction */
            case 1:        /* just computed morning visibility */
               t += synodic / 2. - conj_to_elong;
               break;
            }
         }
      loop_count++;
      }
   if( logfile)
      fclose( logfile);
   free( buff);
   return( 0);
}
