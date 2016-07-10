/* phases.cpp: creates tables of lunar phase dates/times

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

/* 2013 Jun 30:  'long' integers were being written to the binary
output file,  which should be 32-bit ints. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "watdefs.h"
#include "date.h"
#include "lunar.h"
#include "afuncs.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define VSOP_CHUNK 2767U

#define NEW_MOON 0
#define FIRST_QUARTER 1
#define FULL_MOON 2
#define LAST_QUARTER 3

const double days_per_julian_century = 36525.;

/* Some code to compute lunar phases and produce a table.  It can also   */
/* produce a list of the 'principal terms' used in the Chinese calendar, */
/* i.e.,  the times when the mean solar longitude passes a multiple of   */
/* 30 degrees.  When run with the '-c' (chinese calendar) switch,  it    */
/* produced a huge text file listing new moons and Principal Term dates; */
/* sorted and fed through the 'chinese' program,  one could get the      */
/* actual Chinese calendar data,  which was then stored in 'chinese.dat' */
/* for use in the 'date.cpp' routines for the Chinese calendar.          */

static double approx_solar_dist( double t)
{
   double temp, r[3], rval;

   t /= 10.;      /* cvt to julian millennia */
   temp = 6283.07585 * t;
   r[0] = 100013989. +
            1670700. * cos( 3.0984635 + temp) +
              13956. * cos( 3.05525 + temp * 2) +
               3084. * cos( 5.1985 + 77713.7715 * t);
   r[1] =    103019. * cos( 1.107490 + temp) +
               1721. * cos( 1.0644 + temp * 2.);
   r[2] =      4359. * cos( 5.7846 + temp);
   rval = r[0] + t * (r[1] + t * (r[2] + t));
   rval *= 1.e-8;
   rval *= AU_IN_KM;    /* put return value into KILOMETERS */
   return( rval);
}

int main( const int argc, const char **argv)
{
   const double j2000 = 2451545.;    /* 1.5 Jan 2000 = JD 2451545 */
   int i, julian = 0, verbose = 0, chinese_calendar = 0;
   double t,  t0, utc_time;
   double m, mp, f, e, max_date = 4000. * 365.25 + j2000;
   double lunar, dist, fund[N_FUND], rate = 29.5306, t_final;
   char buff[80];
   double k;
   FILE *log_file = NULL, *vsop_file, *data_file = NULL;
   char *vsop_tbuff, FAR *vsop_data;
   time_t curr_time = time( NULL);

   vsop_file = fopen( "vsop.bin", "rb");
   if( !vsop_file)
      {
      printf( "Couldn't open vsop.bin");
      return( -1);
      }
   vsop_tbuff = (char *)malloc( VSOP_CHUNK);
   vsop_data = (char *)malloc( VSOP_CHUNK * 22U);
   for( i = 0; i < 22; i++)
      {
      if( !fread( vsop_tbuff, VSOP_CHUNK, 1, vsop_file))
         {
         printf( "Couldn't read VSOP data\n");
         free( vsop_tbuff);
         free( vsop_data);
         return( -2);
         }
      FMEMCPY( vsop_data + (unsigned)i * VSOP_CHUNK, vsop_tbuff, VSOP_CHUNK);
      }
   fclose( vsop_file);
   free( vsop_tbuff);

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'j': case 'J':
               julian = 1;
               break;
            case 'c': case 'C':
               chinese_calendar = 1;
               break;
            case 'v': case 'V':
               verbose = 1;
               break;
            case 'l': case 'L':
               log_file = fopen( argv[i] + 2, "wb");
               break;
            case 'd': case 'D':
               data_file = fopen( argv[i] + 2, "wb");
               break;
            case 'm': case 'M':
               max_date = (atof( argv[i] + 2) - 2000.) * 365.25 + j2000;
               break;
            default:
               break;
            }
   t0 = t = j2000 + 365.25 * (atof( argv[1]) - 2000.);
   k = floor( (atof( argv[1]) - 2000.) * 12.3685);
   if( data_file)
      {
      const int32_t int32_t_to_write = (int32_t)k;

      fwrite( &int32_t_to_write, 1, sizeof( int32_t), data_file);
      }
   t_final = max_date - 1.;
   while( t_final < max_date)
      for( i = 0; i < 4; i++)
         {
         double t2, dlon_1, dlon_2, phase_angle, solar_lon, time_lag;
         static const char *phase_name[4] = {
                  "New moon ",
                  "1st qtr. ",
                  "Full moon",
                  "last qtr." };
         double t_centuries, t_cen2, t_cen3, t_cen4;

         t_centuries = k / 1236.85;       /* first approx */
         t = 2451550.09765 + 29.530588853 * k
               + (1.337e-4 - 1.5e-7 * t_centuries) * t_centuries * t_centuries;
         t_centuries = (t - j2000) / 36525.;
         t_cen2 = t_centuries * t_centuries;
         t_cen3 = t_cen2 * t_centuries;
         t_cen4 = t_cen3 * t_centuries;
         m = 2.5534 + 29.10535669 * k
                    -  2.18e-5 * t_cen2
                    -  1.1e-7 * t_cen3;
         mp = 201.5643 + 385.81693528 * k
                       +    .0107438 * t_cen2
                       +   1.239e-5 * t_cen3
                       -   5.8e-8 * t_cen4;
         f = 160.7108 + 390.67050274 * k
                      -   1.6541e-3 * t_cen2
                      -   2.27e-6 * t_cen3
                      +   1.1e-8 * t_cen4;
         m *= PI / 180.;
         f *= PI / 180.;
         mp *= PI / 180.;
         e = 1. - .002516 * t_centuries - 7.4e-6 * t_cen2;
         switch( i)
            {
            case NEW_MOON:
               t +=   -.40720 * sin( mp)
                      +.17241 * sin( m) * e
                      +.01608 * sin( mp + mp)
                      +.01039 * sin( f + f)
                      +.00739 * sin( mp - m) * e
                      -.00514 * sin(  mp + m) * e
                      +.00208 * sin( m + m) * e * e
                      -.00111 * sin( mp - f - f);
               break;
            case FIRST_QUARTER:
            case LAST_QUARTER:
               t +=   -.62801 * sin( mp)
                      +.17172 * sin( m) *  e
                      -.01183 * sin( mp + m) * e
                      +.00862 * sin( mp + mp)
                      +.00804 * sin( f + f)
                      +.00454 * sin( mp - m) * e
                      +.00204 * sin( m + m) * e * e
                      -.00180 * sin( mp - f - f);
               t += ((i == FIRST_QUARTER) ? .00306 : -.00306);
               break;
            case FULL_MOON:
               t +=   -.40614 * sin( mp)
                      +.17302 * sin( m) * e
                      +.01614 * sin( mp + mp)
                      +.01043 * sin( f + f)
                      +.00734 * sin( mp - m) * e
                      -.00515 * sin(  mp + m) * e
                      +.00209 * sin( m + m) * e * e
                      -.00111 * sin( mp - f - f);
               break;
            default:
               break;
            }
         phase_angle = (double)i * 90.;
         t_centuries = (t - j2000) / days_per_julian_century;
         time_lag = approx_solar_dist( t_centuries) / SPEED_OF_LIGHT;
         time_lag /= seconds_per_day * days_per_julian_century;
         lunar_fundamentals( vsop_data, t_centuries, fund);
         lunar_lon_and_dist( vsop_data, fund, &lunar, &dist, 0L);
         solar_lon = calc_vsop_loc( vsop_data, 3, 0, t_centuries - time_lag, 0.);
         solar_lon = solar_lon * 180. / PI - 180.;
         dlon_1 = lunar - solar_lon - phase_angle;
         while( dlon_1 < -180.) dlon_1 += 360.;
         while( dlon_1 >  180.) dlon_1 -= 360.;
         if( verbose)
            {
            full_ctime( buff, t, julian);
            printf( "   first  time is %s;  lon diff %f\n", buff, dlon_1);
            }

         t2 = t - rate * dlon_1 / 360.;
         t_centuries = (t2 - j2000) / 36525.;
         lunar_fundamentals( vsop_data, t_centuries, fund);
         lunar_lon_and_dist( vsop_data, fund, &lunar, &dist, 0L);
         solar_lon = calc_vsop_loc( vsop_data, 3, 0, t_centuries - time_lag, 0.);
         solar_lon = solar_lon * 180. / PI - 180.;
         dlon_2 = lunar - solar_lon - phase_angle;
         while( dlon_2 < -180.) dlon_2 += 360.;
         while( dlon_2 >  180.) dlon_2 -= 360.;
         if( verbose)
            {
            full_ctime( buff, t2, julian);
            printf( "   second time is %s;  lon diff %f\n", buff, dlon_2);
            }
         t_final = (t * dlon_2 - t2 * dlon_1) / (dlon_2 - dlon_1);
         utc_time = t_final - td_minus_utc( t_final) / seconds_per_day;
         full_ctime( buff, utc_time, julian);
         if( log_file)
            {
            if( chinese_calendar)
               fprintf( log_file, "%7ld    %s\n",
                                (long)floor( utc_time - 1. / 6.), buff);
            else
               fprintf( log_file, "%s: %s\n", phase_name[i], buff);
            }
         buff[17] = '\0';        /* trim the seconds */
         if( !chinese_calendar)
            {
            printf( "%s  ", buff);
            if( i == 3)
               printf( "\n");
            }
         if( data_file)
            {
            double diff = t_final - ( 2451550.09765 + 29.530588853 * k);
            int32_t int32_t_to_write;

            diff *= seconds_per_day;
            if( verbose)
               printf( "Time diff %f\n", diff);
            int32_t_to_write = (int32_t)diff;
            fwrite( &int32_t_to_write, 1, sizeof( int32_t), data_file);
            }
         if( chinese_calendar)
            {
            k++;
            i = 4;
            }
         else        /* just moving to next _phase_ */
            k += .25;
         if( curr_time != time( NULL))
            {
            curr_time = time( NULL);
            printf( "JD %.3f (%.3f)\r", t_final,
                     (t_final - j2000) / 365.25 + 2000.);
            }
         }

   if( chinese_calendar)
      while( t0 < max_date)
         {
         double t_centuries, time_lag, delta_t = 1.;
         double solar_lon;
         const double thirty_deg = PI / 6.;
         long solar_month = 0L, year;

         while( fabs( delta_t) > .00001)  /* resolution a little better than 1s */
            {
            t_centuries = (t0 - j2000) / 36525.;
            time_lag = approx_solar_dist( t_centuries) / SPEED_OF_LIGHT;
            time_lag /= seconds_per_day * days_per_julian_century;
            solar_lon = calc_vsop_loc( vsop_data, 3, 0, t_centuries - time_lag, 0.);
            solar_month = (long)floor( solar_lon / thirty_deg + .5);
            solar_lon -= (double)solar_month * thirty_deg;
            delta_t = solar_lon * 365.25 / (2. * PI);
            t0 -= delta_t;
            if( verbose)
               printf( "delta_t: %.5f   JD %.5f\n", delta_t, t0);
            }
         solar_month = (solar_month + 7) % 12;  /* flip around to winter sol */
         year = (long)
                  floor( (t0 - (double)solar_month * 30.5 - 100.) / 365.25);
         year -= 2074L;
         utc_time = t0 - td_minus_utc( t0) / seconds_per_day;
         full_ctime( buff, utc_time, julian);
         if( log_file)
            fprintf( log_file, "%7ld z  %s  %2ld %5ld\n",
                          (long)floor( utc_time - 1. / 6.),
                          buff, solar_month + 1, year);
         if( chinese_calendar && curr_time != time( NULL))
            {
            curr_time = time( NULL);
            printf( "%7ld z  %s   %2ld %5ld: %.3f\n",
                          (long)floor( utc_time - 1. / 6.),
                          buff, solar_month + 1, year,
                          (utc_time - j2000) / 365.25 + 2000.);
            }
         t0 += 365.25 / 12.;
         }
   if( data_file)
      fclose( data_file);
   if( log_file)
      fclose( log_file);
   return( 0);
}
