/* persian.cpp: test code for computing the Persian (Jalali) calendar

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

#include <math.h>
#include "watdefs.h"
#include "afuncs.h"

/*
See this code also in DATE.CPP.  It determines the date of the beginning
of the Persian year,  in JD form.  It so happens that this can be expressed,
in somewhat "condensed" form,  using the following not-too-terribly-complex
algorithm.  The main( ) code provides a test routine comparing the results
of this algorithm to those derived the "right" way (computing the time of
the March equinox,  subtracting Delta-T,  adding the longitude of Teheran.)
*/

#define JALALI_ZERO 1947954L
#define LOWER_PERSIAN_YEAR -1096

static long jalali_jd0( const int jalali_year)
{
   static const short breaks[12] = { -708, -221,   -3,    6,  394,  720,
                                      786, 1145, 1635, 1701, 1866, 2328 };
   static const short deltas[12] = { 1108, 1047,  984, 1249,  952,  891,
                                      930,  866,  869,  844,  848,  852 };
   int i;
   long rval;

   if( jalali_year < LOWER_PERSIAN_YEAR)
      return( -1L);           /* out of valid range */
   for( i = 0; i < 12; i++)
      if( jalali_year < breaks[i])
         {
         rval = JALALI_ZERO + (long)jalali_year * 365L +
                     (long)( deltas[i] + jalali_year * 303L) / 1250L;
         if( i < 3)  /* zero point drops one day in first three blocks */
            rval--;
         return( rval);
         }
   return( -1L);           /* out of valid range */
}

#include <stdio.h>
#include <stdlib.h>

/*
   I found a value for the longitude of Teheran of E 51 26'.  I'm not
absolutely sure,  though,  that this coincides with the longitude of
the observatory used for the Persian calendar.  If it's off by a few
arcminutes,  then it is conceivable that,  in situations where the
equinox is almost exactly at noon Teheran time,  the calendar would be
off by one day during that year. */

#define TEHERAN_LONGITUDE ((51. + 26. / 60.) / 360.)

/* The longitude of Paris is that given for MPC station 007: */

#define PARIS_LONGITUDE (2.3371 / 360.)

static int revolutionary = 1;

double get_solstice_equinox_date( double jd);

static double jd_equinox( long year)
{
   double jd0, y = (double)( year - 2000L) / 1000.;

   if( revolutionary)      /* we're after the autumnal equinox */
      jd0 = 2451810.217 + y * 365242.01767;
   else                    /* after the spring equinox */
      jd0 = 2451623.809 + y * 365242.37404;

   return( get_solstice_equinox_date( jd0));
}

int main( int argc, const char **argv)
{
   long year = atol( argv[1]), n_years, year0;
   int i, j;
   double longitude;

   for( i = 3; i < argc; i++)
      {
      if( argv[i][0] == '-')
         {
         switch( argv[i][1])
            {
            case 'j':
               revolutionary = 0;
               printf( "Switching to Persian (Jalaali) calendar\n");
               break;
            default:
               printf( "Unknown option %s\n", argv[i]);
               exit( -1);
            }
         for( j = i; j < argc - 1; j++)      /* remove the switch */
            argv[j] = argv[j + 1];
         i--;
         argc--;
         }
      }
   longitude = (revolutionary ? PARIS_LONGITUDE : TEHERAN_LONGITUDE);
   year0 = (revolutionary ? 1791L : 621L);
   if( argc == 2)
      printf( "%.5f\n", jd_equinox( year));
   if( argc == 3)
      {
      n_years = atol( argv[2]);

      while( n_years--)
         {
         long greg_year = year + year0;
         double jd1 = jd_equinox( greg_year     ) + longitude;
         double jd2 = jd_equinox( greg_year + 1L) + longitude;
         double delta_t = td_minus_ut( 2451545. +
                               (double)(greg_year - 2000L) * 365.25);
         long remains, n_days, is_leap;

               /* Subtract delta_t,  to convert the ET (Ephemeris Time) */
               /* value given by jd_equinox( ) into a UT value.   */

         jd1 -= delta_t / seconds_per_day;
         jd2 -= delta_t / seconds_per_day;

         n_days = (long)( floor( jd2) - floor( jd1));
         printf( "%ld, %ld: %ld   ", greg_year, year, n_days);
                        /* 'Is_leap' is determined using the 33-year */
                        /* cycle used in some other programs.        */
         remains = year % 33;
         if( remains > 18)
            remains--;
         is_leap = (remains % 4L == 1L);
         printf( is_leap ? "Is leap\n" : "Is normal\n");
                        /* Check to see if the 33-year cycle method  */
                        /* matches the 'accurate' method:            */
         if( n_days - is_leap != 365L)
            printf( "DANGER! ALERT!\n");
                        /* Now check to see if my own algorithm matches */
                        /* the 'accurate' method (it does,  over the    */
                        /* years in question):                          */
         if( (long)jd1 != jalali_jd0( (int)year))
            printf( "PANIC: %f, %ld\n", jd1, jalali_jd0( (int)year));
         year++;
         }
      }

                        /* The following code was used in creating the */
                        /* list of 'breaks' and 'deltas' used in the   */
                        /* jalali_jd0( ) function.  It computes the JD */
                        /* of the start of the Persian calendar over a */
                        /* given range of years;  then it tries to find */
                        /* 'delta' value(s) for which the simplified    */
                        /* method of the jalali_jd0( ) returns correct  */
                        /* values.  It's of basically historical        */
                        /* interest now.                                */
                        /* e.g: 'persian 1145 1635 303 1250' verifies   */
                        /* a date range used in the jalali_jd0()        */
                        /* function in date.cpp (and shown above)       */
   if( argc >= 4)
      {
      long y1 = atol( argv[1]);
      long y2 = atol( argv[2]);
      long *jds = (long *)calloc( (size_t)(y2 - y1), sizeof( long));
      long offset, max_year = -100000;
      long const1 = atol( argv[3]), const2 = 10000L;

      if( !jds)
         {
         printf( "Out of memory\n");
         exit( 0);
         }
      if( !const1)
         {
         const1 = 683L;       /* usual rule */
         const2 = 2820L;
         }
      else if( argc == 5)
         const2 = atol( argv[4]);
      for( year = y1; year < y2; year++)
         {
         long greg_year = year + year0;
         double jd1 = jd_equinox( greg_year) + longitude;

         double delta_t = td_minus_ut( 2451545. +
                               (double)(greg_year - 2000L) * 365.25);

         jd1 -= delta_t / seconds_per_day;
         if( revolutionary)         /* Jalaali calendar is noon-based,  */
            jd1 += .5;              /* Revolutionary is midnight-based  */
         jds[year - y1] = (long)floor( jd1);
         }

      for( offset = 0; offset < const2; offset++)
         {
         long jd0 = jds[0] - y1 * 365L -
                            ( y1 * const1 + offset) / const2;
         int mismatch = 0;

         for( year = y1; !mismatch && year < y2; year++)
            {
            long jd_algo = jd0 + year * 365L +
                            ( year * const1 + offset) / const2;

            if( jd_algo != jds[year - y1])
               {
               mismatch = 1;
               if( year > max_year)
                  max_year = year;
               }
            }
         if( !mismatch)
            printf( "Offset: %ld JD0: %ld\n", offset, jd0);
         }
      printf( "Max year reached: %ld\n", max_year);
      free( jds);
      }
   return( 0);
}
