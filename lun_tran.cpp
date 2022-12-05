/* lun_tran.cpp: computes lunar transit times

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
#include "watdefs.h"
#include "lunar.h"
#include "date.h"
#include "afuncs.h"
#include "riseset3.h"
#include "stringex.h"

const static double pi =
     3.1415926535897932384626433832795028841971693993751058209749445923078;

double get_lunar_transit_time( const int year, const int month,
   const int day, const double latitude, const double longitude,
   const int time_zone, const int dst, const int real_transit);
int get_zip_code_data( const int zip_code, double *latitude,
      double *longitude, int *time_zone, int *use_dst,
      char *place_name);                              /* lun_tran.cpp */
void format_hh_mm( char *buff, const double time);    /* lun_tran.cpp */

static double look_for_transit_time( const int planet_no,
                  double jd,
                  const double observer_lat, const double observer_lon,
                  const char *vsop_data, const int real_transit)
{
   int iterations = 10;
   double delta = 1.;

   do
      {
      PLANET_DATA pdata;

      fill_planet_data( &pdata, planet_no, jd,
                          observer_lat, observer_lon, vsop_data);

      delta = atan2( pdata.altaz_loc[1], pdata.altaz_loc[0]);
      if( !real_transit)
         delta += pi;
      while( delta < -pi)
         delta += pi + pi;
      while( delta > pi)
         delta -= pi + pi;
      if( planet_no == 10)       /* moon moves a little faster */
         delta *= 30.5 / 29.5;
      delta *= sqrt( 1. - pdata.altaz_loc[2] * pdata.altaz_loc[2]);
      jd += delta / (2. * pi);
      }
      while( fabs( delta) > .0001 && iterations--);
   return( jd);
}

/* 2006 and before:  US Daylight "Saving" Time changes are on the
first Sunday in April,  last Sun in October.
   2007 and later:  second Sunday in March,  first Sunday in November. */

static int is_dst( const int minutes, const int day_of_year, const int year)
{
   int rval = 0;
   const int day1 = year + year / 4 - year / 100 + year / 400;
   const int is_leap = !(year % 4) && (year % 100 || !(year % 400));
   int start_day,  end_day;

   if( year <= 2006)
      {
      start_day = 31 + 28 + is_leap + 31 + 7 - (day1 + 5) % 7 - 1;
      end_day = 31 + 28 + is_leap + 31 + 30 + 31 + 30 + 31 + 31 + 30 +
                      31 - (day1 + 2) % 7 - 1;
      }
   else
      {
      start_day = 31 + 28 + is_leap + 14 - (day1 + 2) % 7 - 1;
      end_day = 31 + 28 + is_leap + 31 + 30 + 31 + 30 + 31 + 31 + 30 +
                      31 + 7 - (day1 + 2) % 7 - 1;
      }

   if( day_of_year > start_day && day_of_year < end_day)
      rval = 1;
   if( day_of_year == start_day && minutes >= 120)        /* after 2 PM */
      rval = 1;
   if( day_of_year == end_day && minutes < 120)        /* before 2 AM */
      rval = 1;
   return( rval);
}

/* Finds the time of transit (if real_transit == 1) or of "anti-transit",
   when the moon is lowest in the sky (if real_transit == 0).  For example,

   get_lunar_transit_time( 2011, 7, 29, 44.01, -69.9, -5, 1, 1);

   would return the time of the lunar transit for 2011 July 29,  as seen
   from latitude +44.01, longitude -69.9 (Bowdoinham, Maine),  in zone -5
   (Eastern US Time),  with daylight "saving" time included,  and looking
   for the actual transit.

   Return values:

      -2.         The 'vsop.bin' file wasn't found.
      <0. or >=1. The event occurred slightly before or slightly after
                  that calendar day.  Lunar transits average about 25 hours
                  apart,  so this will happen about once a month for each
                  type of event.
      >=0 and <1. Fraction of a day when the event occurred.  Thus, 0=start
                  of day at midnight,  .25 = 6:00 AM,  .7 = 4:48 PM.  See
                  the 'format_hh_mm( )' function below.

      Limitations:
                  Daylight "Saving" Time throws sand in the works.  On the
                  change day in autumn,  the day is 25 hours long,  such that
                  you can have two transits in one day;  this code will find
                  only one of them.  Also,  in autumn,  times such as 2:33 are
                  not determinate;  they could be "2:33 before clocks were set
                  back" or "2:33 after clocks were set back".  This is very
                  rarely a problem,  but one should be aware of it.
*/

double get_lunar_transit_time( const int year, const int month,
   const int day, const double latitude, const double longitude,
   const int time_zone, const int dst, const int real_transit)
{
   const long jd0 = dmy_to_day( day, month, year, 0);
   const long day_of_year = jd0 - dmy_to_day( 1, 1, year, 0);
   double jd = (double)jd0;
   static char *vsop_data = NULL;
   int dst_hours = (dst ? is_dst( 720, (int)day_of_year, year) : 0);

   if( !vsop_data)
      vsop_data = load_file_into_memory( "vsop.bin", NULL);
   if( !vsop_data)
      return( -2.);
   jd -= (double)time_zone / 24.;      /* convert local to UT */
   jd -= dst_hours / 24.;

   jd = look_for_transit_time( 10, jd,
            latitude * pi / 180., longitude * pi / 180.,
            vsop_data, real_transit);
   jd += (double)time_zone / 24.;      /* convert UT back to local */
   jd += dst_hours / 24.;
   return( jd - (double)jd0 + .5);
}

/* Given a "time" from 0 to 1,  referring to a fraction of a day, */
/* format_hh_mm( ) produces a time string in hours and minutes.   */
/* If the time is out of range -- as happens routinely with the   */
/* transit times -- the time is shown as "--:--".                 */

void format_hh_mm( char *buff, const double time)
{
   if( time < 0. || time >= 1.)
      strcpy( buff, "--:--");
   else
      {
      const int minutes = (int)( time * 24. * 60.);

      snprintf_err( buff, 6, "%02d:%02d", minutes / 60, minutes % 60);
      }
}

int get_zip_code_data( const int zip_code, double *latitude,
      double *longitude, int *time_zone, int *use_dst,
      char *place_name)
{
   FILE *ifile = fopen( "zips5.txt", "rb");
   char buff[100];
   int rval = -1;

   if( !ifile)
      return( -2);
   while( rval && fgets( buff, sizeof( buff), ifile))
      if( atoi( buff) == zip_code)
         {
         *longitude = atof( buff + 22);
         *latitude = atof( buff + 11);
         *time_zone = atoi( buff + 35);
         *use_dst = atoi( buff + 39);
         if( place_name)
            {
            int i;

            for( i = 0; buff[i + 48] >= ' '; i++)
               place_name[i] = buff[i + 48];
            buff[8] = '\0';
                  /* Add on two-letter state abbr: */
            snprintf_err( place_name + i, 5, ", %s", buff + 6);
            }
         rval = 0;
         }
   fclose( ifile);
   return( rval);
}
