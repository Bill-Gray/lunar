/* jd.cpp: example of use of date/time functions & calendar conversions

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

#define __USE_MINGW_ANSI_STDIO 1
         /* above causes MinGW to use "real" long doubles */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "watdefs.h"
#include "date.h"
#include "afuncs.h"

#ifdef __WATCOMC__
#define floorl floor
#define sinl sin
#endif

#ifdef OBSOLETE_FOR_REFERENCE_ONLY
/* At one point,  I had the months shown by name,  not number.  This is
no longer true,  but I saw no need to chop out the month names for
those who might be interested: */

extern const char *islamic_month_names[12] = {
      "Muharram", "Safar", "Rabi'a I",
      "Rabi'a II", "Jumada I", "Jumada II",
      "Rajab", "Sha'ban", "Ramadan",
      "Shawwal", "Dhu al-Q'adah", "Dhu al-Hijjah" };

extern const char *hebrew_month_names[13] = {
      "Tishri", "Heshvan", "Kislev", "Tevet", "Shevat", "Adar",
      "Adar II", "Nisan", "Iyar", "Sivan", "Tammuz", "Av", "Elul" };

extern const char *french_month_names[12] = {
  "Vendmiaire",           "Germinal",
  "Brumaire",              "Floral",
  "Frimaire",              "Prairial",
  "Nivse",                "Messidor",
  "Pluvise",              "Thermidor",
  "Ventse",               "Fructidor" };

extern const char *french_extra_day_names[6] = {
        "Jour de la vertu (Virtue Day)",
        "Jour du genie (Genius Day)",
        "Jour du travail (Labour Day)",
        "Jour de l'opinion (Reason Day)",
        "Jour des recompenses (Rewards Day)",
        "Jour de la revolution (Revolution Day)" };
#endif


/* The Chinese calendar involves added complications for two reasons.
First,  instead of being computed algorithmically (as all the other
calendars are),  it's computed using a pre-compiled data table.  So you
have to load the file into memory,  as shown in the following code;  the
data is stored using the set_chinese_calendar_data( ) function.

   Second,  as will be seen in the main( ) code,  you have to watch out
for the intercalary months.  */

static int load_chinese_calendar_data( const char *filename)
{
   static char *calendar_data = NULL;
   int rval = 0;

   if( !filename)
      {
      if( calendar_data)
         free( calendar_data);
      calendar_data = NULL;
      }
   else
      {
      FILE *ifile = fopen( filename, "rb");

      rval = -1;
      if( ifile)
         {
         size_t filesize;

         fseek( ifile, 0L, SEEK_END);
         filesize = (size_t)ftell( ifile);
         calendar_data = (char *)malloc( filesize);
         if( calendar_data)
            {
            fseek( ifile, 0L, SEEK_SET);
            if( !fread( calendar_data, filesize, 1, ifile))
               rval = -2;
            else
               {
               set_chinese_calendar_data( calendar_data);
               rval = 0;
               }
            }
         fclose( ifile);
         }
      }
   return( rval);
}

#if !defined( _MSC_VER) && !defined( __WATCOMC__)
static void error_exit( void)  __attribute__ ((noreturn));
#endif

static void error_exit( void)
{
   printf( "jd takes either a Julian Day or a year/month/day as command\n");
   printf( "line arguments.  It then expresses that date in several calendar\n");
   printf( "systems (Julian, Gregorian, Hebrew, Islamic, Persian, Chinese,\n");
   printf( "French Revolutionary).  Example usages:\n\n");
   printf( "jd 2450000.5   (to get calendar data for JD 2450000.5)\n");
   printf( "jd 2007 3 30   (to get calendar data for 2007 March 30)\n");
   printf( "jd 2007-mar-30   (same as above; 30-mar-2007 would also work)\n");
   exit( -1);
}

int main( int argc, char **argv)
{
   int calendar, err_code, i, is_ut;
   const double tdt_minus_tai = 32.184;
   const double tai_minus_gps = 19.;
   const long double j2000 = 2451545.;
   const long double jan_1970 = 2440587.5 - j2000;
               /* Set default initial time to "right now": */
   long double t2k = jan_1970 + (long double)time( NULL) / (long double)seconds_per_day;
   double jd;
   char buff[90];

   if( argc < 2)
      strcpy( buff, "+0");          /* show current time */
   else
      {
      strcpy( buff, argv[1]);
      for( i = 2; i < argc; i++)
         {
         strcat( buff, " ");
         strcat( buff, argv[i]);
         }
      }

   t2k = get_time_from_stringl( t2k, buff,
        CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD | FULL_CTIME_TWO_DIGIT_YEAR, &is_ut);

   if( is_ut < 0)
      printf( "Error parsing string: %d\n", is_ut);
   err_code = load_chinese_calendar_data( "chinese.dat");
   if( err_code)
      printf( "WARNING:  Chinese calendar data not loaded: %d\n", err_code);

   if( t2k + j2000 == 0.)    /* no date found in command line;  show an error message: */
      error_exit( );
   else
      {
      full_ctimel( buff, t2k, CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD
                   | FULL_CTIME_DAY_OF_WEEK_FIRST | FULL_CTIME_12_PLACES);
      printf( "%s = JD %.8Lf\n", buff, t2k + j2000);
      for( calendar = 0; calendar < 9; calendar++)
         {
         int day, month;
         long year;
         const long ljd = (long)floorl( t2k + .5 + j2000);
         static const char *calendar_names[9] = {
                "Gregorian", "Julian", "Hebrew", "Islamic", "Revolutionary",
                "Persian (Jalali)", "Greg/Jul", "Chinese",  "Modern Persian" };
         char is_intercalary = ' ';

         day_to_dmy( ljd, &day, &month, &year, calendar);
         if( calendar == CALENDAR_CHINESE)
            {
            const int chinese_intercalary_month = get_chinese_intercalary_month( );

            if( chinese_intercalary_month)
               {
               if( month == chinese_intercalary_month)
                  is_intercalary = 'i';
               if( month >= chinese_intercalary_month)
                  month--;
               }
            }

         printf( "%-20s %4ld%3d%c%3d\n", calendar_names[calendar],
                        year, month, is_intercalary, day);
         }
      }
   jd = (double)( t2k + j2000);
   printf( "Delta-T = TD - UT1 = %.4f; TD - UTC = %.4f; UT1 - UTC = DUT1 = %.4f\n",
                            td_minus_ut( jd),
                            td_minus_utc( jd),
                            td_minus_utc( jd) - td_minus_ut( jd));
   printf( "TDB - TDT = %f milliseconds   TAI-UTC = %.3f    GPS-UTC = %.3f\n",
         (double)tdb_minus_tdt( t2k / 36525.) * 1000.,
         td_minus_utc( jd) - tdt_minus_tai,
         td_minus_utc( jd) - tdt_minus_tai - tai_minus_gps);
   load_chinese_calendar_data( NULL);
   return( 0);
}
