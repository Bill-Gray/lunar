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

/* This is code written to create the pre-generated table required for
Chinese calendar conversion in 'date.cpp' (see the 'get_chinese_year_data'
function).  The actual "astronomy" bit (computing times of New Moons
and principal terms) takes place in 'phases.cpp';  that program
emits a list of times of those events,  along with the (integer)
JD of the time after adjusting for longitude 120 East.

   One then sorts that log file and feeds it to this program,
which crunches it into a binary form in which each year requires
only three bytes of data.  The bit-packing may be considered extreme,
but was significant when this code was written circa 1995.  It
meant that a table for ten millennia could be fitted into a 30K
file (see 'chinese.dat',  which is that file for Gregorian years
-3000 to +7000).

   Also see https://projectpluto.com/calendar.htm#chinese for a
better explanation of how the Chinese calendar and this system work.

   Note also that this is an observational calendar,  and is susceptible
to uncertainties in Delta-T.  I computed 'chinese.dat' using the best
available TD-UT formula in 1995;  running again in 2020,  I see the
first discrepancy in the future at (Gregorian) 2173,  with the number
of discrepancies gradually increasing as you go forward.  Going
backward,  I see discrepancies in 1570.  (With the added warning that
I seriously doubt the Chinese calendar was measured precisely enough
back then to allow my computation to match the historical record.
Though actually,  it'd probably be good for any month where New Moon
didn't occur really close to the start/end of the day.)

   To make 'chinese.dat',  I ran (and you can run)

./phases -3000 -m8000 -c -l/tmp/log.txt

   I edited /tmp/log.txt and sorted it in straight ASCII order,
causing principal term data and new moon times to become interleaved,
and saved the result.  Then :

./chinese /tmp/log.txt 10000

   generated the actual 'chinese.dat' file.        */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "stringex.h"

#define YEAR_DATA struct year_data

#pragma pack(1)

YEAR_DATA
   {
   uint16_t mask;
   uint8_t intercalary_month, offset;
   } *ydata;

#pragma pack()

int16_t n_years, year0;
uint16_t bitmask = 1u;
long prev_jd;

static void dump_data( char *text[], const int n_found)
{
   int i, month = 12, year = atoi( text[0] + 35), intercalary_month = 0;

   if( text[0][8] != 'z')
      printf( "NOT SOLSTICE\n");
   if( n_found < 24 || n_found > 25)
      {
      printf( "n_found = %d\n", n_found);
      exit( 0);
      }
   if( !year0)               /* Remember,  we get the last month of the */
      year0 = (short)( year + 1);    /* previous year... skip it */

   printf( "\n");
   for( i = 0; i < n_found; i++)
      {
      if( text[i][8] == ' ')        /* it's a month */
         {
         const long jd = atol( text[i]);
         const int prev_month = (month + 10) % 12 + 1;
         const int prev_year = year - (prev_month == 12);

         if( jd - prev_jd == 30L)
            if( prev_year >= year0 && prev_year < year0 + n_years)
               ydata[prev_year - year0].mask |= bitmask;
         bitmask <<= 1;

         if( month == 1 && intercalary_month != 1)
            {           /* it's the for-real first month of the year: */
            const long offset = jd - (year * 365L + year / 4L + 757861L);

            bitmask = 1;
            if( year >= year0 && year < year0 + n_years)
               ydata[year - year0].offset = (unsigned char)offset;
            }

         if( n_found == 25 && !intercalary_month && text[i + 1][8] == ' ')
            {
            snprintf_err( text[i] + 31, 12, "%4di%5d\n", prev_month, prev_year);
            intercalary_month = prev_month;
            if( prev_year >= year0 && prev_year < year0 + n_years)
               ydata[prev_year - year0].intercalary_month =
                                                 (unsigned char)prev_month;
            }
         else           /* it was just a normal month... move on by one */
            {
            snprintf_err( text[i] + 31, 12, "%4d%6d\n", month, year);
            month++;
            if( month == 13)
               {
               month = 1;
               year++;
               }
            }

         prev_jd = jd;
         }
      printf( "%s", text[i]);
      }
}

int main( const int argc, const char **argv)
{
   char *buff, *text[30], ibuff[60];
   int i, n_found = 0;
   FILE *ifile = fopen( argv[1], "rb");

   if( !ifile)
      {
      printf( "%s not opened\n", argv[1]);
      return( -1);
      }
   if( argc == 3)
      {
      n_years = (short)atoi( argv[2]);
      ydata = (YEAR_DATA *)calloc( n_years, sizeof( YEAR_DATA));
      printf( "Looking for %d years of data\n", n_years);
      }
   buff = (char *)malloc( 30 * 60);
   for( i = 0; i < 30; i++)
      text[i] = buff + i * 60;
   while( fgets( ibuff, 60, ifile))
      {
      if( ibuff[26] == ':')      /* workaround for years before -999 */
         memmove( ibuff + 17, ibuff + 18, strlen( ibuff + 17));
      strcpy( text[n_found], ibuff);
      if( ibuff[33] == '1' && ibuff[34] == '1')
         {
         if( n_found == 24 || n_found == 25)
            dump_data( text, n_found);
         n_found = 0;
         }
      strcpy( text[n_found++], ibuff);
      memset( ibuff, 0, 60);
      if( n_found > 25)
         {
         printf( "ERROR in input file\n");
         for( i = 0; i < n_found; i++)
            printf( "%s", text[i]);
         exit( 0);
         }
      }
   if( ydata)
      {
      FILE *ofile = fopen( "chinese.dat", "wb");

      printf( "Writing %d years,  year0 = %d\n", n_years, year0);
      fwrite( &n_years, sizeof( short), 1, ofile);
      fwrite( &year0, sizeof( short), 1, ofile);
      for( i = 0; i < n_years; i++)
         {
         unsigned long oval = (unsigned long)ydata[i].mask +
                (unsigned long)ydata[i].intercalary_month * 8192L +
                (unsigned long)ydata[i].offset * 8192L * 14L;
         fwrite( &oval, 3, 1, ofile);

         printf( "  Year %4d: offset %d, intercalary %2d, mask %04x\n",
            i + year0,
            ydata[i].offset,
            ydata[i].intercalary_month,
            ydata[i].mask);
         }
      fclose( ofile);
      }
   return( 0);
}
