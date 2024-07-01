/* easter.cpp: functions for computing date of Easter, plus test code

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

/* See Meeus,  _Astronomical Algorithms_,  p 67.  He in turn states that
   "the following method has been given by Spencer Jones in his book
   _General Astronomy_ (p 73-74 of the 1922 edition).  It has been
   published again in the _Journal of the British Astronomical Association_,
   vol 88, page 91 (December 1977),  where it is said that it was devised
   in 1876 and appeared in Butcher's _Astronomical Calendar._

      Unlike the formula given by Gauss,  this method has no exception and
   is valid for all years in the Gregorian calendar,  hence from 1583 on."

   Note that this formula results from the Gregorian calendar change.
Before 1583,  different rules were used to determine the date of Easter.
You can use the formula for dates before 1583,  but the results will
not match the dates when Easter was actually observed.

   The relevant pages from _General Astronomy_ can be found at

https://archive.org/details/generalastronomy0000hspe/page/72/mode/2up

   The test program can be run with zero,  one,  or two command line
arguments,  to get three very different types of output.

   Easter repeats over a cycle of 5.7 million years.  If run without
command-line arguments,  this program computes the day of Easter over an
entire cycle.  That lets it show the frequency on which Easter occurs on
a given day;  the result looks like

Mar 22: 0.48333 = 29/60    Apr  3: 3.38333 = 203/60   Apr 15: 3.38333 = 203/60
Mar 23: 0.95000 = 19/20    Apr  4: 3.26666 = 49/15    Apr 16: 3.26666 = 49/15
Mar 24: 1.42500 = 57/40    Apr  5: 3.38333 = 203/60   Apr 17: 3.38333 = 203/60
Mar 25: 1.93333 = 29/15    Apr  6: 3.32500 = 133/40   Apr 18: 3.46315 = 329/95
Mar 26: 2.33333 = 7/3      Apr  7: 3.32500 = 133/40   Apr 19: 3.86666 = 58/15
Mar 27: 2.90000 = 29/10    Apr  8: 3.38333 = 203/60   Apr 20: 3.32500 = 133/40
Mar 28: 3.26666 = 49/15    Apr  9: 3.26666 = 49/15    Apr 21: 2.85000 = 57/20
Mar 29: 3.38333 = 203/60   Apr 10: 3.38333 = 203/60   Apr 22: 2.41666 = 29/12
Mar 30: 3.32500 = 133/40   Apr 11: 3.26666 = 49/15    Apr 23: 1.86666 = 28/15
Mar 31: 3.32500 = 133/40   Apr 12: 3.38333 = 203/60   Apr 24: 1.45000 = 29/20
Apr  1: 3.38333 = 203/60   Apr 13: 3.32500 = 133/40   Apr 25: 0.73684 = 14/19
Apr  2: 3.26666 = 49/15    Apr 14: 3.32500 = 133/40

   Run with a single command-line argument of a year,  you get the dates of
Easter for a 75-year time span.  For example,  'easter 2010' gives

2010 Apr  4    2025 Apr 20    2040 Apr  1    2055 Apr 18    2070 Mar 30
2011 Apr 24    2026 Apr  5    2041 Apr 21    2056 Apr  2    2071 Apr 19
2012 Apr  8    2027 Mar 28    2042 Apr  6    2057 Apr 22    2072 Apr 10
2013 Mar 31    2028 Apr 16    2043 Mar 29    2058 Apr 14    2073 Mar 26
2014 Apr 20    2029 Apr  1    2044 Apr 17    2059 Mar 30    2074 Apr 15
2015 Apr  5    2030 Apr 21    2045 Apr  9    2060 Apr 18    2075 Apr  7
2016 Mar 27    2031 Apr 13    2046 Mar 25    2061 Apr 10    2076 Apr 19
2017 Apr 16    2032 Mar 28    2047 Apr 14    2062 Mar 26    2077 Apr 11
2018 Apr  1    2033 Apr 17    2048 Apr  5    2063 Apr 15    2078 Apr  3
2019 Apr 21    2034 Apr  9    2049 Apr 18    2064 Apr  6    2079 Apr 23
2020 Apr 12    2035 Mar 25    2050 Apr 10    2065 Mar 29    2080 Apr  7
2021 Apr  4    2036 Apr 13    2051 Apr  2    2066 Apr 11    2081 Mar 30
2022 Apr 17    2037 Apr  5    2052 Apr 21    2067 Apr  3    2082 Apr 19
2023 Apr  9    2038 Apr 25    2053 Apr  6    2068 Apr 22    2083 Apr  4
2024 Mar 31    2039 Apr 10    2054 Mar 29    2069 Apr 14    2084 Mar 26

   Run with two command line arguments of a month (3=March or 4=April) and
a day,  one gets a list of years from 1583 AD to 10000 AD when Easter would
fall on that day.  For example,  'easter 3 22' would give the following
list of years when Easter fell on 22 March,  the earliest day possible:

 1598 1693 1761 1818 2285 2353 2437 2505 2972 3029 3401 3496 3564 3648 3716
 4308 5299 5671 6043 6195 6263 6415 6635 6703 6798 6882 6950 7322 7474 7542
 7637 7789 7914 8161 8533 8685 8753 8848 8905 9125 9220 9372 9440 9812 9964

45 found over 8417 years
*/

void easter_date( const long year, int *month, int *day)
{
   const long a = year % 19L, b = year / 100L, c = year % 100L;
   const long d = b / 4L, e = b % 4L, f = (b + 8L) / 25L;
   const long g = (b - f + 1L) / 3L, h = (19L * a + b - d - g + 15L) % 30L;
   const long i = c / 4L, k = c % 4L, l = (32L + e + e + i + i - h - k) % 7L;
   const long m = (a + 11L * h + 22L * l) / 451L, tval = h + l - 7L * m + 114L;

   *month = (int)( tval / 31L);
   *day = (int)( tval % 31L) + 1;
}

/* From Meeus' _Astronomical Algorithms_.  This older method has
the 19-year Metonic cycle of lunar phases combined with the 28-year
cycle of the Julian calendar,  for a periodicity of 532 years.
I think it initially got use after the Council of Nicea in 325 AD
(there seems to be some confusion about this),  was then used by the
Catholic Church up to the Gregorian switchover in 1582,  and it
continues to be in use by the Orthodox Churches (so you can use
it to compute the date of Orthodox Easter). */

void easter_date_julian( const long year, int *month, int *day)
{
   const long a = year % 4L, b = year % 7L, c = year % 19L;
   const long d = (19L * c + 15L) % 30L, e = (2L * a + 4L * b + 34L - d) % 7L;
   const long tval = d + e + 114L;

   *month = (int)( tval / 31L);
   *day = (int)( tval % 31L) + 1;
}

#ifdef TEST_CODE

static unsigned gcd( unsigned a, unsigned b)
{
   while( a && b)
      {
      if( a > b)
         a %= b;
      else
         b %= a;
      }
   return( a | b);
}

#include <stdio.h>
#include <stdlib.h>

#if defined(_MSC_VER) && _MSC_VER < 1900
                      /* For older MSVCs,  we have to supply our own  */
                      /* snprintf().  See snprintf.cpp for details.  */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

int main( int argc, char **argv)
{
   int month, day;
   long year;

   if( argc == 2)
      {
      int i;
      const int n_across = 5, n_down = 15;

      for( i = 0; i < n_across * n_down; i++)
         {
         year = atol( argv[1]) + (i % n_across) * n_down + i / n_across;
         if( year >= 1583L)
            easter_date( year, &month, &day);
         else
            easter_date_julian( year, &month, &day);
         printf( "%ld %s %2d%s", year, (month == 3 ? "Mar" : "Apr"), day,
                  ((i + 1) % n_across) ? "    " : "\n");
         }
      }
   else if( argc == 3)
      {
      int n_found = 0;

      for( year = 1583; year < 10000; year++)
         {
         easter_date( year, &month, &day);
         if( month == atoi( argv[1]) && day == atoi( argv[2]))
            {
            printf( "%5ld", year);
            if( n_found++ % 15 == 14)
               printf( "\n");
            }
         }
      printf( "\n%d found over 8417 years\n", n_found);
      }
   else
      {
      unsigned march[32], april[32];
      unsigned i;

      printf( "Percentage frequency of Easter in the Gregorian calendar:\n");
      for( day = 0; day < 32; day++)
         april[day] = march[day] = 0;
      for( year = 0; year < 5700000; year++)
         {
         easter_date( year, &month, &day);
         if( month == 3)
            march[day]++;
         else
            april[day]++;
         }
      for( i = 0; i < 35; i++)
         {
         const int n_cols = 3;
         const int n_rows = (35 + n_cols - 1) / n_cols;
         unsigned remap = (i % n_cols) * n_rows + i / n_cols;
         unsigned tval, denom = 57000, g;
         const char *month_name;
         char buff[80];

         if( remap < 10)
            {
            remap += 22;
            month_name = "Mar";
            tval = march[remap];
            }
         else
            {
            remap -= 9;
            month_name = "Apr";
            tval = april[remap];
            }
         g = gcd( tval, denom);
         snprintf( buff, sizeof( buff), "%s %2d: %7.5f = %u/%u", month_name, remap,
                     (double)tval / (double)denom,
                     tval / g, denom / g);
         if( i % n_cols == n_cols - 1)
            printf( "%s\n", buff);
         else
            printf( "%-26s ", buff);
         }

      printf( "\nRun 'easter' with a year on the command line to get the date\n");
      printf( "of Easter for that year.  For example,  'easter 2008' will\n");
      printf( "get the output 'March 23'.  Alternatively,  give a month and\n");
      printf( "day on the command line to get the years between 1583 and 10000\n");
      printf( "when Easter will fall on that day.  For example,  'easter 3 23'\n");
      printf( "will produce a list of all years when Easter is on 23 March.\n");
      printf( "\nNote that Easter cannot occur before 22 March or after 25 April.\n");
      }
   return( 0);
}
#endif
