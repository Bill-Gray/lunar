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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "watdefs.h"
#include "date.h"

/* Code to test a (relatively) simple,  quick way to determine in which half
year a given MJD falls.  This is of use in 'delta_t.cpp' in figuring out
guesstimates of when future leap seconds should be inserted.  See the
'td_minus_utc()' function.

   Having been written and used to test the handling of future leap seconds,
this code is really of only historical use.  I'm leaving it around in order
to document how this rather strange algorithm works.

   The basic idea is this.  400 Gregorian years will contain 400*365 "regular"
days,  plus 397 leap days,  for a total of 146097 days.  Therefore,  given a
number of days past "year zero",  the year will be approximately

year = day * 400 / 146097

   Note that,  in what follows,  the above math is promoted to 64-bit ints
to avoid overflow.  1 Jan 0 (first day of the year zero) is JD 1721057.5.
Round this up,  and one finds that day = mjd + 2400000 - 1721058.  Using
this,  the above formula can sometimes yield a year that is one too high
for dates in late December;  hence the computation of 'low' as the MJD of
the corresponding year and the test to see "if( mjd < low)".

   So now we know what year corresponds to the given MJD.  Next problem is,
are we before or after 1 July of that year?  To find this out,  we compute
the MJD for 1 January of the following year,  then compute the MJD of the
preceding 1 July.  There can't be leap days over that span,  so all we
need do is to subtract 31 + 31 + 30 + 31 + 30 + 31 = 183 days to back up
from 1 Jan to 1 July.

   We then know in which half-year we are.  In the actual Delta-T code,
we can then average the MJDs of the ends of that time span to get the MJD
for the _middle_ of the time span,  compute Delta-T at that point,  and come
up with a UTC offset that'll keep |UTC-UT| < .5 second;  i.e.,  keep UTC
and UT as close together as possible.

   This is admittedly something of a kludge.  But anything done to estimate
future leap seconds has to be a kludge;  there's no "good" way to do it.

   See 'delta_t.cpp' for an explanation of the following macro. */

#define BASE_YEAR 19999999999L
#define JAN_1( YEAR) (((YEAR) * 365L + ((YEAR) + BASE_YEAR) / 4L - ((YEAR) + BASE_YEAR) / 100L \
                         + ((YEAR) + BASE_YEAR) / 400L) - 678940L \
                         - ((BASE_YEAR + 1L) / 400L) * 97L)

// unsigned day = mjd + 1931367u + 2400000u;
// int year = (int)( day * 400u / 146097u) - 10000;

int main( const int argc, const char **argv)
{
   int mjd = atoi( argv[1]);

   if( argc == 2)
      {
      int day = mjd + 2400000 - 1721058;
      int year = (int)( (int64_t)day * (int64_t)400 / (int64_t)146097);
      int low, high, july_1;

      printf( "Year = %d\n", year);
      low = JAN_1( year);
      printf( "MJD 1 Jan %d = %ld\n", year, JAN_1( year));
      if( mjd < low)
         {
         year--;
         low = JAN_1( year);
         printf( "MJD 1 Jan %d = %ld\n", year, JAN_1( year));
         }
      high = JAN_1( year + 1);
      printf( "MJD 1 Jan %d = %d\n", year + 1, high);
                   /*  jul  aug  sep  oct  nov  dec */
      july_1 = high - (31 + 31 + 30 + 31 + 30 + 31);
      printf( "MJD 1 Jul %d = %d\n", year, july_1);
      if( mjd < july_1)
         printf( "In first half of %d\n", year);
      else
         printf( "In second half of %d\n", year);
      }
   else                             /* Test to see:  over a given range of  */
      {                             /* MJDs,  how often does the above year */
      int n_less = 0, n_more = 0;   /* determination land "spot on"?  How   */
                                    /* often must we go to the next year?   */
      while( mjd < atoi( argv[2]))
         {
         int day = mjd + 2400000 - 1721058;
         int year = (int)( (int64_t)day * (int64_t)400 / (int64_t)146097);

         if( JAN_1( year) > mjd)
            n_more++;
         if( JAN_1( year + 1) <= mjd)
            n_less++;
         mjd++;
         }
      printf( "%d cases went over;  %d went under\n", n_more, n_less);
      }
   return( 0);
}
