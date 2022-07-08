#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "watdefs.h"
#include "mpc_func.h"
#include "date.h"
#include "afuncs.h"

/*
   The obsID field in ADES starts with eight mutant-hex (base 62)
digits that give the time the observation was received at MPC.
Mike Rudenko pointed me to

https://minorplanetcenter.net/decode_obsid.py

   from which I learned that :

   The first three digits,  when decoded,  give the year,  month,
and day of submission,  packed as

value(0,1,2) = (year - 1800) * 12 * 31 + (month - 1) * 31 + (day-1)

   The next three digits,  when decoded,  give

value(3,5) = seconds since 00:00 UT

   and therefore ranges from 0 to 86399 (MTX in base 62).  And
finally,  the last two digits give

value(6,7) = milliseconds

   and is therefore from 0 to 999 (G7 in base 62).  Valid times can
run from

00000000 = 1800-Jan-01 00:00:00.000
zzzMTXG7 = 2440-Aug-31 23:59:59.999

   unless that day has a leap second (it _is_ at the end of a
month,  and we'll be having leap seconds monthly by then,  so
that will probably happen.)   In that case,  the scheme can be
extended by one second :

zzzMTYG7 = 2440-Aug-31 24:00:00.999

   The need to handle leap second madness may be why MPC went with
this scheme,  instead of the straightforward 'milliseconds since
an epoch' one I'd expected.

   Some test cases :

< Date ><SubNum>sb<ObsNu>
K6yCeR0100005WG401000000A        2007-10-30T13:30:35.001_00005WG4_01
LP6CWP340000E5ZY010000QZe        2021-03-07T13:22:17.190_0000E5ZY_01
LQBB0Q8k0000EAnd0100006Ne        2021-05-12T11:45:10.542_0000EAnd_01

   David Bell provided some additional info :

1-8 - ISO date to millisecond in base-62

9-16 - base62 counter that increases for each submission

17-18 - base62 sub-batch number within the main batch, usually “01”
unless a single submission included multiple header blocks (different
telescopes, different observers for each subbatch, etc.)

19-25 - base62 counter for the observation number within the sub-batch */

static int show_usage( void)
{
   printf( "usage : ./mpc_time (obsID)\n"
           "Can be full obsID or,  if you just want the submission time,\n"
           "supply the first eight characters of the obsID.\n");
   return( -1);
}

int main( const int argc, const char **argv)
{
   const char *id = (argc == 2 ? argv[1] : "            ");
   const int ymd = get_mutant_hex_value( id, 3);
   const int hms = get_mutant_hex_value( id + 3, 3);
   const int millisec = get_mutant_hex_value( id + 6, 2);
   const int year = ymd / (31 * 12) + 1800;
   const int month = (ymd / 31) % 12 + 1;
   const int day = ymd % 31 + 1;
   const int n_days = days_in_month( month, year, CALENDAR_GREGORIAN);
   char buff[80];

   if( argc != 2)
      return( show_usage( ));
   snprintf( buff, sizeof( buff), "%04d-%02d-%02d %02d:%02d:%02d.%03d",
            year, month, day,
            hms / 3600, (hms / 60) % 60, hms % 60,
            millisec);
   printf( "%s = JD %f\n", buff,
               get_time_from_string( 0., buff, FULL_CTIME_YMD, NULL));
   if( strlen( id) >= 16)
      printf( "Submission counter: %.8s (mutant hex)\n", id + 8);
   if( strlen( id) >= 18)
      printf( "Sub-batch number within main batch counter: %.2s (mutant hex)\n", id + 16);
   if( strlen( id) >= 25)
      printf( "Observation number within the sub-batch: %.7s (mutant hex)\n", id + 18);
   if( day > n_days)
      printf( "WARNING: There are only %d days in that month\n", n_days);
   if( millisec >= 1000)
      printf( "WARNING: milliseconds should be 999 or less\n");
   if( hms > 86400)        /* allow for one possible leap second */
      printf( "WARNING: time of day is out of range\n");
   if( hms == 86400)
      {
      if( (month % 6) || day != n_days)
         printf( "WARNING: This is an invalid leap second.  It's not on the\n"
                 "last day of June or December.\n");
      else
         printf( "WARNING: This is a leap second.  Theoretically,  it's allowed,\n"
                 "but my code and MPC's will probably break if this submission\n"
                 "time is used.\n");
      }
         /* Further checking could be done.  Times in the future or
         before the time of the actual observation should be flagged. */
   return( 0);
}
