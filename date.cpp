/* date.cpp: date/time/calendar conversions

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

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include "watdefs.h"
#include "get_bin.h"
#include "date.h"

/* General calendrical comments:

This code supports conversions between JD and eight calendrical systems:
Julian,  Gregorian,  Hebrew,  Islamic,  Jalaali (Persian),  Chinese,
(French) revolutionary,  and Modern Persian.  Comments pertaining to
specific calendars are found near the code for those calendars.

    For each calendar,  there is a "get_(calendar_name)_year_data( )"
function, used only within this source code. This function takes a particular
year number,  and computes the JD corresponding to "new years day" (first day
of the first month) in that calendar in that year.  It also figures out the
number of days in each month of that year,  returned in the array
month_data[].  There can be up to 13 months,  because the Hebrew and Chinese
calendars can include an "intercalary month";  also,  the five or six extra
days at the end of the French Revolutionary calendar are considered as "the
thirteenth month" in this code.

   If a month doesn't exist,  then the month_data[] entry for it will be zero.
 Thus,  in the Gregorian and Julian and Islamic calendars, month_data[12] is
always zero,  since these calendars have only 12 months. The same will happen
in Chinese and Hebrew years that lack intercalary months,  meaning about 63%
of the time.  So far,  all calendars involve twelve or thirteen months.
(Some involve more;  the Baha'i calendar,  for example,  involves 19,  and
the Celtic "year-and-a-day" calendar involves 13 months of 28 days plus one
or two leftover days that would be treated as a fourteenth month.  If the
code is ever revised to handle this,  N_MONTHS will have to be re-defined,
and other changes probably made.)

The next level up is the get_calendar_data( ) function,  which (through
the wonders of a switch statement) can get the JD of New Years Day and
the array of months for any given year for any calendar.  Above this
point,  all calendars can be treated in a common way;  one is shielded
from the oddities of individual calendrical systems.  The
get_calendar_data( ) function is still local to this file.

Finally,  at the top level,  we reach the only two functions that are
exported for the rest of the world to use:  dmy_to_day( ) and day_to_dmy( ).
The first takes a day,  month,  year,  and calendar system.  It calls
get_calendar_data( ) for the given year,  adds in the days in the months
intervening New Years Day and the desired month,  and adds in the day
of the month,  and returns the resulting Julian Day.

day_to_dmy( ) reverses this process.  It finds an "approximate" year
corresponding to an input JD,  and calls get_calendar_data( ) for
that year.  By adding all the month_data[] values for that year,  it
can also find the JD for the _end_ of that year;  if the input JD is
outside that range,  it may have to back up a year or add in a year.
Once it finds "JD of New Years Day < JD < JD of New Years Eve",  it's
a simple matter to track down which month and day of the month corresponds
to the input JD.
*/

/* The following mod( ) function returns the _positive_ remainder after */
/* a division.  Annoyingly,  if x < 0,  then x % y <= 0;  thus,  this   */
/* function is needed for things such as determining a day of the week. */

static long mod( const long x, const long y)
{
   long rval = x % y;

   if( rval < 0L)
      rval += y;
   return( rval);
}

/* Begin:  Gregorian and Julian calendars (combined for simplicity) */

/* It's common to implement Gregorian/Julian calendar code with the     */
/* aid of cryptic formulae,  rather than through simple lookup tables.  */
/* For example,  consider this formula from Fliegel and Van Flandern,   */
/* to convert Gregorian (D)ay, (M)onth, (Y)ear to JD:                   */

/* JD = (1461*(Y+4800+(M-14)/12))/4+(367*(M-2-12*((M-14)/12)))/12       */
/*       -(3*((Y+4900+(M-14)/12)/100))/4+D-32075                        */

/* The only way to verify that they work is to feed through all possible    */
/* cases.  Personally,  I like to be able to look at a chunk of code and    */
/* see what it means.  It should resemble the Reformation view of the       */
/* Bible:  anyone can read it and witness the truth thereof.                */

/* Several of these calendars have intercalary months,  so we gotta allow */
/* for up to thirteen months.                                             */

/* 3 Jan 2004:  Modified 'set_month_name( )' so that it resets a pointer to  */
/* an externally-provided month name,  rather than providing the storage for */
/* the names in locally static space.                                        */

/* 20 Aug 2008:  made 'month_names' and 'day_of_week_names' static const.   */

#define N_MONTHS 13

static const char *month_names[N_MONTHS] = { "Jan", "Feb", "Mar",
                          "Apr", "May", "Jun",
                          "Jul", "Aug", "Sep",
                          "Oct", "Nov", "Dec", NULL };

static const char *day_of_week_names[7] = { "Sun", "Mon", "Tue", "Wed",
                               "Thu", "Fri", "Sat"};

const char * DLL_FUNC set_month_name( const int month, const char *new_name)
{
   assert( month >= 1 && month <= N_MONTHS);
   if( new_name)
      month_names[month - 1] = new_name;
   return( month_names[month - 1]);
}

const char * DLL_FUNC set_day_of_week_name( const int day_of_week,
                                            const char *new_name)
{
   assert( day_of_week >= 0 && day_of_week < 7);
   if( new_name)
      day_of_week_names[day_of_week] = new_name;
   return( day_of_week_names[day_of_week]);
}

#define JUL_GREG_CALENDAR_EPOCH 1721060L

static void get_jul_greg_year_data( const long year, long *days,
                                       char *month_data, const int julian)
{
   static const char months[13] =
                   { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 0 };
                 /*  Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec */

   if( year >= 0L)
      {
      *days = year * 365L + year / 4L;
      if( !julian)
         *days += -year / 100L + year / 400L;
      }
   else
      {
      *days = year*365L + (year-3L) / 4L;
      if( !julian)
         *days += - (year-99L) / 100L + (year-399L) / 400L;
      }

   if( julian)
      *days -= 2L;
   memcpy( month_data, months, 13);
   if( !(year % 4))
      if( (year % 100L) || !(year % 400L) || julian)
         {
         month_data[1] = 29;
         (*days)--;
         }
   *days += JUL_GREG_CALENDAR_EPOCH + 1;
}

/* End:  Gregorian and Julian calendars */

/* Begin:  Islamic calendar */

/* 23 Feb 2004:  Revised both Islamic and Hebrew calendars to remove the */
/* tables of leap years,  relying instead on algorithmic versions.  That */
/* let me remove a few dozen lines of code.                              */

#define ISLAMIC_CALENDAR_EPOCH 1948086L

static void get_islamic_year_data( const long year, long *days,
                                                      char *month_data)
{
   const long thirty_islamic_years = 10631L;
   const long year_within_cycle = mod( year, 30L);
   const long thirty_year_cycles = (year - year_within_cycle) / 30L;
   long rval;
   const long tval = year_within_cycle * 11 + 3;

   rval = ISLAMIC_CALENDAR_EPOCH +
            thirty_year_cycles * thirty_islamic_years +
                               year_within_cycle * 354L;
   month_data[12] = 0;
   month_data[11] = (char)( 29 + ((tval % 30) > 18));
   rval += tval / 30;
   *days = rval;
            /* The Islamic calendar alternates between 30-day and 29-day */
            /* months for the first eleven months;  the twelfth is 30    */
            /* days in a leap year,  29 otherwise (see above).           */
   for( unsigned i = 0; i < 11; i++)
      month_data[i] = (char)( 30 - (i % 2));
}

/* End:  Islamic calendar */

/* Begin:  Hebrew calendar */

/* See p 586,  _Explanatory Supplement_,  for explanation. */
/* There are 1080 Halakim,  or 'parts',  each 3.33 seconds long,  in  */
/* an hour.  So: */

#define HALAKIM_IN_DAY (24L * 1080L)
#define HEBREW_CALENDAR_EPOCH 347996L

static long lunations_to_tishri_1( const long year)
{
   const long year_within_cycle  = mod( year - 1, 19L);
   const long full_nineteen_year_cycles = (year - 1 - year_within_cycle) / 19L;
   long rval;

   rval = full_nineteen_year_cycles * 235L + year_within_cycle * 12L;
   rval += (year_within_cycle * 7 + 1) / 19;
   return( rval);
}

   /* One lunation is 29 days, 13753 halakim long,  which equals 765433
halakim. In theory,  one could just write

total_halakim = 765433 * lunations;
*days    += total_halakim / HALAKIM_IN_DAY;
*halakim += total_halakim % HALAKIM_IN_DAY;

   The only problem with this is that (assuming 32-bit integers)
'total_halakim' overflows after 2805 lunations,  or about 220 years.
Ideally,  I'd have switched to 64-bit integers,  but (a) such are not
always available and (b) it's easy enough to get around this problem.
The trick is to recognize that 25920 lunations is exactly equal to
765433 days.  This cycle has no name that I know of,  and no real
significance outside of this particular function... but it _does_ let us
write simpler code that won't get wrong answers for large or negative
numbers of lunations.  Let's call 25920 lunations a "glumph."  We figure
out how many glumphs have passed and our location within that glumph,
and the rest is easy.

   As a side effect of all this : since there are 235 lunations
every 19 years,  and 765433 halakim per lunation,  there are exactly
179876755 halakim every 19 years,  or 179876755/(1080*24) days =
35975351/(216*24) days.  Thus,  35975351 days are exactly 19*216*24
= 98496 years;  that is to say,  the Hebrew calendar repeats itself
exactly every 98496 years.  Put another way,  the average length of
a year is exactly 35975351/98496 = 365.24682220598... days.  */

static void lunations_to_days_and_halakim( const long lunations, long *days,
               long *halakim)
{
   const long lunation_within_glumph = mod( lunations, 25920L);
   const long curr_glumph = (lunations - lunation_within_glumph) / 25920L;

   *days += curr_glumph * 765433L + lunation_within_glumph * 29L;
   *halakim += lunation_within_glumph * 13753L;

               /* Now make sure excess halakim carry over correctly: */
   *days += *halakim / HALAKIM_IN_DAY;
   *halakim %= HALAKIM_IN_DAY;
}

static void find_tishri_1( const long year, long *days, long *halakim)
{
            /* Set days and halakim to the 'epoch':  1 Tishri 1 = 2 5604 */
   *days = 2L;
   *halakim = 5604L;
   lunations_to_days_and_halakim( lunations_to_tishri_1( year), days, halakim);
}

static int is_hebrew_leap_year( const long year)
{
   return( mod( year * 7 - 6, 19) >= 12);
}

/* Certain aspects of get_hebrew_year_data( ) will definitely fail for */
/* years before zero... something will have to be done about that.     */

static void get_hebrew_year_data( const long year, long *days, char *month_data)
{
   for( int i = 0; i < 2; i++)
      {
      long day, halakim;

      find_tishri_1( year + i, &day, &halakim);
               /* Check dehiyyah (c): */
      if( mod( day, 7L) == 3 && halakim >= 9L * 1080L + 204L &&
                    !is_hebrew_leap_year( year + i))
         day += 2;
      else               /* Check dehiyyah (d): */
         if( mod( day, 7L) == 2 && halakim >= 15L * 1080L + 589L &&
                                is_hebrew_leap_year( year - 1 + i))
            day++;
      else
         {
         if( halakim > 18L * 1080L)
            day++;
         if( mod( day, 7L) == 1 || mod( day, 7L) == 4 || mod( day, 7L) == 6L)
            day++;
         }
      days[i] = day + HEBREW_CALENDAR_EPOCH;
      }
   int year_length = (int)( days[1] - days[0]);
   if( month_data)
      {
      for( int i = 0; i < 6; i++)                 /* "normal" lengths */
         month_data[i] = month_data[i + 7] = (char)( 30 - (i & 1));
      if( is_hebrew_leap_year( year))
         {
         month_data[5] = 30;     /* Adar I is bumped up a day in leap years */
         month_data[6] = 29;
         }
      else                       /* In non-leap years,  Adar II doesn't    */
         month_data[6] = 0;      /* exist at all;  set it to zero days     */
      if( year_length == 353 || year_length == 383)      /* deficient year */
         month_data[2] = 29;
      if( year_length == 355 || year_length == 385)      /* complete year  */
         month_data[1] = 30;
      }
}

/*  Some test cases:  16 Av 5748 AM (16 12 5748) = 30 Jul 1988 Gregorian */
/*                 14 Nisan 5730 AM (14 8 5730) = 20 Apr 1970 Gregorian */
/*                  1 Tishri 5750 AM (1 1 5750) = 30 Sep 1989 Gregorian */

/* End:  Hebrew calendar */

/* Begin:  (French) Revolutionary calendar */

/*
The French Revolutionary calendar is simplest,  in some respects;
you just have twelve months,  each of 30 days,  with five or six
"unattached" days at the end of the year.  The only real problem
is in handling leap years.  There are no fewer than four possible
ways of doing this.

   The calendar was originally defined to have 1 Vendemiare,  "New
Years Day",  to be the date on which the autumnal equinox occurred
as seen from Paris.  It was also defined to occur once every four
years.  You can't have both,  and the contradiction was quickly
noticed and resulted in the four schemes:

   (1) Keep it matched to the autumnal equinox. In this case,  leap
years are usually four years apart but are sometimes five (similar
to the Jalaali calendar).  It also becomes an "observational" calendar;
extension into the distant future (and past) becomes difficult,
because the uncertainty in the earth's rotation means we can't be
totally sure as to the date on which the equinox will fall.  The pattern
of leap years is far from intuitively obvious (you can't just look at
a year and see if it's divisible by four).  And at the time,  computing
the date of the autumnal solstice wasn't a trivial exercise. However,
this _was_ the scheme in actual use during the short active lifetime of
this calendar, resulting in Years 3, 7, 11,  and 15 AR being leap years.

   (2) The Gregorian scheme:  leap years are those divisible by four,
except for those divisible by 100,  except for those divisible by
400.  Alternatively...

   (3) ...the same thing,  with the addition of "except for those
divisible by 4000",  something occasionally suggested for the Gregorian
calendar as well.

   (4) A scheme in which leap years are those divisible by four,
except for those divisible by 128.  This slight deviation from the
Gregorian scheme,  of "divisible by four,  unless divisible by
100,  unless divisible by 400",  is slightly simpler and gives a
calendar that is _much_ closer to the true tropical year.

   In practice,  the only real certainty appears to be that years
3, 7, 11,  and 15 were leap years.  After that,  one of the above
schemes was to be implemented,  but the calendar was abolished
before that happened.

   I really doubt the likelihood of schemes (2) and (3).  A revolution
so devoted to revising every aspect of human existence that it
changed names of all months,  "regularized" each to be 30 days,  and
made a week ten days long,  probably went out of its way not to
resemble earlier calendars proposed by a Pope.  The defects of (1)
being somewhat apparent,  my tendency is to go with the 4/128 rule.
I've found code and authors supporting all four schemes.  By default,
the following code uses the 4/128 rule,  but as you'll see,  all
four schemes can be implemented by #ifdef'fing chunks of code.

   A 'proleptic' calendar wasn't defined,  to my knowledge...
however,  "BR" (Before the Revolution) years are handled logically,
for all four schemes,  in this code.
*/

#ifdef REVOLUTIONARY_128_RULE

#define REVOLUTIONARY_CALENDAR_EPOCH 2375475L

static long jd_of_french_rev_year( long year)
{
   long rval = REVOLUTIONARY_CALENDAR_EPOCH + year * 365L;

   if( year >= 20)
      year--;
#ifdef GREGORIAN_REVOLUTIONARY
   rval += (long)(year / 4 - year / 100 + year / 400);
#else
   rval += (long)(year / 4 - year / 128);
#endif

   if( year <= 0L)
      rval--;
   return( rval);
}

#else    /* REVOLUTIONARY_128_RULE */

#define REVOLUTIONARY_CALENDAR_EPOCH 2375475L
/* Following macros are still valid,  but not used;
#define LOWER_REVOLUTIONARY_YEAR -1007
#define UPPER_REVOLUTIONARY_YEAR 611L
it's included for reference,  and commented out */

static long jd_of_french_rev_year( const long revolutionary_year)
{
   static const short breaks[9] = { -814, -492, -331, -108,   0, 144,
                                      301, 487, 611 };
   static const short deltas[9] = {  405,  439,  498,  469, 419, 393,
                                      322, 184,  92 };
   long rval = 0;

   for( int i = 0; !rval; i++)
      if( revolutionary_year < (long)breaks[i] || i == 8)
         {
         rval = REVOLUTIONARY_CALENDAR_EPOCH + revolutionary_year * 365L +
                     ((long)deltas[i] + revolutionary_year * 683L) / 2820L;
         if( i < 5)  /* zero point drops one day in first five blocks */
            rval--;
         }
   return( rval);
}
#endif    /* REVOLUTIONARY_128_RULE */

static void get_revolutionary_year_data( const long year, long *days,
                                   char *month_data)
{
   days[0] = jd_of_french_rev_year( year);
   days[1] = jd_of_french_rev_year( year + 1);
   memset( month_data, 30, 12);
            /* There are twelve months of 30 days each,  followed by  */
            /* five (leap years,  six) days;  call 'em an extra       */
            /* thirteenth "month",  containing all remaining days:    */
   month_data[12] = (char)( days[1] - days[0] - 360L);
}

/* End:  (French) Revolutionary calendar */

/* Begin:  Persian (Jalaali) calendar */

#define JALALI_ZERO 1947954L
#define LOWER_PERSIAN_YEAR -1096
#define UPPER_PERSIAN_YEAR 2327

static long jalali_jd0( const int jalali_year)
{
   static const short breaks[12] = { -708, -221,   -3,    6,  394,  720,
                                      786, 1145, 1635, 1701, 1866, 2328 };
   static const short deltas[12] = { 1108, 1047,  984, 1249,  952,  891,
                                      930,  866,  869,  844,  848,  852 };
   long rval;

   if( jalali_year < LOWER_PERSIAN_YEAR)
      return( -1L);           /* out of valid range */
   for( int i = 0; i < 12; i++)
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

/* 7 April 2004:  added the 'modern Persian calendar',  which follows the  */
/* pattern of the astronomically-based Persian (Jalaali) calendar closely, */
/* but not exactly.  The 'modern' (algorithmic) flavor has a pattern of    */
/* 683 leap years over a 2820-year cycle,  in a manner that collapses      */
/* nicely into a few lines of code.                                        */

static long persian_modern_jd0( const long year)
{
   const long persian_epoch = 1948320L;
   const long epbase = year - 474L;
   const long epyear = 474L + mod( epbase, 2820L);

   return( (epyear * 31 - 5) / 128 + (epyear - 1) * 365
               + ((year - epyear) / 2820) * 1029983 + persian_epoch);
}

/* 7 April 2004:  modified 'get_jalali_year_data()' so that it can use   */
/* either the traditional astronomical Jalaali-based start of the year,  */
/* or the 'modern' (algorithmic) Persian calendar.  If we're outside the */
/* range covered by the Jalaali algorithm in this code,  we fall back on */
/* use of the algorithmic rule.                                          */

static void get_jalali_year_data( const long year, long *days,
                                   char *month_data, int is_modern)
{
   if( year < LOWER_PERSIAN_YEAR || year > UPPER_PERSIAN_YEAR)
      is_modern = 1;
   if( is_modern)
      {
      days[0] = persian_modern_jd0( year) + 1L;
      days[1] = persian_modern_jd0( year + 1L) + 1L;
      }
   else
      {
      days[0] = jalali_jd0( year) + 1L;
      days[1] = jalali_jd0( year + 1L) + 1L;
      }
            /* The first six months have 31 days.  The next five have 30  */
            /* days.  The last month has 29 days in ordinary years,  30   */
            /* in leap years.                                             */
   memset( month_data, 31, 6);
   memset( month_data + 6, 30, 5);
   month_data[11] = (char)( days[1] - days[0] - 336L);
   month_data[12] = 0;      /* always a twelve-month/year calendar */
}

/* End:  Persian (Jalali) calendar */

/* Start:  Chinese calendar */

/* The Chinese calendar poses some particularly sticky problems,  and is
a real mess to compute.  Therefore,  I wrote a separate program to
"pre-compile" a calendar and store it in a file,  CHINESE.DAT.  To handle
the Chinese calendar,  one must load that file into a buffer,
chinese_calendar_data.  The following code then simply "un-crunches" the
data.

   NOTE:  The odd arrangement of intercalary months in the Chinese calendar
means that one must be careful about month numbering.  The solution I've
used in this code is to return a twelve or thirteen-month calendar,  much
as with the Chinese calendar;  the external int chinese_intercalary_month
is used to indicate which of the 13 is intercalary (it's left at zero if
there is no intercalary month.)  See JD.CPP for an example of its use.

   NOTE:  There is no 'official' numbering scheme for the Chinese calendar.
Traditionally,  it just ran on a sixty-year cycle.  There are three
different schemes in use.  I've used one in which the Gregorian year 2000
corresponds to a Chinese year 4637.  Some people think the Chinese
calendar was created sixty years earlier,  and that therefore 2000 (Greg.)
= 4397 (Chinese).  Just to make things still more confusing,  some people
add another year to this: 2000 (Greg.) = 4398 (Chinese).  So you may see
Chinese "year" numbers that are 60 or 61 greater than those used in this
code.

   NOTE:  Some parts of the following code extract integer data without
proper regard to byte-order.  If you compile this on "wrong-endian"
machines,  this will have to be fixed!  (Not a big deal;  three lines,
marked below,  need changes.)

   The "packing" works as follows.  The data for each year is stored in
24 bits.  The 13 least significant bits are set (or unset) to indicate
30-day (or 29-day) months.  Shifting the "packed value" down by 13 thus
leaves an 11-bit quantity.

   That value,  modulo 14,  tells you which month in that year is the
intercalary one.  (If it's zero,  there is no intercalary month,  and
this is a twelve-month year.)  Dividing the value by 14 gives you an
"offset",  and the JD of the Chinese New Year is then given as

jd = 365 * year + year / 4 + CHINESE_CALENDAR_EPOCH + offset

Thus,  the offset can be a value from 0 to (2^11 / 14) = 186. */

#ifdef LOAD_CHINESE_CALENDAR_DATA_FROM_FILE
static const unsigned char *chinese_calendar_data = NULL;

void DLL_FUNC set_chinese_calendar_data( const void *cdata)
{
   chinese_calendar_data = (const unsigned char *)cdata;
}
#else
   #include "chinese.h"
#endif    /* #ifdef LOAD_CHINESE_CALENDAR_DATA_FROM_FILE */

static int chinese_intercalary_month = 0;

#define CHINESE_CALENDAR_EPOCH 757862L

static int get_chinese_year_data( const long year, long *days,
                                   char *month_data)
{
   int32_t packed_val;
   char tbuff[4];

#ifdef LOAD_CHINESE_CALENDAR_DATA_FROM_FILE
   if( !chinese_calendar_data)
      return( -1);
#endif

   int index = (int)year - get16sbits( chinese_calendar_data + 2);
   const int n_years = get16sbits( chinese_calendar_data);

   if( index < 0 || index >= n_years)
      return( -2);
   memcpy( tbuff, chinese_calendar_data + 4 + 3 * index, 3);
   tbuff[3] = 0;
   packed_val = get32sbits( tbuff);
   for( int i = 0; i < 13; i++)
      month_data[i] = (char)(((packed_val >> i) & 1L) ? 30 : 29);
   chinese_intercalary_month = (int)( (packed_val >> 13) % 14L);
   if( chinese_intercalary_month)
      chinese_intercalary_month++;
   else                       /* if there's no intercalary month... */
      month_data[12] = 0;     /* ...then the '13th month' has zero days */
   *days = year * 365L + year / 4L + CHINESE_CALENDAR_EPOCH +
                               (packed_val >> 13) / 14L;

   return( 0);
}

/* End:  Chinese calendar */

static int get_calendar_data( const long year, long *days, char *month_data,
               const int calendar)
{
   int rval = 0;

   memset( month_data, 0, N_MONTHS);
   switch( calendar)
      {
      case CALENDAR_GREGORIAN:
      case CALENDAR_JULIAN:
         get_jul_greg_year_data( year, days, month_data,
                            (calendar == CALENDAR_JULIAN));
         break;
      case CALENDAR_HEBREW:
         get_hebrew_year_data( year, days, month_data);
         break;
      case CALENDAR_ISLAMIC:
         get_islamic_year_data( year, days, month_data);
         break;
      case CALENDAR_REVOLUTIONARY:
         get_revolutionary_year_data( year, days, month_data);
         break;
      case CALENDAR_PERSIAN:
         get_jalali_year_data( year, days, month_data, 0);
         if( year < LOWER_PERSIAN_YEAR || year > UPPER_PERSIAN_YEAR)
            rval = -1;
         break;
      case CALENDAR_CHINESE:
         rval = get_chinese_year_data( year, days, month_data);
         break;
      case CALENDAR_MODERN_PERSIAN:
         get_jalali_year_data( year, days, month_data, 1);
         break;
      default:
         rval = -1;
         break;
      }              /* days[1] = JD of "New Years Eve" + 1;  that is,    */
   if( !rval)        /* New Years Day of the following year.  If you have */
      {              /* days[0] <= JD < days[1],  JD is in the current year. */
      days[1] = days[0];
      for( int i = 0; i < N_MONTHS; i++)
         days[1] += month_data[i];
      }
   return( rval);
}

int DLL_FUNC get_chinese_intercalary_month( void)
{
   return( chinese_intercalary_month);
}

#define RETURN_DAYS_IN_MONTH    -999

/* dmy_to_day( ) just gets calendar data for the current year,  including
the JD of New Years Day for that year.  After that,  all it has to do is
add up the days in intervening months,  plus the day of the month,  and
it's done:  */

long DLL_FUNC dmy_to_day( const int day, const int month, const long year,
                            const int calendar)
{
   char mdata[N_MONTHS];
   long jd;
   long year_ends[2];
   int calendar_to_use = calendar;
   int rval;

   if( calendar == CALENDAR_JULIAN_GREGORIAN)
      {
      if( year > 1582 || (year == 1582 &&
                  (month > 10 || (month == 10 && day > 5))))
         calendar_to_use = CALENDAR_GREGORIAN;
      else
         calendar_to_use = CALENDAR_JULIAN;
      }
   rval = get_calendar_data( year, year_ends, mdata, calendar_to_use);
   if( !rval)
      {
      if( day == RETURN_DAYS_IN_MONTH)
         return( mdata[month - 1]);
      jd = year_ends[0];
      for( int i = 0; i < month - 1; i++)
         jd += mdata[i];
      jd += (long)day - 1;
      }
   else
      jd = 0;
   return( jd);
}

int DLL_FUNC days_in_month( const int month, const long year,
                            const int calendar)
{
   return( dmy_to_day( RETURN_DAYS_IN_MONTH, month, year, calendar));
}

/* This usually gets you the correct year for a given JD,  but is
sometimes off by one around the New Year of the calendar in question.
Which is why the subsequent day_to_dmy( ) function sometimes finds it
has to move ahead or back up by one year.

   The French Revolutionary and both Persian calendars have (over the
long range) years of 365 + 683/2820 days;  i.e.,  2820 years have
2820 * 365 + 683 days.  The other calendars have similar long-range
exact recurrences in which n1 years will have n2 days.  The exact
recurrences for the Persian,  Jalali,  Hebrew,  and French calendars
would overflow 32-bit arithmetic,  so "almost" recurrences are used
if longs aren't 64 bits.  The error would become noticeable after
about a billion years,  but on 32-bit systems,  we're limited to
+/- 2^31 days = about 5.8 million years anyway.     */

/* #define LONGS_ARE_64_BITS */

static long approx_year( long jd, const int calendar)
{
   long year, n1 = 0, n2 = 0, calendar_epoch, day_in_cycle;

   switch( calendar)
      {
      case CALENDAR_GREGORIAN:
         calendar_epoch = JUL_GREG_CALENDAR_EPOCH;
         n1 = 400;               /* 400 Gregorian years contain 400 */
         n2 = 400 * 365 + 97;    /* 'normal' years plus 97 leap days */
         break;
      case CALENDAR_JULIAN:
         calendar_epoch = JUL_GREG_CALENDAR_EPOCH - 2;
         n1 = 4;                 /* The Julian calendar just repeats */
         n2 = 365 * 4 + 1;       /* every four years */
         break;
      case CALENDAR_HEBREW:
         calendar_epoch = HEBREW_CALENDAR_EPOCH - 235;
#ifdef LONGS_ARE_64_BITS
         n1 = 98496;          /* exact values which overflow on 32 bits */
         n2 = 35975351;       /*                                        */
#else
         n1 = 944;
         n2 = 344793;
#endif
         break;
      case CALENDAR_ISLAMIC:
         calendar_epoch = ISLAMIC_CALENDAR_EPOCH - 1;
         n1 = 30;             /* 30 Islamic years = 10631 days,  exactly */
         n2 = 10631;
         break;
      case CALENDAR_REVOLUTIONARY:
         calendar_epoch = REVOLUTIONARY_CALENDAR_EPOCH - 1;
#ifdef LONGS_ARE_64_BITS
         n1 = 2820;           /* exact values which overflow on 32 bits */
         n2 = 2820 * 365 + 683;
#else
         n1 = 2147;
         n2 = 784175;
#endif
         break;
      case CALENDAR_PERSIAN:
      case CALENDAR_MODERN_PERSIAN:
         calendar_epoch = JALALI_ZERO + 1;
#ifdef LONGS_ARE_64_BITS
         n1 = 2820;           /* exact values which overflow on 32 bits */
         n2 = 2820 * 365 + 683;
#else
         n1 = 2147;
         n2 = 784175;
#endif
         break;
      case CALENDAR_CHINESE:
         calendar_epoch = CHINESE_CALENDAR_EPOCH + 90;
         n1 = 128;
         n2 = 46751;
         break;
      default:       /* undefined calendar */
         assert( 1);
         return( -1);
      }
   jd -= calendar_epoch;
   day_in_cycle = mod( jd, n2);
   year = n1 * (( jd - day_in_cycle) / n2);
   year += day_in_cycle * n1 / n2;
   return( year);
}

/* day_to_dmy( ) first estimates the year corresponding to an input JD,
and calls get_calendar_data( ) for that year.  Occasionally,  it will
find that the guesstimate was off;  in such cases,  it moves ahead or
back a year and tries again.  Once it's done,  jd - year_ends[0] gives
the number of days since New Years Day;  by subtracting month_data[]
values,  we quickly determine which month and day of month we're in.  */

void DLL_FUNC day_to_dmy( const long jd, int DLLPTR *day,
                  int DLLPTR *month, long DLLPTR *year, const int calendar)
{
   long year_ends[2];
   long curr_jd;
   char month_data[N_MONTHS];
   int calendar_to_use = calendar;

   if( calendar == CALENDAR_JULIAN_GREGORIAN)
      calendar_to_use =
          ((jd > GREGORIAN_SWITCHOVER_JD) ? CALENDAR_GREGORIAN : CALENDAR_JULIAN);

   *year = approx_year( jd, calendar_to_use);
   *day = -1;           /* to signal an error */
   do
      {
      if( get_calendar_data( *year, year_ends, month_data, calendar_to_use))
         return;
      if( year_ends[0] > jd)
         (*year)--;
      if( year_ends[1] <= jd)
         (*year)++;
      }
   while( year_ends[0] > jd || year_ends[1] <= jd);

   curr_jd = year_ends[0];
   *month = -1;
   for( int i = 0; i < N_MONTHS; i++)
      {
      *day = (int)( jd - curr_jd);
      if( *day < (int)month_data[i])
         {
         *month = i + 1;
         (*day)++;
         return;
         }
      curr_jd += (long)month_data[i];
      }
   return;
}
