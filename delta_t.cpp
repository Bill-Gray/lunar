/* delta_t.cpp: computes Delta-T (UT-TD) for time system conversions

Copyright (C) 2015, Project Pluto

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
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include "watdefs.h"
#include "afuncs.h"
#include "mjd_defs.h"

#ifdef __WATCOMC__
#define sinl(x) ((long double)sin((double)x))
#endif

/* 2013 Feb 1:  (BJG) Corrected some comments;  also revised td_minus_ut( )
to check default_delta_t_string only for years before 1620.  This should
have zero effect on the results,  but should save some CPU cycles.   */

/* The following method for computing Delta-T works as follows.  For dates
between 1620 and 2016,  the value comes from the 'delta_t_table' (copied
in turn from Meeus' _Astronomical Algorithms_),  extended using data from
the USNO sites referenced below.  The table gives delta-t at two-year
intervals;  between those values,  a linear interpolation is used.  If
the year is before 1620,  one of four polynomial approximations is used
(see references below).   If the date is after the end of the table,
a linear extrapolation is used,  with a quadratic term added to that
which assumes an acceleration of 32.5 seconds/century^2.

   Updated at two-year intervals to add a new Delta-T table entry from

ftp://maia.usno.navy.mil/ser7/deltat.data
ftp://maia.usno.navy.mil/ser7/deltat.preds

   which is no longer available;  2020 and later uses

ftp://ftp.iers.org/products/eop/rapid/standard/finals.all
*/

static const short delta_t_table[] =
 {12400, 11500, 10600, 9800, 9100, 8500, 7900,  /*  1620-1632 */
   7400, 7000, 6500, 6200, 5800, 5500, 5300,    /*  1634-1646 */
   5000, 4800, 4600, 4400, 4200, 4000, 3700,    /*  1648-1660 */
   3500, 3300, 3100, 2800, 2600, 2400, 2200,    /*  1662-1674 */
   2000, 1800, 1600, 1400, 1300, 1200, 1100,    /*  1676-1688 */
   1000,  900,  900,  900,  900,  900,  900,    /*  1690-1702 */
    900,  900, 1000, 1000, 1000, 1000, 1000,    /*  1704-1716 */
   1100, 1100, 1100, 1100, 1100, 1100, 1100,    /*  1718-1730 */
   1100, 1200, 1200, 1200, 1200, 1200, 1300,    /*  1732-1744 */
   1300, 1300, 1300, 1400, 1400, 1400, 1500,    /*  1746-1758 */
   1500, 1500, 1500, 1600, 1600, 1600, 1600,    /*  1760-1772 */
   1600, 1700, 1700, 1700, 1700, 1700, 1700,    /*  1774-1786 */
   1700, 1700, 1600, 1600, 1500, 1400, 1370,    /*  1788-1800 */
   1310, 1270, 1250, 1250, 1250, 1250, 1250,    /*  1802-1814 */
   1250, 1230, 1200, 1140, 1060,  960,  860,    /*  1816-1828 */
    750,  660,  600,  570,  560,  570,  590,    /*  1830-1842 */
    620,  650,  680,  710,  730,  750,  770,    /*  1844-1856 */
    780,  790,  750,  640,  540,  290,  160,    /*  1858-1870 */
   -100, -270, -360, -470, -540, -520, -550,    /*  1872-1884 */
   -560, -580, -590, -620, -640, -610, -470,    /*  1886-1898 */
   -270,    0,  260,  540,  770, 1050, 1340,    /*  1900-1912 */
   1600, 1820, 2020, 2120, 2240, 2350, 2390,    /*  1914-1926 */
   2430, 2400, 2390, 2390, 2370, 2400, 2430,    /*  1928-1940 */
   2530, 2620, 2730, 2820, 2910, 3000, 3070,    /*  1942-1954 */
   3140, 3220, 3310, 3400, 3500, 3650, 3830,    /*  1956-1968 */
   4020, 4220, 4450, 4650, 4850, 5050, 5220,    /*  1970-1982 */
   5380, 5490, 5580,                            /*  1984-1988 */
  5686,    /* 1990  1 1    56.8554     .3286        -24.6714 */
  5831,    /* 1992  1 1    58.3093    -.1253        -26.1253 */
  5998,    /* 1994  1 1    59.9847     .1993        -27.8007 */
  6163,    /* 1996  1 1    61.6287     .5553        -29.4447 */
  6297,    /* 1998  1 1    62.9659     .2181        -30.7819 */
  6383,    /* 2000  1 1    63.8285268  .3554732     -31.6445268 */
  6430,    /* 2002  1 1    64.2998152 -.1158152     -32.1158152 */
  6457,    /* 2004  1 1    64.5736400 -.3896400     -32.3896400 */
  6485,    /* 2006  1 1    64.8452                              */
  6546,    /* 2008  1 1:   65.4574                              */
  6607,    /* 2010  1 1:   66.0699                              */
  6660,    /* 2012  1 1:   66.6030                              */
  6728,    /* 2014  1 1:   67.2810                              */
  6810,    /* 2016  1 1:   68.1024                              */
  6897,    /* 2018  1 1:   68.9677                              */
  6936,    /* 2020  1 1:   UT1 - UTC = -0.1772554; 69.3612554   */
  6929,    /* 2022  1 1:   UT1 - UTC = -0.1104988; 69.2944988   */
  6918 };  /* 2024  1 1:   UT1 - UTC =  0.0087688; 69.1752068   */

/* 8 Aug 2000:  Some people have expressed an interest in being able to
   insert their own formulae for Delta-T while running Guide.  I've
   enabled this by letting people add a string of coefficients in
   GUIDE.DAT;  gory details of how the string works are at

http://www.projectpluto.com/update7.htm#delta_t_def

   The following functions parse out the string of coefficients and
   let you set which string is used. */

static int evaluate_delta_t_string( const double year, double *delta_t,
                              const char *dt_string)
{
   int rval = -1;

   for( const char *tptr = dt_string; rval && *tptr; tptr++)
      if( tptr == dt_string || tptr[-1] == ';')
         {
         double year1, year2;
         int bytes_read;

         if( sscanf( tptr, "%lf,%lf:%n", &year1, &year2, &bytes_read) == 2
                  && year >= year1 && year < year2)
            {
            double offset = 2000., x, power = 1., coeff;

            tptr += bytes_read;
            if( *tptr == 'o')
               {
               tptr++;
               sscanf( tptr, "%lf,%n", &offset, &bytes_read);
               tptr += bytes_read;
               }
            x = (year - offset) / 100.;
            *delta_t = 0.;
            rval = 0;
            while( tptr[-1] != ';' && *tptr
                     && sscanf( tptr, "%lf%n", &coeff, &bytes_read))
               {
               *delta_t += power * coeff;
               power *= x;
               tptr += bytes_read;
               if( *tptr == ',')
                  tptr++;
               }
            tptr--;
            }
         }
   return( rval);
}

static const char *default_delta_t_string =
"-1e+308,-500:o1820,-20,0,32;\
-500,500:o0,10583.6,-1014.41,33.78311,-5.952053,-.1798452,.022174192,.0090316521;\
500,1600:o1000,1574.2,-556.01,71.23472,.319781,-.8503463,-.005050998,.0083572073;\
1600,1620:o1600,120,-98.08,-153.2,140.272";
// static const char *td_minus_dt_string = default_delta_t_string;
static const char *td_minus_dt_string = NULL;

/* If a user attempts to set a NULL or "blank" Delta-T definition,   */
/* we fall back on the above default_delta_t_string.                 */

void DLL_FUNC reset_td_minus_dt_string( const char *string)
{
   td_minus_dt_string = (string && *string ? string :
                        default_delta_t_string);
}


/* Niels Bohr noted that "prediction is hard,  especially about the future."
This is especially true for predicting Delta-T.  Morrison and Stephenson
arrived at a formula which is intended to match the long-term behavior of
Delta-T,  but it doesn't match current observations,  so it's pretty bad
for near-term results (say,  the next decade or two).  This code assumes
that the rate of change of Delta-T can be determined from the last two
entries in the above 'delta_t_table' array, with a quadratic term of 32.5
seconds/century^2 added to this.  This is likely to provide good
"near-term" results.  I don't know of anything likely to provide good
longer-term results.

   Prior to March 2012,  I used quadratic approximations for Delta-T for
years before 948 and years between 948 to 1620,  from F.R. Stephenson and
M. A. Houlden, "Atlas of Historical Eclipse Maps",  Cambridge University
Press,  England (1986), page x.  These can also be found on page 73 of
Meeus' _Astronomical Algorithms_,  1st edition.  They correspond to the
Delta-T string:

-10000,948:2715.6,573.36,46.5;948,1620:50.6,67.5,22.5

   Since then,  I've switched to quadratic and higher-order polynomials
from the _Five Millennium Catalogue of Solar Eclipses_ by Espenak and
Meeus,  as defined in the above default_delta_t_string.  */

#define J2000 2451545.

double default_td_minus_ut( const double jd)
{
   const double year = 2000. + (jd - J2000) / 365.25;
   double rval;
   /* first convert t from JD to years */

   if( td_minus_dt_string)
      if( !evaluate_delta_t_string( year, &rval, td_minus_dt_string))
         return( rval);
   if( year < 1620. && !evaluate_delta_t_string( year, &rval,
                      default_delta_t_string))
         return( rval);
#ifdef NOW_OBSOLETE_MORRISON_STEPHENSON_FORMULA
   dt = (year - 2000.) / 100.;
   if( year < 948.)
      rval = 2715.6 + dt * (573.36 + 46.5 * dt);
   else if( year < 1620.)
      rval = 50.6 + dt * (67.5 + 22.5 * dt);
   else     /* from 1620 to +infinity */
#endif
      {
      const double index_loc = (year - 1620.) / 2.;
      double dt;
      int index = (int)index_loc;
      const short *tptr;
      const int delta_t_table_size =
            sizeof( delta_t_table) / sizeof( delta_t_table[0]);

      if( index > delta_t_table_size - 2 || jd > 3e+7)   /* running off end of table */
         index = delta_t_table_size - 2;
      dt = index_loc - (double)index;
      tptr = delta_t_table + index;
      rval = (double)tptr[0] + (double)(tptr[1] - tptr[0]) * dt;
      rval /= 100.;
      if( dt > 1.)            /* again, past end of table */
         {
         dt = (dt - 1.) / 50.;    /* cvt to centuries past end,  and add the */
         rval += 32.5 * dt * dt;  /* same 32.5 sec/cy^2 used by Stephenson   */
         }
      }
#ifdef TRY_OMITTING
   if( year < 1620.)       /* apply correction from _Astro Almanac_, */
      {                    /* 1991,  page K8: */
      const double n = -23.8946;   /* corrected lunar secular acceleration */

      rval -= 0.000091 * (n + 26.) * (year - 1955.) * (year - 1955.);
      }
#endif
   return( rval);
}

/* If an Earth Orientation Parameters (EOP) file has been loaded,  and
its time span covers the input JD,  you'll get a very precise Delta-T
value.  (See eop_prec.cpp.)  If that doesn't work out,  it'll fill in
a Delta-T using the above method.   */

double DLL_FUNC td_minus_ut( const double jd)
{
   earth_orientation_params eop;

   get_earth_orientation_params( jd, &eop, 4);
   return( eop.tdt_minus_ut1);
}

/* Some notes for other time systems that I may,  someday,  get around
to implementing:

TDT = TAI + 32.184 seconds = TAI + 0.0003725 days,  exactly

TAI & UTC differ by an integer number of seconds;  after 2012 Jul 1,
for example,  TAI-UTC = 35 seconds,  TDT-UTC = 67.184 seconds.
USNO & IERS bulletins often give TDT-UT1;  from this you can get
UT1-UTC (always between -1 and +1 second)

TCB - TDB = L * (jd - 2443144.5) * seconds_per_day,
   where L=1.550505e-08

TAI-GPS = 19 seconds ("fixed" as of 1980;  at that point,  GPS time was
identical to UTC.  But since then,  leap seconds have been added to
UTC,  so that GPS and UTC have slipped apart.)

Note also that this records only announced leap seconds.  UTC can only
really be determined up to about six months in advance.  Beyond that,
it's indeterminate.

The 'leap_intervals' array will work up to 65535 days past the start of
use of leap seconds in 1972,  or until the January 2151 leap second (if
any).  At that point,  it will have to switch from 16-bit integers to
32 bits.  Fortunately,  the compiler will emit a warning.  I hope my
descendants 136 years hence will regard leap seconds as a quaint
historical aberration,  but one never knows.

Previously,  this function assumed no leap seconds except those which had
been officially announced.  The problem is that this causes UTC to deviate
from (predicted) UT.  Thus,  in computing TD-UTC for future dates (those
beyond the end of the leap second table),  we find out which in which
Gregorian half-year jd_utc falls,  then compute the current prediction for
Delta-T in the middle of that year.  Then

TD-UTC = td_minus_tai + floor( td_minus_ut + .5 - td_minus_tai)

   ...i.e.,  choose the value for TD-UTC for that half-year to match
Delta-T as closely as possible (within half a second at the center of
that half-year).

   Note that eventually,  sometime around 2270,  this can result in two
leap seconds being added at the end of June or December.  As currently
defined,  we would instead start inserting leap seconds at the end of
March or September.  Eventually,  if that wasn't enough,  we'd insert a
leap second at the end of every month.  My only excuse is that I expect
to be dead by 2270.

   You may not like this so-called "solution".  I don't,  either.  I'm
open to ideas,  but doubt there is a better "solution".  If you don't
like it,  set the global variable 'mjd_end_of_predictive_leap_seconds'
to (MJD) 59215 = 2021 Jan 1.  (Note that this can even be set to shut
off leap seconds that have already elapsed,  though I can't come up
with a situation where you'd want to do that.)

   Also,  note that before 1961 Jan 1 = MJD 37400, there was no "real" UTC,
and we assume UTC=UT.  From then until 1972 Jan 1 = MJD 41317, the idea
was that UTC followed atomic clocks,  but at a slightly different rate to
keep UTC mostly in sync with the earth's rotation.  The end result is that
TD-UTC could be expressed as a linear relationship (or a series of 13
linear relationships given below,  each covering a given time span.)

   After 1972 Jan 1,  things switched over to the current system:
TAI-UTC must always be an integer,  with leap seconds inserted at
irregular intervals,  always at the end of June or December,  to keep
UTC within .9 seconds of UT.   */

/* TDB (barycentric dynamical time) and TDT (terrestrial dynamic time)
differ by periodic terms.  The largest has an annual frequency and
amplitude of 1.656 milliseconds;  the rest are so small that they are
often neglected.  However,  I've implemented the first six terms from
1988IAUS..128..419F .  That paper gives 44 terms;  the remaining terms
are below the two microsecond level.

   At present,  I'm ignoring this less-than-two-millisecond difference.

   An earlier version of this code had the following values.  I dunno
whence I got them (failure to document my own code properly!);  I include
them here for reference.

   const double deg2rad = pi / 180.;
   static const double amplitude[4] = {
           1656.675, 22.418, 13.84, 10.216 };
   static const double freq[4] = { 35999.4729 * deg2rad,
            32964.467 * deg2rad, 71998.746 * deg2rad, 3599.373 * deg2rad };
   static const double phase[4] = { 357.5287 * deg2rad,
             246.199 * deg2rad, 355.057 * deg2rad, 243.451 * deg2rad };  */

long double DLL_FUNC tdb_minus_tdt( const long double t_centuries)
{
   static const long double amplitude[6] = { 1656.6894e-6,
            22.4175e-6, 13.8399e-6, 4.7701e-6, 4.6767e-6, 2.2566e-6 };
   static const long double freq[6] = { 628.30758494, 575.33848843,
            1256.61516988, 52.96909651, 606.97767539, 21.32990954 };
   static const long double phase[6] = { 6.2400497, 4.2969771,
            6.1968995, 0.4444038, 4.0211937, 5.5431320 };
   long double rval = 0;
   size_t i;

   for( i = 0; i < sizeof( phase) / sizeof( phase[0]); i++)
      rval += amplitude[i] * sinl( freq[i] * t_centuries + phase[i]);
   return( rval);     /* difference is in _seconds_  */
}

#define utc0  (JAN_1( 1972))
      /*  'utc0' = MJD of date when the UTC leap seconds began */

int mjd_end_of_predictive_leap_seconds = INT_MAX;

double DLL_FUNC td_minus_utc( const double jd_utc)
{
   const double tdt_minus_tai = 32.184;
   const double mjd_utc = jd_utc - 2400000.5;

   if( mjd_utc < (double)utc0)  /* between jan 1961 & jan 1972 */
      {
      static const unsigned short ranges[13] =  { JAN_1( 1961), AUG_1( 1961),
                      JAN_1( 1962), NOV_1( 1963), JAN_1( 1964), APR_1( 1964),
                      SEP_1( 1964), JAN_1( 1965), MAR_1( 1965), JUL_1( 1965),
                      SEP_1( 1965), JAN_1( 1966), FEB_1( 1968) };

      for( int i = 12; i >= 0; i--)
         if( mjd_utc >= (double)ranges[i])
            {
            static const double offset[13] = {
                    1.4228180 - JAN_1( 1961) * 0.0012960,
                    1.3728180 - JAN_1( 1961) * 0.0012960,
                    1.8458580 - JAN_1( 1962) * 0.0011232,
                    1.9458580 - JAN_1( 1962) * 0.0011232,
                    3.2401300 - JAN_1( 1965) * 0.0012960,
                    3.3401300 - JAN_1( 1965) * 0.0012960,
                    3.4401300 - JAN_1( 1965) * 0.0012960,
                    3.5401300 - JAN_1( 1965) * 0.0012960,
                    3.6401300 - JAN_1( 1965) * 0.0012960,
                    3.7401300 - JAN_1( 1965) * 0.0012960,
                    3.8401300 - JAN_1( 1965) * 0.0012960,
                    4.3131700 - JAN_1( 1966) * 0.0025920,
                    4.2131700 - JAN_1( 1966) * 0.0025920  };
            static const short scale[13] =      { 12960, 12960, 11232,
                                    11232, 12960, 12960, 12960, 12960,
                                    12960, 12960, 12960, 25920, 25920 };
            const double tai_minus_utc = offset[i] +
                        mjd_utc * (double)scale[i] * 1.e-7;

            return( tdt_minus_tai + tai_minus_utc);
            }
      }
   else              /* integral leap seconds */
      {
      int imjd_utc = (int)mjd_utc;
      static const unsigned short leap_intervals[] = {
                 JAN_1( 1972) - utc0, JUL_1( 1972) - utc0, JAN_1( 1973) - utc0,
                 JAN_1( 1974) - utc0, JAN_1( 1975) - utc0, JAN_1( 1976) - utc0,
                 JAN_1( 1977) - utc0, JAN_1( 1978) - utc0, JAN_1( 1979) - utc0,
                 JAN_1( 1980) - utc0, JUL_1( 1981) - utc0, JUL_1( 1982) - utc0,
                 JUL_1( 1983) - utc0, JUL_1( 1985) - utc0, JAN_1( 1988) - utc0,
                 JAN_1( 1990) - utc0, JAN_1( 1991) - utc0, JUL_1( 1992) - utc0,
                 JUL_1( 1993) - utc0, JUL_1( 1994) - utc0, JAN_1( 1996) - utc0,
                 JUL_1( 1997) - utc0, JAN_1( 1999) - utc0, JAN_1( 2006) - utc0,
                 JAN_1( 2009) - utc0, JUL_1( 2012) - utc0, JUL_1( 2015) - utc0,
                 JAN_1( 2017) - utc0 };
      const int n_leap_seconds = sizeof( leap_intervals) / sizeof( leap_intervals[0]);

      if( imjd_utc > mjd_end_of_predictive_leap_seconds)
         imjd_utc = mjd_end_of_predictive_leap_seconds;
      if( imjd_utc >= DEC_1( 2021))
         {
         int day = imjd_utc + 2400000 - 1721058;
         int year = (int)( (int64_t)day * (int64_t)400 / (int64_t)146097);
         int low, high, july_1;

         low = (int)JAN_1( year);  /* The above value for 'year' is correct */
         if( imjd_utc < low)     /* more than 99% of the time.  But we    */
            {                    /* may find,  for 31 December,  that     */
            year--;              /* it's too high by one year.            */
            high = low;
            low = (int)JAN_1( year);
            }
         else
            high = (int)JAN_1( year + 1);
                      /*  jul  aug  sep  oct  nov  dec.. jul 1 is exactly 184 */
         july_1 = high - (31 + 31 + 30 + 31 + 30 + 31); /* days before jan 1  */
         if( imjd_utc < july_1) /* in first half of the year */
            high = july_1;
         else                /* in second half of the year */
            low = july_1;
         return( tdt_minus_tai + floor( td_minus_ut( 2400000.5 +
                        (double)( low + high) * .5) + .5 - tdt_minus_tai));
         }
      for( int i = n_leap_seconds - 1; i >= 0; i--)
         if( imjd_utc - utc0 >= (int)leap_intervals[i])
            return( (double)(i + 10) + tdt_minus_tai);
      }
                     /* still here?  Must be before jan 1961,  so UTC = UT1: */
   return( td_minus_ut( jd_utc));
}

double DLL_FUNC tdb_minus_utc( const double jd_utc)
{
   const double t_cen = (jd_utc - J2000) / 36525.;

   return( tdb_minus_tdt( t_cen) / seconds_per_day + td_minus_utc( jd_utc));
}
