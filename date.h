/* date.h: header files for time/date functions
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

#ifdef __cplusplus
extern "C" {
#endif

long DLL_FUNC dmy_to_day( const int day, const int month, const long year,
                            const int calendar);
void DLL_FUNC day_to_dmy( const long jd, int DLLPTR *day,
                  int DLLPTR *month, long DLLPTR *year, const int calendar);
void DLL_FUNC full_ctime( char *buff, double jd, const int format);
void DLL_FUNC full_ctimel( char *buff, long double t2k, const int format);
const char * DLL_FUNC set_month_name( const int month, const char *new_name);
const char * DLL_FUNC set_day_of_week_name( const int day_of_week,
                                            const char *new_name);
double DLL_FUNC decimal_day_to_dmy( double jd, long *year, int *month,
                                 const int calendar);
double DLL_FUNC split_time( double jd, long *year, int *month, int *day,
                                 int *hr, int *min, const int time_format);
long double DLL_FUNC split_timel( long double t2k, long *year, int *month, int *day,
                                 int *hr, int *min, const int time_format);
double DLL_FUNC get_time_from_string( double initial_jd,
         const char *time_str, const int time_format, int *is_ut);
long double DLL_FUNC get_time_from_stringl( long double initial_t2k,
         const char *time_str, const int time_format, int *is_ut);
int DLL_FUNC days_in_month( const int month, const long year,
                            const int calendar);

#define FULL_CTIME_FORMAT_MASK           0x700
#define FULL_CTIME_FORMAT_SECONDS        0x000
#define FULL_CTIME_FORMAT_HH_MM          0x100
#define FULL_CTIME_FORMAT_HH             0x200
#define FULL_CTIME_FORMAT_DAY            0x300
#define FULL_CTIME_FORMAT_YEAR           0x400
#define FULL_CTIME_FORMAT_JD             0x500
#define FULL_CTIME_FORMAT_MJD            0x600
#define FULL_CTIME_YEAR_FIRST            0x800
#define FULL_CTIME_YEAR_LAST                 0
#define FULL_CTIME_MONTH_DAY            0x1000
#define FULL_CTIME_DAY_MONTH                 0
#define FULL_CTIME_MONTHS_AS_DIGITS     0x2000
#define FULL_CTIME_TIME_ONLY            0x4000
#define FULL_CTIME_LEADING_ZEROES       0x8000
#define FULL_CTIME_TWO_DIGIT_YEAR       0x10000
#define FULL_CTIME_NO_YEAR              0x20000
#define FULL_CTIME_DAY_OF_WEEK_FIRST    0x40000
#define FULL_CTIME_DAY_OF_WEEK_LAST     0x80000
#define FULL_CTIME_NO_SPACES            0x100000
#define FULL_CTIME_ROUNDING             0x200000
#define FULL_CTIME_NO_COLONS            0x400000
#define FULL_CTIME_DAY_OF_YEAR          0x800000

         /* Some convenience combos of the above flags: */
#define FULL_CTIME_YMD           (FULL_CTIME_YEAR_FIRST | FULL_CTIME_MONTH_DAY)
#define FULL_CTIME_YDM           (FULL_CTIME_YEAR_FIRST | FULL_CTIME_DAY_MONTH)
#define FULL_CTIME_DMY           (FULL_CTIME_YEAR_LAST | FULL_CTIME_DAY_MONTH)
#define FULL_CTIME_MDY           (FULL_CTIME_YEAR_LAST | FULL_CTIME_MONTH_DAY)

         /* "date only" = show in days with zero decimal places: */
#define FULL_CTIME_DATE_ONLY   FULL_CTIME_FORMAT_DAY

#define FULL_CTIME_TENTHS_SEC           0x10
#define FULL_CTIME_WHOLE_SECONDS           0
#define FULL_CTIME_HUNDREDTH_SEC        0x20
#define FULL_CTIME_MILLISECS            0x30

   /* The above four lines are accurate,  but the meaning of those bits is */
   /* now more general:  bits 4-7 describe the number of decimal places.   */
   /* Thus,  0x30 could mean milliseconds,  or milliminutes,  or millidays. */
   /* And up to fifteen decimal places can be shown... though I'd expect    */
   /* serious loss of precision with 64-bit floating-point math!            */

#define FULL_CTIME_N_PLACES( n)            ((n) << 4)

#define FULL_CTIME_0_PLACES             FULL_CTIME_N_PLACES( 0)
#define FULL_CTIME_1_PLACE              FULL_CTIME_N_PLACES( 1)
#define FULL_CTIME_2_PLACES             FULL_CTIME_N_PLACES( 2)
#define FULL_CTIME_3_PLACES             FULL_CTIME_N_PLACES( 3)
#define FULL_CTIME_4_PLACES             FULL_CTIME_N_PLACES( 4)
#define FULL_CTIME_5_PLACES             FULL_CTIME_N_PLACES( 5)
#define FULL_CTIME_6_PLACES             FULL_CTIME_N_PLACES( 6)
#define FULL_CTIME_7_PLACES             FULL_CTIME_N_PLACES( 7)
#define FULL_CTIME_8_PLACES             FULL_CTIME_N_PLACES( 8)
#define FULL_CTIME_9_PLACES             FULL_CTIME_N_PLACES( 9)
#define FULL_CTIME_10_PLACES            FULL_CTIME_N_PLACES( 10)
#define FULL_CTIME_11_PLACES            FULL_CTIME_N_PLACES( 11)
#define FULL_CTIME_12_PLACES            FULL_CTIME_N_PLACES( 12)
#define FULL_CTIME_13_PLACES            FULL_CTIME_N_PLACES( 13)
#define FULL_CTIME_14_PLACES            FULL_CTIME_N_PLACES( 14)
#define FULL_CTIME_15_PLACES            FULL_CTIME_N_PLACES( 15)

/* Some convenience macros combining units and defined number of decimal */
/* places... also given just as examples of how to mix-and-match your own: */

#define FULL_CTIME_DECIMINUTES      (FULL_CTIME_FORMAT_HH_MM | FULL_CTIME_1_PLACE)
#define FULL_CTIME_CENTIMINUTES     (FULL_CTIME_FORMAT_HH_MM | FULL_CTIME_2_PLACES)
#define FULL_CTIME_MILLIMINUTES     (FULL_CTIME_FORMAT_HH_MM | FULL_CTIME_3_PLACES)
#define FULL_CTIME_MICROMINUTES     (FULL_CTIME_FORMAT_HH_MM | FULL_CTIME_6_PLACES)

#define FULL_CTIME_DECIHOURS      (FULL_CTIME_FORMAT_HH | FULL_CTIME_1_PLACE)
#define FULL_CTIME_CENTIHOURS     (FULL_CTIME_FORMAT_HH | FULL_CTIME_2_PLACES)
#define FULL_CTIME_MILLIHOURS     (FULL_CTIME_FORMAT_HH | FULL_CTIME_3_PLACES)
#define FULL_CTIME_MICROHOURS     (FULL_CTIME_FORMAT_HH | FULL_CTIME_6_PLACES)

#define FULL_CTIME_DECIDAYS      (FULL_CTIME_FORMAT_DAY | FULL_CTIME_1_PLACE)
#define FULL_CTIME_CENTIDAYS     (FULL_CTIME_FORMAT_DAY | FULL_CTIME_2_PLACES)
#define FULL_CTIME_MILLIDAYS     (FULL_CTIME_FORMAT_DAY | FULL_CTIME_3_PLACES)
#define FULL_CTIME_MICRODAYS     (FULL_CTIME_FORMAT_DAY | FULL_CTIME_6_PLACES)
#define FULL_CTIME_NANODAYS      (FULL_CTIME_FORMAT_DAY | FULL_CTIME_9_PLACES)
#define FULL_CTIME_PICODAYS      (FULL_CTIME_FORMAT_DAY | FULL_CTIME_12_PLACES)
#define FULL_CTIME_FEMTODAYS     (FULL_CTIME_FORMAT_DAY | FULL_CTIME_15_PLACES)

/* The calendar is stored in the lowest four bits (nibble),  allowing
for 16 possible calendars;  only 9 are currently assigned,  but there
are plenty of candidates for the remaining seven: */

#define CALENDAR_MASK                  0xf

#define CALENDAR_GREGORIAN             0
#define CALENDAR_JULIAN                1
#define CALENDAR_HEBREW                2
#define CALENDAR_ISLAMIC               3
#define CALENDAR_REVOLUTIONARY         4
#define CALENDAR_PERSIAN               5
#define CALENDAR_CHINESE               7
#define CALENDAR_MODERN_PERSIAN        8

/* The "Julian/Gregorian" calendar uses the Julian system before
   mid-October 1582,  and Gregorian after that.   */

#define CALENDAR_JULIAN_GREGORIAN      6

/* The usual switchover date between the two calendars is JD 2299160.5,   */
/* which corresponds to 15 Oct 1582 (Gregorian) = 5 Oct 1582 (Julian).    */
/* That is to say,  4 October 1582 was followed by 15 October 1582;       */
/* 5 October-14 October 1582 "never happened".                            */

/*           October 1582           */
/*       Su Mo Tu We Th Fr Sa       */
/*           1  2  3  4 15 16       */
/*       17 18 19 20 21 22 23       */
/*       24 25 26 27 28 29 30       */
/*       31                         */

#define GREGORIAN_SWITCHOVER_JD     2299160L

         /* The Chinese calendar routines use a buffer of precomputed data */
         /* indicating when various months occur,  which are considered    */
         /* intercalary,  etc.  This array can come from 'chinese.h' and be */
         /* statically included in the library,  or can be omitted and only */
         /* set when desired using the following function.  See jd.cpp      */
         /* for an example of the use of both of the following functions. */
void DLL_FUNC set_chinese_calendar_data( const void *cdata);
int DLL_FUNC get_chinese_intercalary_month( void);

#ifdef __cplusplus
}
#endif
