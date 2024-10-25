#ifndef MJD_DEFS_H
#define MJD_DEFS_H

#include <stdint.h>
#include <limits.h>

/* These macros determine the MJD of the start of a month in the given
'YEAR' at compile time.  If sizeof( long) == 8,  they are valid for
(Gregorian) years after roughly 20 billion years ago. This was chosen
to comfortably include all dates since the Big Bang.  If sizeof( long)
== 4,  they only work between years -5877630 and +5881469;  beyond
that,  the result overflows the 32-bit long.

   February 1 will always be 31 days after January 1.  March 1 through
December 1 will always be a fixed number of days before January 1 of
the following year.   */

#if LONG_MAX == INT64_MAX
   #define BASE_YEAR 19999999999L
#elif LONG_MAX == INT32_MAX
   #define BASE_YEAR 1999999999L
#else
   #error "Long integers have an unrecognized size"
#endif
#define JAN_1( YEAR) (((YEAR) * 365L + ((YEAR) + BASE_YEAR) / 4L - ((YEAR) + BASE_YEAR) / 100L \
                         + ((YEAR) + BASE_YEAR) / 400L) - 678940L \
                         - ((BASE_YEAR + 1L) / 400L) * 97L)
#define FEB_1( YEAR) (JAN_1( YEAR) + 31)
#define DEC_1( YEAR) (JAN_1( (YEAR)+1) - 31)
#define NOV_1( YEAR) (DEC_1( YEAR) - 30)
#define OCT_1( YEAR) (NOV_1( YEAR) - 31)
#define SEP_1( YEAR) (OCT_1( YEAR) - 30)
#define AUG_1( YEAR) (SEP_1( YEAR) - 31)
#define JUL_1( YEAR) (AUG_1( YEAR) - 31)
#define JUN_1( YEAR) (JUL_1( YEAR) - 30)
#define MAY_1( YEAR) (JUN_1( YEAR) - 31)
#define APR_1( YEAR) (MAY_1( YEAR) - 30)
#define MAR_1( YEAR) (APR_1( YEAR) - 31)

#endif   /* #ifdef MJD_DEFS_H */
