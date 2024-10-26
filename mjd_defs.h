#ifndef MJD_DEFS_H
#define MJD_DEFS_H

/* These macros determine the MJD of the start of a month in the given
(Gregorian) 'YEAR' at compile time.  Separate macros are needed for
positive and non-negative years.

   February 1 will always be 31 days after January 1.  March 1 through
December 1 will always be a fixed number of days before January 1 of
the following year.   */

#define JAN_1_POS( YEAR) ((YEAR) * 365L + ((YEAR) - 1) / 4L - ((YEAR) - 1) / 100L \
                                        + ((YEAR) - 1) / 400L - 678940L)
#define JAN_1_NEG( YEAR) ((YEAR) * 365L + (YEAR) / 4L - (YEAR) / 100L \
                                        + (YEAR) / 400L - 678940L - 1L)
#define JAN_1( YEAR)  ((YEAR) > 0 ? JAN_1_POS( YEAR) : JAN_1_NEG( YEAR))

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
