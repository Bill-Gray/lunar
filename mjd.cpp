#include <stdio.h>
#include <stdlib.h>

#define JAN_1( YEAR) ((YEAR * 365 + (YEAR - 1) / 4 - (YEAR - 1) / 100 + (YEAR - 1) / 400) - 678940)
#define JUL_1( YEAR) ((YEAR * 365 + YEAR / 4 - YEAR / 100 + YEAR /400) - 678759)

#define N_LEAP_SECONDS 28

      static const unsigned short leap_intervals[N_LEAP_SECONDS] = {
                 JAN_1( 1972), JUL_1( 1972), JAN_1( 1973),
                 JAN_1( 1974), JAN_1( 1975), JAN_1( 1976),
                 JAN_1( 1977), JAN_1( 1978), JAN_1( 1979),
                 JAN_1( 1980), JUL_1( 1981), JUL_1( 1982),
                 JUL_1( 1983), JUL_1( 1985), JAN_1( 1988),
                 JAN_1( 1990), JAN_1( 1991), JUL_1( 1992),
                 JUL_1( 1993), JUL_1( 1994), JAN_1( 1996),
                 JUL_1( 1997), JAN_1( 1999), JAN_1( 2006),
                 JAN_1( 2009), JUL_1( 2012),
                 JAN_1( 2016),    /* PREDICTED: not an 'official' leap sec */
                 JAN_1( 2019) };  /* PREDICTED: not an 'official' leap sec */

      static const unsigned short old_intervals[N_LEAP_SECONDS] = {
                 41317, 41499, 41683,    /* jan 1972, jul 1972, jan 1973 */
                 42048, 42413, 42778,    /* jan 1974, jan 1975, jan 1976 */
                 43144, 43509, 43874,    /* jan 1977, jan 1978, jan 1979 */
                 44239, 44786, 45151,    /* jan 1980, jul 1981, jul 1982 */
                 45516, 46247, 47161,    /* jul 1983, jul 1985, jan 1988 */
                 47892, 48257, 48804,    /* jan 1990, jan 1991, jul 1992 */
                 49169, 49534, 50083,    /* jul 1993, jul 1994, jan 1996 */
                 50630, 51179, 53736,    /* jul 1997, jan 1998, jan 2006 */
                 54832, 56109,           /* jan 2009, jul 2012           */
                 57388,                  /* jan 2016:  PREDICTED */
                 58484 };                /* jan 2019:  PREDICTED */

int main( const int argc, const char **argv)
{
   int year = atoi( argv[1]), i;

   printf("Year  Jan1   Jul1\n");
   for( i = 0; i < 30; i++, year++)
      printf( "%d: %d %d\n", year,
            JAN_1( year), JUL_1( year));
   for( i = 0; i < N_LEAP_SECONDS; i++)
      if( old_intervals[i] != leap_intervals[i])
         printf( "Mismatch %d\n", i);
   return( 0);
}
