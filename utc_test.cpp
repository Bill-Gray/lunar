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

/* These macros determine the MJD of the given date in 'YEAR'.         */
/* They're valid for _non-negative_ years in the _Gregorian_ calendar. */
/* Some are not used here,  and therefore commented out,  just to      */
/* avoid compiler warnings.                                            */

#define JAN_1( YEAR) (((YEAR) * 365 + ((YEAR) - 1) / 4 - ((YEAR) - 1) / 100 \
                         + ((YEAR) - 1) / 400) - 678940)
// #define FEB_1( YEAR) (JAN_1( YEAR) + 31)
#define MAR_1( YEAR) (((YEAR)*365 + (YEAR)/4 - (YEAR)/100 + (YEAR)/400) - 678881)
#define APR_1( YEAR) (MAR_1( YEAR) + 31)
#define MAY_1( YEAR) (APR_1( YEAR) + 30)
#define JUN_1( YEAR) (MAY_1( YEAR) + 31)
#define JUL_1( YEAR) (JUN_1( YEAR) + 30)
// #define AUG_1( YEAR) (JUL_1( YEAR) + 31)
// #define SEP_1( YEAR) (AUG_1( YEAR) + 31)
// #define OCT_1( YEAR) (SEP_1( YEAR) + 30)
// #define NOV_1( YEAR) (OCT_1( YEAR) + 31)
// #define DEC_1( YEAR) (NOV_1( YEAR) + 30)

#include <stdio.h>
#include <stdlib.h>
#include "watdefs.h"
#include "afuncs.h"

int main( const int argc, const char **argv)
{
   int year = (argc > 1 ? atoi( argv[1]) : 1970);
   const int end_year = (argc > 2 ? atoi( argv[2]) : 2040);
   unsigned count = 0;

   printf( "Leap seconds for years %d to %d\n", year, end_year);
   printf( "(See the 'official' list at http://maia.usno.navy.mil/ser7/tai-utc.dat)\n");
   printf( "Future leap seconds are predicted using the method described\n");
   printf( "in 'delta_t.cpp',  and come without warranty of any kind.\n");
   for( ; year <= end_year; year++)
      if( year >= 1972)
         for( unsigned pass = 0; pass < 2; pass++)
            {
            const double tai_minus_tdt = -32.184;      /* by definition */
            const double jd = 2400000.5 + (pass ? JUL_1( year) : JAN_1( year));
            const double tai_minus_utc_before =
                       tai_minus_tdt + td_minus_utc( jd - .0001);
            const double tai_minus_utc_after =
                       tai_minus_tdt + td_minus_utc( jd + .0001);

            if( tai_minus_utc_before != tai_minus_utc_after)
               {
               printf( "%d %s : %.3f%s", year, (pass ? "Jul" : "Jan"),
                              tai_minus_utc_after,
                              (count % 3 == 2 ? "\n" : "   "));
               count++;
               }
            }
   printf( "\n");
   return( 0);
}
