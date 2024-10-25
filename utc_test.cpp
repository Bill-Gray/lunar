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
#include "watdefs.h"
#include "afuncs.h"
#include "mjd_defs.h"

int main( const int argc, const char **argv)
{
   int year = 1970, end_year = 2040, i;
   unsigned count = 0;

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'p':
               mjd_end_of_predictive_leap_seconds = atoi( argv[i] + 2);
               printf( "No predicted leap seconds after MJD %d\n",
                           mjd_end_of_predictive_leap_seconds);
               break;
            default:
               printf( "Option '%s' ignored\n", argv[i]);
               break;
            }
      else
         sscanf( argv[i], "%d,%d", &year, &end_year);
   printf( "Leap seconds for years %d to %d\n", year, end_year);
   printf( "(See the 'official' list at https://hpiers.obspm.fr/iers/bul/bulc/UTC-TAI.history)\n");
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
