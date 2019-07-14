/* prectes2.cpp: compares Earth precession done in ecliptic and equatorial
systems.  See 'precess.cpp' for details as to why both methods matter
and when they ought to be used.

Copyright (C) 2017, Project Pluto

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
#include <string.h>
#include <stdio.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"
#include "lunar.h"         /* for obliquity( ) prototype */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int DLL_FUNC setup_equatorial_precession_from_j2000( double DLLPTR *matrix,
                           const double year);

#include <stdio.h>
#include <stdlib.h>

int main( const int argc, const char **argv)
{
   double t1, matrix[9];
   int i, pass;
   const double J2000 = 2451545.;

   t1 = get_time_from_string( 0., (argc == 1 ? "now" : argv[1]),
                                                FULL_CTIME_YMD, NULL);
   printf( "JD %f =", t1);
   t1 = (t1 - J2000) / 365.25 + 2000.;
   printf( " year %f\n", t1);
   printf( "Last digits represent about 0.64 mm on earth's surface\n");
   for( pass = 0; pass < 3; pass++)
      {
      if( !pass)
         {
         printf( "New,  ecliptic-based precession :\n");
         setup_precession( matrix, 2000., t1);
         }
      else if( pass == 1)
         {
         printf( "Original,  rougher precession :\n");
         setup_equatorial_precession_from_j2000( matrix, t1);
         }
      else
         {
         double mat2[9];

         printf( "Differences are (old - new) :\n");
         setup_equatorial_precession_from_j2000( matrix, t1);
         setup_precession( mat2, 2000., t1);
         for( i = 0; i < 9; i++)
            matrix[i] -= mat2[i];
         }
      for( i = 0; i < 9; i++)
         printf( "%14.10f%s", matrix[i], (i % 3 == 2) ? "\n" : " ");
      }
}
