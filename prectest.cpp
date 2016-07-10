/* prectest.cpp: demonstrates/shows results from precession funcs

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

#include <stdio.h>
#include <stdlib.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"         /* for obliquity( ) prototype */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int main( const int argc, const char **argv)
{
   double t1, t2, matrix[9];
   double p[2];
   int i;

   t1 = atof( argv[1]);
   t2 = atof( argv[2]);
   if( argc > 3)
      {
      p[0] = atof( argv[3]) * PI / 180.;
      p[1] = atof( argv[4]) * PI / 180.;
      }
   if( argc < 6)
      setup_precession( matrix, t1, t2);
   else
      setup_ecliptic_precession( matrix, t1, t2);
   for( i = 0; i < 9; i++)            /* print out the precession matrix */
      printf( "%15.11lf%s", matrix[i], (i % 3 == 2) ? "\n" : " ");
   if( argc > 3)        /* precess an example value;  then reverse */
      {                 /* the precession to see if we recover     */
      precess_ra_dec( matrix, p, p, 0);
      printf( "%lf %lf\n", p[0] * 180. / PI, p[1] * 180. / PI);
      precess_ra_dec( matrix, p, p, 1);
      printf( "%lf %lf\n", p[0] * 180. / PI, p[1] * 180. / PI);
      }
}
