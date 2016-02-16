/* oblitest.cpp: tests various obliquity formulae

Copyright (C) 2011, Project Pluto

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
#include "lunar.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

/* By default,  this program shows the obliquity from years 0 to 4000 at
200-year steps,  using a variety of obliquity formulae.  Command-line
arguments can change the range/step.      */

double iau_obliquity( const double t_cen);
double file_obliquity( const double t_cen);
double spline_obliquity( const double t_cen);

int main( const int argc, const char **argv)
{
   double t0 = (argc > 1 ? atof( argv[1]) : 0.);
   const double dt = (argc > 2 ? atof( argv[2]) : 200.);
   unsigned n_steps = (argc > 3 ? (unsigned)atoi( argv[3]) : 21);
   const char *header_trailer =
   "    year  :    Laskar      cubic         file      spline  dSpline   dCubic\n";

   printf( "%s", header_trailer);
   while( n_steps--)
      {
      const double t_cen = (t0 - 2000.) / 100.;
      const double laskar_obliquity = mean_obliquity( t_cen);
      const double iau_value        = iau_obliquity( t_cen);
      const double file_value       = file_obliquity( t_cen);
      const double spline_value     = spline_obliquity( t_cen);

      printf( "%10.2f: %.8f %.8f %.8f %.8f %.3f %.3f\n", t0,
               laskar_obliquity * 180. / PI,
               iau_value        * 180. / PI,
               file_value       * 180. / PI,
               spline_value     * 180. / PI,
               (laskar_obliquity - spline_value) * 3600. * 180. / PI,
               (laskar_obliquity - iau_value) * 3600. * 180. / PI);
      t0 += dt;
      }
   printf( "%s", header_trailer);
   if( file_obliquity( 1.) == 0.)
      {
      printf( "'File' obliquity is interpolated from values found in files at\n\n");
      printf( "http://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/63\n\n");
      printf( "Get those files,  and the 'file' obliquities will be computed and shown.\n");
      }
   return( 0);
}
