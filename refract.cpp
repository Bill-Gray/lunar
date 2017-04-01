/* refract.cpp: functions & test code for low-precision refraction

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

#include <math.h>
#include "watdefs.h"
#include "afuncs.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

/*
   REFRACT.CPP contains functions to convert an observed altitude
(one affected by refraction) to a true altitude (one that would be
seen on an airless planet),  or vice versa.  The first case is
handled using the method of G G Bennett,  'The Calculation of
Astronomical Refraction in Marine Navigation',  _Journal of the
Institute for Navigation_, Vol 35, p 255-259 (1982),  as explained
in Meeus' _Astronomical Algorithms_,  p 102.  The maximum error of
this is stated to be .9 arcsecond,  for the range 0-90 degrees.

   For the inverse problem,  handled in the function reverse_refraction(),
we first get a "pretty good" answer using the formula of Saemundsson,
_Sky & Telescope_,  p 70,  July 1986,  again as explained in _AA_.
The only problem is that this formula is accurate to within about 4".
One iteration using the refraction( ) function repairs this problem.

   All inputs and outputs from these functions are in radians.
*/

double DLL_FUNC refraction( const double observed_alt)
{
   double rval, ang = observed_alt;

   ang += (7.31 * PI / 180.) / (ang * 180. / PI + 4.4);
   rval = cos( ang) / sin( ang);       /* from Meeus,  _AA_,  p 102 */
   rval -= .06 * sin( (14.7 * rval + 13.) * PI / 180.);
   return( rval * (PI / 180.) / 60.);        /* cvt to radians */
}

double DLL_FUNC reverse_refraction( const double true_alt)
{
   double rval, delta = 1.;
   const double tolerance = .01 * (PI / 180.) / 3600.;
   int n_iter = 10;

   rval = true_alt + (10.3 * PI / 180.) / (true_alt * 180. / PI + 5.11);
   rval = cos( rval) / sin( rval);  /* from Meeus,  _AA_,  p 102 */
   rval *= 1.02 * (PI / 180.) / 60.;        /* cvt to radians */
                  /* The above gives a good first approximation */
                  /* Now improve that answer as follows: */
   while( n_iter-- && (delta > tolerance || delta < -tolerance))
      {
      delta = rval;
      rval = refraction( true_alt + rval);
      delta -= rval;
      }
   return( rval);
}

double DLL_FUNC saasta_refraction( const double observed_alt,
         const double pressure_mb, const double temp_kelvin,
         const double relative_humidity)
{
   double xi;
   const double tan_z0 = cos( observed_alt) / sin( observed_alt);
   const double tan_z0_2 = tan_z0 * tan_z0;
   const double delta = 18.36;
   const double pw0 =
             relative_humidity * exp( delta * log( temp_kelvin / 247.1));
   const double q = (pressure_mb - .156 * pw0) / temp_kelvin;

   xi = 16.271 * q * tan_z0 * (1. + .0000394 * q * tan_z0_2) -
               .0000749 * pressure_mb * tan_z0 * (1. + tan_z0_2);
         /* The above refraction is in _arcseconds_... */
   return( xi * (PI / 180.) / 3600.);        /* cvt to radians */
}

double DLL_FUNC reverse_saasta_refraction( const double true_alt,
         const double pressure_mb, const double temp_kelvin,
         const double relative_humidity)
{
   double rval;

   rval = true_alt + (10.3 * PI / 180.) / (true_alt * 180. / PI + 5.11);
   rval = cos( rval) / sin( rval);  /* from Meeus,  _AA_,  p 102 */
   rval *= 1.02 * (PI / 180.) / 60.;        /* cvt to radians */
                  /* The above gives a good first approximation */
                  /* Now improve that answer as follows: */
   rval = saasta_refraction( true_alt + rval,
                        pressure_mb, temp_kelvin, relative_humidity);
   return( rval);
}

#ifdef TEST_PROGRAM
#include <stdio.h>
#include <stdlib.h>

int main( const int argc, const char **argv)
{
   int i;
   double ang, ref, scale = 1.;

   if( argc == 3)
      scale = atof( argv[2]);
   for( i = 1; i <= 90; i++)
      {
      ang = (double)i * PI / 180;
      printf( "%2d:%9.5lf %9.5lf    %9.5lf %9.5lf\n", i,
                  refraction( ang) * 60. * 180. / PI,
                  reverse_refraction( ang) * 60. * 180. / PI,
                  saasta_refraction( ang, 1013., 293., .20 ) * 60. * 180. / PI,
                  reverse_saasta_refraction( ang, 1013., 293., .20) * 60. * 180. / PI);
      }
   printf( "Polar refraction: %.4lf arcsec\n",
                  refraction( PI / 2) * 3600. * 180. / PI);
   if( argc > 1)
      {
      ang = atof( argv[1]) * PI / 180.;
      ref = refraction( ang);
      printf( "M:   %9.5lf %9.5lf %9.5lf\n", ref * 60. * 180. / PI,
                  reverse_refraction( ang - ref) * 60. * 180. / PI,
                  reverse_refraction( ang) * 60. * 180. / PI);
      ref = saasta_refraction( ang, 1013., 293., .20);
      printf( "S:   %9.5lf %9.5lf %9.5lf\n", ref * 60. * 180. / PI,
             reverse_saasta_refraction( ang - ref, 1013., 293., .20)
                                                    * 60. * 180. / PI,
             reverse_saasta_refraction( ang, 1013., 293., .20)
                                                    * 60. * 180. / PI);
      ref = 0.;
      for( i = 0; i < 9; i++)
         {
         ref = refraction( ang + ref) * scale;
         printf( "     %9.5lf %9.5lf\n", ref * 60. * 180. / PI,
                                         ref * 180. / PI);
         }
      }
   return( 0);
}
#endif
