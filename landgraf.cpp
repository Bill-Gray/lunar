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
#include <math.h>

/* Heavily modified version of method from chapter 35 of Meeus' _Astronomical
Algorithms_,  "Near-Parabolic Motion",  which describes a method due to
Werner Landgraf.  WARNING:  this was of intellectual interest.  As will
be described,  I don't see this method as being practically useful,  even
with the enhancements I've made.

   Specifically,  I modified the way in which (35.1) is solved.  That
equation,  rearranged slightly and using 'y' for 'gamma', reads

f(s) = Qt - s - (1-2y)s^3/3 + y(2-3y)s^5/5 - y^2(3-4y)s^7/7 + ....

   ...and we need to do a root-finding such that f(s) = 0.   Meeus
suggests a pretty simple root-finder that has some convergence issues.
In the following,  Newton-Raphson is used.  Fortunately, it's quite easy
to compute f(s) and f'(s) at the same time.

   Theoretically,  the result will always converge for all hyperbolic
cases and for all elliptical cases until you reach the ends of the
semimajor axes (equivalently,  the eccentric anomaly is +/- 90 degrees),
which occurs at 'max_t',  or if s > 1/sqrt(y).  In practice,  a tremendous
number of iterations is required as you approach that point.         */

/* Paul Schlyter's approximation for near-parabolic orbits,  from

http://stjarnhimlen.se/comp/ppcomp.html#19

revised to do one cube root instead of two.  While it doesn't always
provide a good "end solution",  it does seem to always provide a good
starting value for the Landgraf iteration used below.  Though since it
only eliminates one or two iterations compared to the simpler parabolic
starting guess,  it's probably not really worth it.  */

static double paul_schlyter_soln( const double e, const double q, const double dt)
{
   const double k = 0.01720209895;
   const double a = 0.75 * dt * k * sqrt( (1 + e) / (q*q*q) );
   const double tval = cbrt( sqrt( 1 + a*a ) + a);
   const double W = tval - 1. / tval;
   const double W2 = W * W;
   const double f = (1 - e) / (1 + e);

   const double a1 = (2/3) + .4 * W2;
   const double a2 = 7./5. + W2 * (33./35. + 37./175. * W2);
   const double a3 = W2 * ( (432/175) + W2 * ((956/1125) + (84/1575) * W2));

   const double C = W2 / (1 + W2);
   const double g = f * C*C;
   const double w = W * ( 1 + f * C * ( a1 + a2*g + a3*g*g ) );

   const double v = 2 * atan(w);
//    r = q * ( 1 + w*w ) / ( 1 + w*w * f )

   return( v);
}

int main( const int argc, const char **argv)
{
   const double k = 0.01720209895;
   const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923;
   const double tolerance = 1.e-9;
   const double ecc = atof( argv[1]);
   const double q = atof( argv[2]);
   const double t = atof( argv[3]);
   const double Qt = fabs( t) * k * sqrt( (1 + ecc) / q) / (2. * q);
   const double gamma = (1. - ecc) / (1. + ecc);
   const double max_s = 1. / sqrt( fabs( gamma));
   const double g = 1.5 * Qt;         /* (34.6) in Meeus */
   const double y = cbrt( g + sqrt( g * g + 1.));
   double true_anomaly, radius, prev_s;
   double s = y - 1. / y;
   const double pauls_v = paul_schlyter_soln( ecc, q, t);

   if( ecc < 1.)
      {
      const double a = q / (1. - ecc);
      const double t0 =  a * sqrt( a) / k;

      printf( "max_t = %f\n", t0 * (pi / 2. - ecc));
      }

   printf( "Paul Schlyter true anomaly: %f\n", (180. / pi) * pauls_v);
   if( argc > 4)
      {
      s = tan( pauls_v / 2.);
      printf( "Using Paul Schlyter's initial value\n");
      }
   printf( "Max s = %f; initial s %g\n", max_s, s);
   if( s > max_s * .999)
      s = max_s * .999;
   if( s)
      do
         {
         double s_power = -s * s * s, f = Qt - s, f_prime = -1.;
         const double multiplier = -s * s * gamma;
         int i;
         const int max_iter = 1000000;
         bool keep_iterating = true;

         prev_s = s;
         for( i = 1; i < max_iter && keep_iterating; i++)
            {
            double term = ((double)i - (double)(i + 1) * gamma) * s_power;

            f_prime += term / s;
            term /= (double)( i + i + 1);
            f += term;
            s_power *= multiplier;
            if( fabs( term) < tolerance)
               keep_iterating = false;
            if( i > 3 && fabs( term / f) < .001)
               keep_iterating = false;
            }
         s -= f / f_prime;
         printf( "prev_s = %g, s = %g, f = %g, f_prime = %g, i = %d\n",
                        prev_s, s, f, f_prime, i);
         if( s > (8. * max_s + prev_s) / 9.)
            s = (8. * max_s + prev_s) / 9.;
         }
      while( fabs( prev_s - s) > tolerance);
   if( t < 0.)
      s = -s;
   true_anomaly = 2. * atan( s);
   radius = q * (1. + ecc) / (1. + ecc * cos( true_anomaly));
   printf( "True anomaly: %f\n", true_anomaly * 180. / pi);
   printf( "Radius vector (AU): %f\n", radius);
   return( 0);
}
