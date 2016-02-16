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

   I also noticed that the above series will diverge for s > 1/sqrt(|y|),
and as you approach that limit,  convergence will be absurdly slow.
And there's no way around _that_ convergence issue.  Which is why I don't
see this method as being actually useful;  it's slow and doesn't always
work.  Maybe it was once useful,  long ago (pre-1985 or so)  on machines
without built-in math coprocessors to handle trig functions?         */

static double paul_schlyter_soln( const double e, const double q, const double dt)
{
   const double k = 0.01720209895;
   const double a = 0.75 * dt * k * sqrt( (1 + e) / (q*q*q) );
   const double b = sqrt( 1 + a*a );
   const double W = cbrt(b + a) - cbrt(b - a);
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
   const double q1 = k * sqrt( (1 + ecc) / q) / (2. * q);
   const double gamma = (1. - ecc) / (1. + ecc);
   const double q2 = q1 * t;
   const double max_s = 1. / sqrt( fabs( gamma));
   double true_anomaly, radius, prev_s;
   double s = 2. / (3. * fabs( q2));
   double pauls_v = paul_schlyter_soln( ecc, q, t);

   printf( "Paul Schlyter soln: %lf\n", (180. / pi) * pauls_v);
   if( argc == 4)
      {
      s = 2. / tan( 2. * atan( exp( log( tan( atan( s) / 2.)) / 3.)));
      if( t < 0)
         s = -s;
      }
   else
      s = tan( pauls_v / 2.);
   printf( "Max s = %lf; initial s %lf\n", max_s, s);
   if( s > max_s * .999)
      s = max_s * .999;
   if( ecc != 1.)
      do
         {
         double term = 1., s_power = -s * s * s, f, f_prime;
         const double multiplier = -s * s * gamma;
         int i;
         const int max_iter = 10000;

         prev_s = s;
         f = q1 * t - s;
         f_prime = -1.;
         for( i = 1; i < max_iter && fabs( term) > tolerance; i++)
            {
            term = ((double)i - (double)(i + 1) * gamma) * s_power;
            f_prime += term / s;
            term /= (double)( i + i + 1);
            f += term;
            s_power *= multiplier;
            }
         if( i == max_iter)
            f *= .3;
         s -= f / f_prime;
         printf( "prev_s = %lf, s = %lf, f = %lf, i = %d\n", prev_s, s, f, i);
         if( i == max_iter)
            printf( "Last term = %lf\n", term);
         if( s > (max_s + prev_s) / 2.)
            s = (max_s + prev_s) / 2.;
         }
      while( fabs( prev_s - s) > tolerance);
   true_anomaly = 2. * atan( s);
   radius = q * (1. + ecc) / (1. + ecc * cos( true_anomaly));
   printf( "True anomaly: %lf\n", true_anomaly * 180. / pi);
   printf( "Radius vector (AU): %lf\n", radius);
   return( 0);
}
