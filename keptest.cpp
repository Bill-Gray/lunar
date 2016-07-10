/* keptest.cpp: test code for Kepler-solving functions

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
#include <stdio.h>
#include <stdlib.h>

#define THRESH 1.e-8
#define PI 3.141592653589793238462643383279502884197169399375105
#define CUBE_ROOT( X)  (exp( log( X) / 3.))


/* 'asinh' = 'arc-hyperbolic sine.'  Most compilers now implement this. */

#ifdef _MSC_VER
static double asinh( const double z)
{
   return( log( z + sqrt( z * z + 1.)));
}
#endif

/* If the eccentricity is very close to parabolic,  and the eccentric
anomaly is quite low,  you can get an unfortunate situation where
roundoff error keeps you from converging.  Consider the just-barely-
elliptical case,  where in Kepler's equation,

M = E - e sin( E)

   E and e sin( E) can be almost identical quantities.  To
around this,  near_parabolic( ) computes E - e sin( E) by expanding
the sine function as a power series:

E - e sin( E) = E - e( E - E^3/3! + E^5/5! - ...)
= (1-e)E + e( -E^3/3! + E^5/5! - ...)

   It's a little bit expensive to do this,  and you only need do it
quite rarely.  (I only encountered the problem because I had orbits
that were supposed to be 'pure parabolic',  but due to roundoff,
they had e = 1+/- epsilon,  with epsilon _very_ small.)  So 'near_parabolic'
is only called if we've gone seven iterations without converging. */

static double near_parabolic( const double ecc_anom, const double e)
{
   const double anom2 = (e > 1. ? ecc_anom * ecc_anom : -ecc_anom * ecc_anom);
   double term = e * anom2 * ecc_anom / 6.;
   double rval = (1. - e) * ecc_anom - term;
   int n = 4;

   while( fabs( term) > 1e-15)
      {
      term *= anom2 / (double)(n * (n + 1));
      rval -= term;
      n += 2;
      }
   return( rval);
}

/* For a full description of this function,  see KEPLER.HTM on the Guide
Web site,  http://www.projectpluto.com.  There was a long thread about
solutions to Kepler's equation on sci.astro.amateur,  and I decided to
go into excruciating detail as to how it's done below. */

#define MAX_ITERATIONS 7

static double kepler( const double ecc, double mean_anom)
{
   double curr, err, thresh, offset = 0., delta_curr = 1.;
   int is_negative = 0, n_iter = 0;

   if( !mean_anom)
      return( 0.);

   if( ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
      {
      curr = atan2( sin( mean_anom), cos( mean_anom) - ecc);
            /* two correction steps,  and we're done */
      for( n_iter = 2; n_iter; n_iter--)
         {
         err = curr - ecc * sin( curr) - mean_anom;
         curr -= err / (1. - ecc * cos( curr));
         }
      return( curr);
      }

   if( ecc < 1.)
      if( mean_anom < -PI || mean_anom > PI)
         {
         double tmod = fmod( mean_anom, PI * 2.);

         if( tmod > PI)             /* bring mean anom within -pi to +pi */
            tmod -= 2. * PI;
         else if( tmod < -PI)
            tmod += 2. * PI;
         offset = mean_anom - tmod;
         mean_anom = tmod;
         }

   if( mean_anom < 0.)
      {
      mean_anom = -mean_anom;
      is_negative = 1;
      }

   curr = mean_anom;
   thresh = THRESH * fabs( 1. - ecc);
   if( thresh < 1e-15)
      thresh = 1e-15;
   if( ecc > .8 && mean_anom < PI / 3. || ecc > 1.)    /* up to 60 degrees */
      {
      double trial = mean_anom / fabs( 1. - ecc);

      if( trial * trial > 6. * fabs(1. - ecc))   /* cubic term is dominant */
         {
         if( mean_anom < PI)
            trial = CUBE_ROOT( 6. * mean_anom);
         else        /* hyperbolic w/ 5th & higher-order terms predominant */
            trial = asinh( mean_anom / ecc);
         }
      curr = trial;
      }
   if( ecc > 1. && mean_anom > 4.)    /* hyperbolic, large-mean-anomaly case */
      curr = log( mean_anom);

   if( ecc < 1.)
      while( fabs( delta_curr) > thresh)
         {
         if( n_iter++ > MAX_ITERATIONS)
            err = near_parabolic( curr, ecc) - mean_anom;
         else
            err = curr - ecc * sin( curr) - mean_anom;
         delta_curr = -err / (1. - ecc * cos( curr));
         curr += delta_curr;
         printf( "Iter %d: %.16lf %.16lf %.16lf\n",
                           n_iter, curr, err, delta_curr);
         }
   else
      while( fabs( delta_curr) > thresh)
         {
         if( n_iter++ > MAX_ITERATIONS)
            err = -near_parabolic( curr, ecc) - mean_anom;
         else
            err = ecc * sinh( curr) - curr - mean_anom;
         delta_curr = -err / (ecc * cosh( curr) - 1.);
         curr += delta_curr;
         printf( "Iter %d: %.16lf %.16lf %.16lf\n",
                           n_iter, curr, err, delta_curr);
         }
   return( is_negative ? offset - curr : offset + curr);
}

void main( int argc, char **argv)
{
   double ecc = atof( argv[1]);
   double mean_anom = atof( argv[2]);

   printf( "E=%lf\n", kepler( ecc, mean_anom));
}
