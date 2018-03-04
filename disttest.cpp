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
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"

/* Test routine for position angle/distance code.

   I'm not using the library rand( ) function,  because I want to be
confident that the results will be the same when compiled on different
systems. Following is the MMIX linear congruential pseudo-random number
generator, copied from Donald Knuth.  It's not a great PRNG,  but it's
fine for the current humble purpose.  */

uint64_t pseudo_random( const uint64_t prev_value)
{
   uint64_t a = UINT64_C( 6364136223846793005);
   uint64_t c = UINT64_C( 1442695040888963407);

   return( a * prev_value + c);
}

static double get_pseudorandom_double( void)
{
   static uint64_t rval = 1;
   const double two_to_the_53rd_power = 9007199254740992.0;

   rval = pseudo_random( rval);
   return (int64_t)(rval >> 11) * (1.0 / two_to_the_53rd_power);
}

/* Should provide RA/dec values uniformly distributed over the */
/* celestial sphere. */

static void get_random_ra_dec( double *ra, double *dec)
{
   double x, y, z, dist;

   do
      {
      x = get_pseudorandom_double( ) - .5;
      y = get_pseudorandom_double( ) - .5;
      z = get_pseudorandom_double( ) - .5;
      }
      while( (dist = x * x + y * y + z * z) > .125);
   *ra = atan2( y, x);
   *dec = asin( z / sqrt( dist));
}

static double calc_dist_with_acos( double *p1, double *p2)
{
   double d_ra = p1[0] - p2[0];
   double cos_dist = sin( p1[1]) * sin( p2[1]) + cos( p1[1]) * cos( p2[1]) * cos( d_ra);

   return( acos( cos_dist));
}

int main( const int argc, const char **argv)
{
   int i;

   for( i = 0; i < 10000; i++)
      {
      double p1[2], p2[2], dist1, dist2;
      const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923;

      get_random_ra_dec( p1, p1 + 1);
      get_random_ra_dec( p2, p2 + 1);
      if( !i && argc == 5)
         {
         p1[0] = atof( argv[1]) * pi / 180.;
         p1[1] = atof( argv[2]) * pi / 180.;
         p2[0] = atof( argv[3]) * pi / 180.;
         p2[1] = atof( argv[4]) * pi / 180.;
         }
      calc_dist_and_posn_ang( p1, p2, &dist1, NULL);
      dist2 = calc_dist_with_acos( p1, p2);
      printf( "%12.8lf %12.8lf %lg\n",
                  dist1 * 180. / pi,
                  dist2 * 180. / pi,
                  (dist1 - dist2) * 180. / pi);
      }
   return( 0);
}
