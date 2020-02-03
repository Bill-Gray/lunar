#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "brentmin.h"

/* Code for testing 'brentmin.cpp',  a modified version of Brent's
single-variable minimization routine.  Note that the starting
values specified on the command line _must_ bracket a minimum.  */

double func( const double x)
{
   return( -x * exp( -x));    /* minimum at 1 */
// return( sin( x));          /* minima at (2n+1.5) * pi */
}

int main( const int argc, const char **argv)
{
   double x[3];
   int i, is_done = 0;
   brent_min_t b;

   assert( argc >= 4);
   for( i = 0; i < 3; i++)
      x[i] = atof( argv[i + 1]);
   brent_min_init( &b, x[0], func( x[0]), x[1], func( x[1]),
                                          x[2], func( x[2]));
   if( argc == 5)
      b.tolerance = atof( argv[4]);
   while( !is_done && b.n_iterations < 30)
      {
      const double new_x = brent_min_next( &b);
      const double new_y = func( new_x);
      const char *types[5] = { "Done", "Golden", "cubic", "shrink", "quadratic"  };

      assert( new_x == b.next_x);
      assert( new_x != b.xmax);
      assert( new_x != b.xmin);
      assert( new_x < b.xmax);
      assert( new_x > b.xmin);
      is_done = brent_min_add( &b, new_y);
      printf( "%2d %f %f (%f to %f, range %.12f) %s; gold=%f\n",
               b.n_iterations, new_x, new_y, b.xmin, b.xmax,
               b.xmax - b.xmin, types[b.step_type], b.gold_ratio);
      printf( "%.15f %.15f %.15f\n", b.x[1] - b.x[0], b.x[2] - b.x[0], b.x[3] - b.x[0]);
      }
   if( b.tolerance)
      printf( is_done ? "Converged\n" : "FAILED TO CONVERGE TO TOLERANCE\n");
   return( 0);
}
