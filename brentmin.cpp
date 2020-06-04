#include <assert.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include "brentmin.h"

/* A modified version of Brent's minimization algorithm.  You
first call brent_min_init( ) with three points bracketing the
minimum (bracketing is _required_).  Call brent_min_next( )
to determine the algorithm's suggested next point to look at;
evaluate the function there;  then call brent_min_add( ) to
add that result back in.  Repeat until the bracket is narrowed
to your satisfaction.  */

static void bubble_down( brent_min_t *b, int count)
{
   while( count--)
      if( b->y[count + 1] < b->y[count])
         {
         double tval = b->y[count + 1];

         b->y[count + 1] = b->y[count];
         b->y[count] = tval;
         tval = b->x[count + 1];
         b->x[count + 1] = b->x[count];
         b->x[count] = tval;
         }
}

#define PHI 0.6180339887498948482045868343656381177203091798057628621354486227052605

#define STEP_TYPE_INITIALIZED -1
#define STEP_TYPE_DONE        0
#define STEP_TYPE_GOLDEN      1
#define STEP_TYPE_CUBIC       2
#define STEP_TYPE_QUADRATIC   4

void brent_min_init( brent_min_t *b, const double x1, const double y1,
                                     const double x2, const double y2,
                                     const double x3, const double y3)
{
   b->x[0] = x1;
   b->y[0] = y1;
   b->x[1] = x2;
   b->y[1] = y2;
   b->x[2] = x3;
   b->y[2] = y3;
   b->xmin = b->xmax = x1;
   if( b->xmin > x2)
      b->xmin = x2;
   if( b->xmin > x3)
      b->xmin = x3;
   if( b->xmax < x2)
      b->xmax = x2;
   if( b->xmax < x3)
      b->xmax = x3;
   bubble_down( b, 2);     /* full bubble sort of three items */
   bubble_down( b, 2);
   assert( b->x[0] > b->xmin);
   assert( b->x[0] < b->xmax);
   b->gold_ratio = PHI;
   b->n_iterations = 0;
   b->prev_range = b->prev_range2 = 0.;
   b->tolerance = b->ytolerance = 0.;
   b->step_type = STEP_TYPE_INITIALIZED;
}

/* Fits a parabola y = ax^2 + bx + c to x[0, 1, 2] and y[0, 1, 2],  and
returns a,  which is always computed.  If b is non-NULL,  it's computed.
And if c is non-NULL,  it's computed also.  For some purposes -- specifically
those of determining the minimum point in a parabola -- you don't need c.
But I think this function may prove more generally useful (which is why
it's not static to this file). */

double fit_parabola( const double *x, const double *y, double *b, double *c)
{
   const double x21 = x[2] - x[1], x01 = x[0] - x[1], x20 = x[2] - x[0];
   const double y21 = y[2] - y[1], y01 = y[0] - y[1];
   double a;

   assert( x21 && x20 && x01);
   a = (y21 / x21 - y01 / x01) / x20;
   if( b)
      {
      *b = y01 / x01 - a * (x[1] + x[0]);
      if( c)
         *c = y[0] - x[0] * (*b + a * x[0]);
      }
   return( a);
}

/* This function takes four points,  x and y[0..3],  and finds the value
of x for which the cubic polynomial passing through all four points is a
minimum.  Similar to the use made above of fit_parabola( ),  except all we
care about here is the minimum x.

   To do this,  we shift the origin to x[0], y[0],  so that the cubic
passing through the origin and the remaining three (shifted) points is

      3     2
y = ax  + bx  + cx
 i    i     i     i

   (i=1, 2, 3).  So we have three equations and three unknowns.  A bit
of algebra takes advantage of the fact that this isn't just a general case
of three linear equations and three unknowns;  the coefficients are
squares and cubes of one another.  Eliminating c gets us

k12 = (y1 * x2 - y2 * x1) / (x1 * x2 * (x1 - x2)) = a * (x1 + x2) + b
k13 = (y1 * x3 - y3 * x1) / (x1 * x3 * (x1 - x3)) = a * (x1 + x3) + b

(k12 - k13) = a * (x2 - x3)

   The slope of y (remember,  we're after the minimum of the cubic here)
is dy/xy = 3ax^2 + 2bx + c,  and we're looking for zeroes of that.  Triple
a and double b,  and the quadratic we need to solve is plain ol ax^2+bx+c.

   There will actually be zero or two minima.  In the first case,  the
cubic is monotonic and we fail safely.  In the second,  we take the
minimum with smallest absolute value (i.e.,  closest to x0.)  */

static double cubic_min( const double *x, const double *y, int *err)
{
   const double x1 = x[1] - x[0], x2 = x[2] - x[0], x3 = x[3] - x[0];
   const double y1 = y[1] - y[0], y2 = y[2] - y[0], y3 = y[3] - y[0];
   const double k12_denom = x1 * x2 * (x1 - x2);
   const double k13_denom = x1 * x3 * (x1 - x3);
   double k12, k13, a, b, c, discr, rval;

   *err = 0;
   assert( k12_denom && k13_denom && x2 != x3);
   k12 = (y1 * x2 - y2 * x1) / k12_denom;
   k13 = (y1 * x3 - y3 * x1) / k13_denom;
   a = (k12 - k13) / (x2 - x3);
   b = k12 - a * (x1 + x2);
   c = y1 / x1 - x1 * (a * x1 + b);
   a *= 3.;
   b += b;
   discr = b * b - 4. * a * c;
   if( discr < 0.)      /* cubic is monotonically increasing or decreasing */
      {
      *err = 1;
      return( 0);
      }
   discr = sqrt( discr);
   if( b < 0.)                   /* quadratic formula rearranged */
      rval = -b + discr;         /* slightly to reduce loss of precision */
   else
      rval = -b - discr;
   return( x[0] + 2. * c / rval);
}

static int is_done( const brent_min_t *b)
{
   if( b->n_iterations > 3 && b->y[3] - b->y[0] <= b->ytolerance)
      return( 1);
   return( b->xmax - b->x[0] < b->tolerance && b->x[0] - b->xmin < b->tolerance);
}

double brent_min_next( brent_min_t *b)
{
   double rval;
   const double right = b->xmax - b->x[0], left = b->x[0] - b->xmin;
   const double range = b->xmax - b->xmin;

   assert( right);
   assert( left);
   if( is_done( b))
      {
      b->step_type = STEP_TYPE_DONE;
      return( b->x[0]);
      }
   b->step_type = STEP_TYPE_GOLDEN;
   if( b->n_iterations)
      {
      int err;

      rval = cubic_min( b->x, b->y, &err);
      if( !err)
         b->step_type = STEP_TYPE_CUBIC;
      }
   else     /* with only three points,  we'll try a quadratic step */
      {
      double quad, linear;

      quad = fit_parabola( b->x, b->y, &linear, NULL);
      rval = -linear * .5 / quad;
      b->step_type = STEP_TYPE_QUADRATIC;
      }
   b->prev_range2 = b->prev_range;
   b->prev_range = range;
         /* Unless it's a golden section step,  make sure we're at least
          b->tolerance away from x[0],  xmin,  and ymin.  (With a GS step,
          we'll be safe anyway.)  */
   if( b->step_type != STEP_TYPE_GOLDEN)
      {
      const double tol = b->tolerance * .9;

      if( rval > b->x[0] - tol && rval < b->x[0] + tol)
         rval = b->x[0] + (right > left ? tol : -tol);
      else if( rval < b->x[0])
         {
         if( rval < b->xmin)
            b->step_type = STEP_TYPE_GOLDEN;
         else if( rval < b->xmin + tol)
            rval = b->xmin + tol;
         }
      else           /* mirror image of preceding section */
         {
         if( rval > b->xmax)
            b->step_type = STEP_TYPE_GOLDEN;
         else if( rval > b->xmax - tol)
            rval = b->xmax - tol;
         }
      if( b->n_iterations > 30 && (b->n_iterations & 1))
         b->step_type = STEP_TYPE_GOLDEN;
      }
   if( b->step_type == STEP_TYPE_GOLDEN)
      rval = b->x[0] + (right > left ? right : -left) * (1. - b->gold_ratio);
   assert( rval >= b->xmin);
   assert( rval <= b->xmax);
   assert( rval != b->x[0]);
   b->next_x = rval;
   return( rval);
}

int brent_min_add( brent_min_t *b, const double next_y)
{
   int idx = (b->n_iterations ? 4 : 3);

   if( next_y <= b->y[0])      /* we have a new minimum */
      {
      if( b->next_x < b->x[0])
         b->xmax = b->x[0];
      else
         b->xmin = b->x[0];
      if( b->step_type == STEP_TYPE_GOLDEN)
         b->gold_ratio = (1. + b->gold_ratio) / 2.;
      }
   else                       /* no new minimum,  but the brackets can */
      {                       /* still be moved in */
      if( b->next_x < b->x[0])
         b->xmin = b->next_x;
      else
         b->xmax = b->next_x;
      if( b->step_type == STEP_TYPE_GOLDEN)
         b->gold_ratio = (PHI + b->gold_ratio) / 2.;
      }
   b->n_iterations++;
   while( idx && next_y <= b->y[idx - 1])
      {
      b->x[idx] = b->x[idx - 1];
      b->y[idx] = b->y[idx - 1];
      idx--;
      }
   b->x[idx] = b->next_x;
   b->y[idx] =    next_y;
   assert( b->y[0] <= b->y[1]);
   assert( b->y[1] <= b->y[2]);
   assert( b->y[2] <= b->y[3]);
   assert( b->x[0] > b->xmin);
   assert( b->x[0] < b->xmax);
   return( is_done( b));
}
