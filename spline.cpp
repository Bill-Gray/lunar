#include <math.h>       /* for floor() prototype */

double cubic_spline_interpolate_within_table(      /* spline.cpp */
         const double *table, const int n_entries, double x, int *err_code);

/* The following cubic_spline_interpolate_within_table( ) function
assumes you have a table of n_entries values in an array table[]
and wish to get an interpolated value at x.  The function is
continuous and so is its first derivative.

   First,  we'll want to establish which points in table[] will
be used.  If xn = (int)floor(x) - 1,  then table[xn]...table[xn+3]
will bracket x,  with two points to either side of x.  (If xn < 0,
we set xn=0;  if xn > n_entries - 3,  we set xn=n_entries - 3;
in either case,  we end up extrapolating rather than interpolating.
'err_code' is zero if we've two points on each side of x;  it's
-1 if we've only one point on one side and three on the other;
and it's -2 if we're completely outside the table.)

   With that accomplished,  we need to set up a cubic spline between
the points table[xn...xn+3].  To simplify life a little,  we add
xn to table,  so that values table[0...3] will be used instead;
and we subtract xn+1.5 from x,  so that -.5 < x < .5 if we aren't
extrapolating.

   The spline has to match the values in table[] at integer
values,  and the first derivative has to match from one spline to
the next.  To accomplish that,  we say that the first derivative
at a grid point x=n is equal to (table[n+1] - table[n-1]) / 2.

   So now we're creating a cubic spline of a function where

   y0 = f(-1.5) = table[0]
   y1 = f(-0.5) = table[1]
   y2 = f( 0.5) = table[2]
   y3 = f( 1.5) = table[3]

   ...where the cubic is required to match at .5 and -.5,  and the
slope has to be (y2-y0)/2 at x=-.5 and (y3-y1) at x=.5.  This lets
us have a continuous function over -.5 to +.5,  with the slopes matching
up across grid points.  Given the coefficients of that cubic,  we
can compute an interpolated value for -.5 < x < .5.

   So... cubic(x) = ax^3 + bx^2 + cx + d,  so
         cubic'(x) = 3ax^2 + 2bx + c

   cubic'(-.5) = 3a/4 - b + c = (y2 - y0) / 2       (1)
   cubic'( .5) = 3a/4 + b + c = (y3 - y1) / 2       (2)
         -> b = (y3 - y1 - y2 + y0) / 4             (3) = (2) - (1), halved
         and 3a/2 + 2c = (y2 - y0 + y3 - y1) / 2    (4) = (1) + (2)

   cubic( 0.5) =  a/8 + b/4 + c/2 + d = y2          (5)
   cubic(-0.5) = -a/8 + b/4 - c/2 + d = y1          (6)
         -> b/2 + 2d = y2 + y1                      (7) = (5) + (6)
         -> d = (y2 + y1) / 2. - b / 4              (8) = (7) rearranged
         and a/4 + c = y2 - y1                      (9) = (5) - (6)

   OK.  Double (9) and subtract it from (4) and you get

   a = -1.5y2 - y0/2 + y3/2 + 1.5y1 = 1.5(y1-y2) + (y3-y0)/2  (10)

   ...and feed a back into (9) and get, pretty trivially,

   c = y2 - y1 - a / 4                              (11)

   OK,  I admit,  it's all pretty simple math.  I just wrote it
out so I didn't lose track of what I'd done.  You'll notice that
this is all accomplished in about half as many lines as it takes
to describe it...  */

double cubic_spline_interpolate_within_table(
         const double *table, const int n_entries, double x, int *err_code)
{
   int idx = (int)floor( x) - 1;

   if( idx < 0)      /* extrapolate from front of table */
      {
      if( err_code)
         *err_code = (idx < -1 ? -2 : -1);
      idx = 0;
      }
   else if( idx > n_entries - 4)
      {             /* extrapolate beyond end of table */
      if( err_code)
         *err_code = (idx > n_entries - 3 ? -2 : -1);
      idx = n_entries - 2;
      }
   else             /* no extrapolation involved */
      if( err_code)
         *err_code = 0;
   table += idx;
   x -= (double)idx + 1.5;

   const double y1_minus_y2 = table[1] - table[2];
   const double a = 1.5 * y1_minus_y2 + .5 * (table[3] - table[0]);  /* (10) */
   const double c = -.25 * a - y1_minus_y2;                          /* (11) */
   const double y1_plus_y2 = table[1] + table[2];
   const double b = (table[3] + table[0] - y1_plus_y2) / 4.;         /* (3) */
   const double d = .5 * y1_plus_y2 - b * .25;                       /* (8) */
   return( d + x * (c + x * (b + x * a)));
}
