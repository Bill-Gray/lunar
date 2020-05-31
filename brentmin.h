#ifndef BRENTMIN_H_INCLUDED
#define BRENTMIN_H_INCLUDED

/* brentmin.h: header file for Brent minimization without derivatives

Copyright (C) 2018, Project Pluto

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

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

typedef struct
{
   double x[5], y[5];         /* lowest values thus far */
   double xmin, xmax;         /* bracketing values */
   double next_x;
   double gold_ratio;         /* variable;  see brentmin.cpp */
   double tolerance, ytolerance;
   double prev_range, prev_range2;
   int step_type, n_iterations;
} brent_min_t;

void brent_min_init( brent_min_t *b, const double x1, const double y1,
                                     const double x2, const double y2,
                                     const double x3, const double y3);
double brent_min_next( brent_min_t *b);
int brent_min_add( brent_min_t *b, const double next_y);
double fit_parabola( const double *x, const double *y, double *b, double *c);

/* This is a somewhat modified version of Brent's algorithm for finding the
minimum of a single-variable function.  brent_min_init( ) takes three points
with distinct x values,  such that the middle x value has the lowest y value;
i.e.,  it's a bracketed minimum.  (The points need not be passed in any
particular order).

   One then enters a loop wherein one calls brent_min_next( ) to get the
next suggested point to check;  you evaluate the function there and call
brent_min_add( ) to update the algorithm.  Repeat until xmax-xmin is below
a desired threshhold.

   Brent's method combines "golden search" (a not very fast method,  but
with guaranteed linear convergence) and parabolic interpolation : find
the quadratic fitting through the three lowest points seen thus far,  and
look at the point that is the minimum of that quadratic.  The latter
method is much faster,  _if_ the function really resembles a parabola
within the bracket.  Some bookkeeping is required to detect cases where
convergence is slow;  if that happens,  we switch back to the reliable
golden section method for a step or two.  Thus,  even in worst-case
scenarios,  the algorithm _will_ still converge linearly.

   I've made the following modifications :

   -- This code keeps track of the _four_ lowest points and fits a _cubic_
polynomial to them,  and determines a minimum from that.

   -- Convergence to the minimum tends to be extremely fast for sufficiently
"well-behaved" functions.  But even with those,  you can have problems because
you've reduced the (say) left-hand side of the bracket a thousandfold,  but the
right-hand side is as large as ever.  In such cases, we're probably very close
to the minimum,  and don't have to go far into the larger side to find a point
that is,  most likely,  past the maximum.  This "shrink" step can easily reduce
xmax-xmin a hundredfold,  and more than that in the final steps where the
function is effectively a parabola or cubic.

   At present,  if the large side of the bracket is more than ten times
the size of the small side of the bracket,  we take a look at a point that
is at xmin - (xsmall - xmin) * 3,  where xmin is our minimum point and
xsmall is the point at the small end of the bracket.  In other words,
we take the small end's size,  triple it,  and look that far into the
larger bracket.

   If we're lucky (as we usually are),  this drastically reduces the size
of the larger side of the bracket.  If not,  the smaller side is tripled
in size (and the larger size shrinks accordingly).

   This is the result of some empirical testing,  and should be examined
more carefully.   I'd like to have something on a more solid theoretical
ground.  One possibility is to compute both the cubic and parabolic minima;
the difference between them could be a gauge of the error to which we've
determined the minimum.  Or perhaps the magnitude of the leading coefficient
of the cubic polynomial gives us a clue as to how far the function departs
from being a parabola.
*/

#ifdef __cplusplus
}
#endif  /* #ifdef __cplusplus */
#endif  /* #ifndef BRENTMIN_H_INCLUDED */
