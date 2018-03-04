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

#include <math.h>       /* for floor() prototype */
#include <stddef.h>     /* for NULL definition   */
#include <stdio.h>

double cubic_spline_interpolate_within_table(      /* spline.cpp */
         const double *table, const int n_entries, double x, int *err_code);
double lagrange_interpolate_within_table( const double *table,
         const int n_entries, double x, const int n_pts,
         double *deriv);                 /* spline.c */

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
xn + 1 to table,  so that values table[-1...2] will be used instead;
and we subtract xn + 1 from x,  so that 0. < x < 1. if we aren't
extrapolating.

   So we're looking for a cubic spline f(x) with a value of table[0]
for x=0;  a value of table[1] for x=1;  and a continuous first derivative
across grid points.  To accomplish that,  we say that the first derivative
at a grid point x=n is equal to (table[n+1] - table[n-1]) / 2.

   So now we're creating a cubic function f(x) where

   f(0) = table[0]
   f(1) = table[1]
   f'(0) = (table[1] - table[-1]) / 2
   f'(1) = (table[2] - table[0]) / 2

   With this,  we'll have a cubic that matches the grid points
exactly,  and have a continuous first derivative.  (There will
usually be discontinuities in second and higher derivatives.)

   f(x) = ax^3 + bx^2 + cx + table[0]                 (1)
   f(1) - f(0) = a + b + c                            (2)
   f'(0) = (table[1] - table[-1]) / 2 = c             (3)
   f'(1) = (table[2] - table[0]) / 2 = 3a + 2b + c    (4)

   Essentially,  we get the constant term of the cubic at no
charge;  it's table[0].  We get the linear term at not much
charge;  it's the first derivative at x=0.  To get the
cubic term a,  add equations (3) and (4) and subtract
twice equation (2) :

   a = (table[2] - table[0]) / 2 + c - 2. * (table[1] - table[0])   (5)

   Now that we have a and c,

   b = table[1] - table[0] - a - c                   (6)

   We now have all coefficients for the cubic (1) and can evaluate it.
*/

double cubic_spline_interpolate_within_table(
         const double *table, const int n_entries, double x, int *err_code)
{
   int idx = (int)floor( x);

   if( idx < 1)      /* extrapolate from front of table */
      {
      if( err_code)
         *err_code = (idx < 0 ? -2 : -1);
      idx = 1;
      }
   else if( idx > n_entries - 3)
      {             /* extrapolate beyond end of table */
      if( err_code)
         *err_code = (idx > n_entries - 2 ? -2 : -1);
      idx = n_entries - 3;
      }
   else             /* no extrapolation involved */
      if( err_code)
         *err_code = 0;
   table += idx;
   x -= (double)idx;

   const double c = (table[1] - table[-1]) * .5;
   const double y1 = table[1] - table[0];
   const double a = (table[2] - table[0]) * .5 - 2. * y1 + c;
   const double b = y1 - a - c;
   return( table[0] + x * (c + x * (b + x * a)));
}

/* One can pass a NULL 'deriv' if one doesn't actually need the first
derivative of the interpolated polynomial.  For a non-NULL value,  the
usual method for computing that derivative breaks down if you're exactly
on a grid point (divisions by zero occur).  In such cases,  numerical
differentiation is used instead.       */

double lagrange_interpolate_within_table( const double *table,
         const int n_entries, double x, const int n_pts,
         double *deriv)                  /* spline.c */
{
   int idx = (int)floor( x - (double)n_pts / 2.) + 1;
   double t = 1., c = 1., rval = 0., sum_recips = 0., y0;
   int i = (int)( x + .5);

   if( x == (double)i && i >= 0 && i < n_entries)
      {              /* we're exactly on a grid point */
      if( deriv)
         {
         const double epsilon = 1e-5;      /* numerically compute deriv */
         double y1, y2;

         y1 = lagrange_interpolate_within_table( table, n_entries, x - epsilon,
                              n_pts, NULL);
         y2 = lagrange_interpolate_within_table( table, n_entries, x + epsilon,
                              n_pts, NULL);
         *deriv = (y2 - y1) / (2. * epsilon);
         }
      return( table[i]);
      }
   if( idx < 0)      /* extrapolate from front of table */
      idx = 0;
   else if( idx > n_entries - n_pts)
      idx = n_entries - n_pts;   /* extrapolate beyond end of table */
   table += idx;
   x -= (double)idx;
   y0 = table[n_pts / 2];

   for( i = 0; i < n_pts; i++)
      {
      const double dx = x - (double)i;

      if( deriv)
         sum_recips += 1. / dx;
      c *= dx;
      if( i)
         t *= -(double)i;
      }

   if( deriv)
      *deriv = 0.;
   for( i = 0; i < n_pts; i++)
      {
      const double dx = x - (double)i;

      if( i)
         t *= (double)i / (double)( i - n_pts);
      rval += (table[i] - y0) / (t * dx);
      if( deriv)
         *deriv += (sum_recips - 1. / dx) * table[i] / (t * dx);
      }
   rval *= c;
   rval += y0;
   if( deriv)
      *deriv *= c;
   return( rval);
}

#ifdef TEST_CODE

#include <math.h>
#include <stdio.h>

static void show_explanation( void)
{
   printf(
       "See 'spline.cpp' comments for mathematical details.  This test code\n"
       "creates a table of sixteen points spaced at 0.1 radians of the sine\n"
       "function,  and shows the result of cubic splines and Lagrange\n"
       "interpolation within that table.  The interpolation is done at ten\n"
       "degree intervals,  i.e.,  it doesn't line up with the table points.\n"
       "The table covers 0 to 1.5 radians (about 85 degrees).  The interpolation\n"
       "starts at -30 degrees and runs to 260 degrees to show the effects of\n"
       "extrapolation.  As you'd expect,  accuracy is good to excellent within\n"
       "the table and deteriorates as you extrapolate.\n\n");
}


static void show_header( const int order1, const int order2, const int order_step)
{
   int j;

   printf( "Ang    sin(ang)         cubic   ");
   for( j = order1; j <= order2; j += order_step)
      printf( "         %2d-order", j);
   printf( "\n");
}

#define NPTS 16

int main( const int argc, const char **argv)
{
   double table[NPTS];
   int i, j;
   int order1 = 4, order2 = 12, order_step = 2;
   bool show_differences = false, test_derivs = false;
   const double PI =
           3.1415926535897932384626433832795028841971693993751058209749445923;
   const double scale = PI / 30.;

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'd':
               printf( "Showing differences\n");
               show_differences = true;
               break;
            case '1':
               printf( "Testing first derivatives\n");
               test_derivs = true;
               break;
            case 'o':
               sscanf( argv[i] + 2, "%d,%d,%d", &order1, &order2, &order_step);
               break;
            }
   show_explanation( );
   for( i = 0; i < NPTS; i++)
      table[i] = sin( (double)i * scale);
   show_header( order1, order2, order_step);
   for( i = -30; i < 270; i += 10)
      {
      const double angle = (double)i * PI / 180.;
      const double x = angle / scale;
      const double subtract = (show_differences ? sin( angle) : 0.);

      printf( "%3d %16.13f %16.13f", i, sin( angle),
                cubic_spline_interpolate_within_table( table, NPTS, x, NULL)
                                 - subtract);
      for( j = order1; j <= order2; j += order_step)
         {
         double value_to_show, deriv;

         value_to_show = lagrange_interpolate_within_table( table, NPTS, x, j,
                        &deriv) - subtract;
         if( test_derivs)
            value_to_show = deriv / scale - (show_differences ? cos( angle) : 0.);
         printf( " %16.13f", value_to_show);
         }
      printf( "\n");
      }
   show_header( order1, order2, order_step);
}
#endif
