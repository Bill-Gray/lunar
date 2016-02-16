/* dist_pa2: distance/position angle computations
(NOT VERY GOOD;  SEE 'DIST_PA.CPP' FOR A BETTER VERSION!  This
one is mostly for reference use only.)

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

/*
   NOTE that this version is deprecated.  I strongly suggest using
'dist_pa.cpp' instead.  In this older version,  the same problems are
tackled using classical spherical trigonometric methods.  That turns
out to be excessively complicated and prone to lots of "corner cases"
wherein precision is lost.  The 'dist_pa.cpp' is a lot simpler and
more precise.

   calc_dist_and_posn_ang() takes two RA/dec points,  and computes their
angular separation and position angle.  If you don't want the PA,  pass
in NULL for the posn_ang argument.

   After a quick check to make sure p1 != p2,  the code uses a haversine
method to get the distance.  This handles small angles correctly.  The
"usual" formula for such a spherical triangle,

cos(a) = cos(b) * cos(c) + sin(b) * sin(c) * cos( alpha)

   has serious precision problems when cos(a) is close to +/-1 (a=0 or
180 degrees).  Instead,  the haversine formula mentioned in Meeus'
_Astronomical Algorithms_,  page 111,  is used:  (use a fixed-size
font to read the following equations:)

hav( dist) = hav( dec1 - dec2) + cos( dec1) * cos( dec2) * hav( delta_ra)

   A quick review of haversines:  it helps to know that if

                     2
x = hav( theta) = sin ( theta / 2) = (1 - cos( theta) / 2)

   then the "inverse haversine" is computed as

theta = inv_haversine( x) = 2 * arcsine( sqrt( x))

   Once it has the distance,  getting the PA is just about trivial,  and
is again done using haversines,  using the spherical triangle formula

   2                              cos( b - c) - cos( a)
sin ( alpha / 2) = hav( alpha) = -----------------------
                                    2 sin( b) sin( c)

   ...which,  after setting b = 90 - dec1, c = dist, a = 90 - dec2,
alpha = PA (from point 1 to point 2),  gives

            sin( dec1 + dist) - sin( dec2)
hav( PA) = --------------------------------
                 2 cos( dec1) sin( dist)

   In using this,  the code first makes sure that sin(b) * sin(c) isn't
zero (which could happen if the points are superimposed, making dist=0, or
if the starting point is at the pole,  dec1 = +/- 90... in either case,
the PA is undefined).  After doing the division,  it makes sure the result
is <= 1;  mathematically,  it never will be, but roundoff happens.

   The resulting PA lies between 0 and 180,  because the above formula
can't distinguish east from west.  So we do a check of sin( ra1 - ra2)
to fix that.

   (14 Nov 2005)  John Greaves pointed out that the position error
reported using this function is undefined in some cases where the
position angle is zero.  I now see that if PA=0,  then dec1 + dist = dec2
and in the above haversine formula for the position angle,  we're
subtracting identical quantities.  And if there's roundoff (there is
always roundoff),  you might get a negative haversine,  and a square
root error.  So we now check to make sure the value isn't less than zero,
as well as making sure it's less than or equal to one.

   Ideally,  we might make use of the fact that for nearly identical
angles a and b,  it's best to express the difference in their sines
as

sin(a) - sin(b) = sin( (a+b)/2) * cos((a-b)/2)
*/

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define HALF_PI (PI / 2.)

int DLL_FUNC calc_dist_and_posn_ang( const double *p1, const double *p2,
                                       double *dist, double *posn_ang)
{
   double d_ra = p1[0] - p2[0];
   int rval = -1;

   if( d_ra || p1[1] != p2[1])
      {
      double tval;

      tval = (1. - cos( p1[1] - p2[1])) / 2. + cos( p1[1]) * cos( p2[1])
            * (1. - cos( d_ra)) / 2.;
      if( tval >= 1.)        /* points are at 'antipodes'.  Roundoff  */
         {                   /* errors can make tval > 1.             */
         *dist = PI;
         if( posn_ang)
            *posn_ang = 0.;
         }
      else
         {
         *dist = 2. * asin( sqrt( tval));
         if( posn_ang)
            {
            tval = cos( p1[1]) * sin( *dist);
            if( !tval)           /* p1[1] is at a pole... undefined p.a. */
               *posn_ang = 0.;
            else
               {
               const double hav_pa = (sin( p1[1] + *dist) - sin( p2[1]))
                            / (2. * tval);

               if( hav_pa >= 1.)        /* again,  roundoff troubles */
                  *posn_ang = PI;
               else if( hav_pa <= 0.)  /* roundoff: see 14 Nov 05 note */
                  *posn_ang = 0.;
               else
                  *posn_ang = 2. * asin( sqrt( hav_pa));
               if( sin( d_ra) < 0.)         /* Fixed 23 aug 94 to run */
                  *posn_ang = PI + PI - *posn_ang;   /* 90=E,  not W */
               }
            }
         }
      rval = 0;
      }
   else           /* points are identical */
      {
      *dist = 0.;
      if( posn_ang)          /* posn angle is really undefined,  but let's */
         *posn_ang = 0.;     /* zero it to evade "Not A Number" messages */
      }
   return( rval);
}
