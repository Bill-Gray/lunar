/* dist_pa.cpp: functions for distance/position angle computations

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

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

/* Code to compute distance and position angle between two points on a
sphere;  or,  given one point and a distance/PA,  find the second point.
NOTE that this is for a _sphere_.  For an ellipsoidal solution,  see
'dist.cpp'.          */

/* (6 Apr 2008) Added a "reverse" function, reverse_dist_and_posn_ang(),
which takes a point, distance, and PA, and computes the destination point.

   (19 Feb 2006)  This function to compute the distance and position
angle between two points has been completely rewritten relative to
the previous one.  It uses a very different approach intended to
evade some loss of precision problems that arose with that previous
(haversine-based) approach.

   We'll assume two points (theta1, phi1) and (theta2, phi2).  On the
Earth,  theta=lon and phi=lat;  on the sky,  theta=RA and phi=dec.
We can first subtract theta1 from both quantities so that the points
become (0, phi1) and (theta2-theta1, phi2) = (dtheta, phi2);  this
rotation by theta1 around the z-axis won't affect the distance and
position angle between the two points.

   We can then use the following spherical trig identities:

z_component = cos(dist) = cos(dtheta) * cos(phi1) * cos(phi2)
                                      + sin(phi1) * sin(phi2);
x_component = sin(dist)sin(PA) = sin(dtheta) * cos(phi2);
y_component = sin(dist)cos(PA) = sin(phi2) * cos(phi1) -
                                 sin(phi1) * cos(phi2) * cos(dtheta);

   The "component" notation came about because I first came up with this
method by thinking in terms of vectors:  in a coordinate system where
the first RA/dec coordinate points in the z-axis,  and the x-coordinate
is in the equator at right angles to this,  and the y-axis completes
the triplet,  then the second RA/dec coordinate pair has the unit vector
(x_component, y_component, z_component) in that system.  But the above
also corresponds to spherical trigonometric identities.

   Things always seem to break down for small dtheta and phi1 close to
phi2,  so I made use of trig subtraction identities to rearrange the
y and z components to read

z_component = cos( phi2 - phi1) - (1 - cos(dtheta)) * cos(phi1) * cos(phi2)
y_component = sin( phi2 - phi1) + (1 - cos(dtheta)) * sin(phi1) * cos(phi2)

   ...and we can then get the position angle right away,  as the atan2
of the x and y components (suitably moved around to satisfy the usual
conventions that PA=0 is north and so on.)

   If the z_component is between +/-.7,  we can just get the distance
as arccos( z_component).  But for small distances or distances near the
antipodes (z_component near +/- 1),  roundoff will give us grief,  and
we have to use sin( dist) = sqrt( x_component^2 + y_component^2).  Then
we check to see which quadrant we're in (i.e.,  do we have a small distance
or one near the poles?),  and we're done.

   I don't see where a loss-of-precision problem would arise with this
code.  Using the 'test function',  I see that this code remains accurate:
the original calc_dist_and_posn_ang( ) function loses precision for very
small distances and those near to 180 degrees (i.e.,  the points are at
antipodes).  In such cases,  the haversine of the position angle is
poorly rounded off,  and you get an imprecise PA and sometimes distance.
This new code appears to handle both cases with aplomb.   */

int DLL_FUNC calc_dist_and_posn_ang( const double *p1, const double *p2,
                                       double *dist, double *posn_ang)
{
   int rval = 0;
   const double one_minus_cos_dtheta = 1. - cos( p1[0] - p2[0]);
   const double cos_dec2 = cos( p2[1]);
   const double delta_dec = p2[1] - p1[1];
   const double x_component = sin( p2[0] - p1[0]) * cos_dec2;
   const double y_component = sin( delta_dec) +
                  one_minus_cos_dtheta * sin( p1[1]) * cos_dec2;
   const double z_component = cos( delta_dec) -
                  one_minus_cos_dtheta * cos( p1[1]) * cos_dec2;

   if( z_component > -.7 && z_component < .7)
      *dist = acos( z_component);
   else
      {
      const double sin_dist = sqrt( x_component * x_component
                                  + y_component * y_component);

      *dist = asin( sin_dist);
      if( z_component < 0.)
          *dist = PI - *dist;
      }
   if( posn_ang)
      {
      if( x_component != 0. || y_component != 0.)
         {
         *posn_ang = atan2( y_component, x_component) - PI / 2.;
         if( *posn_ang < 0.)
            *posn_ang += PI + PI;
         }
      else                 /* undefined position angle */
         {
         *posn_ang = 0.;
         rval = -1;
         }
      }
   return( rval);
}

void DLL_FUNC reverse_dist_and_posn_ang( double *to, const double *from,
                                 const double dist, const double posn_ang)
{
   const double cos_lat = cos( from[1]), sin_lat = sin( from[1]);
   const double sin_dist = sin( dist), cos_dist = cos( dist);
   const double sin_pa = sin( posn_ang);
   const double sin_dist_cos_pa = sin_dist * cos( posn_ang);
   const double x = cos_lat * cos_dist - sin_lat * sin_dist_cos_pa;
   const double y = sin_dist * sin_pa;
   const double z = sin_lat * cos_dist + cos_lat * sin_dist_cos_pa;

   to[0] = from[0];
   if( x != 0. || y != 0.)
      to[0] -= atan2( y, x);
                       /* See above comments about evading loss of   */
                       /* precision issues when z is close to +/- 1: */
   if( z < .7 && z > -.7)
      to[1] = asin( z);
   else
      {
      to[1] = acos( sqrt( x * x + y * y));
      if( z < 0.)
         to[1] = -to[1];
      }
}

#ifdef TEST_MAIN

#include <stdlib.h>
#include <stdio.h>
#include "watdefs.h"
#include "afuncs.h"

void main( int argc, char **argv)
{
   double dist, posn_ang;
   double points[4];
   int i;

   for( i = 0; i < 4; i++)
      points[i] = atof( argv[i + 1]) * PI / 180.;

   printf( "(New) rval = %d: ", calc_dist_and_posn_ang( points, points + 2,
                                    &dist, &posn_ang));
   printf( "Dist %.10lf, PA %lf\n", dist * 180. / PI, posn_ang * 180. / PI);

   printf( "(Old) rval = %d: ", calc_dist_and_pa( points, points + 2,
                                    &dist, &posn_ang));
   printf( "Dist %.10lf, PA %lf\n", dist * 180. / PI, posn_ang * 180. / PI);
}
#endif
