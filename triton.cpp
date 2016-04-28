/* triton.cpp: low-precision ephems for Triton

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
#include "lunar.h"
#include "afuncs.h"

#define J2000  2451545.
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

/* For the following,  see p 373 of the _Explanatory Supplement_ */
/* Note that 'rocks.cpp' also has code for computing the position of Triton.
The following is,  therefore,  essentially obsolete.   */

void DLL_FUNC calc_triton_loc( const double jd, double *vect)
{
   const double t_cent = (jd - J2000) / 36525.;
   const double n = (359.28 + 54.308 * t_cent) * (PI / 180.);
   const double t0 = 2433282.5;
   const double theta = (151.401 + .57806 * (jd - t0) / 365.25) * (PI / 180.);
            /* Semimajor axis is 488.49 arcseconds at one AU: */
   const double semimajor = 488.49 * (PI / 180.) / 3600.;
   const double longitude =
               (200.913 + 61.2588532 * (jd - t0)) * (PI / 180.);
   const double gamma = 158.996 * (PI / 180.);

         /* Calculate longitude and latitude on invariable plane: */
   const double lon_on_ip = theta + atan2( sin( longitude) * cos( gamma),
                                     cos( longitude));
   const double lat_on_ip = asin( sin( longitude) * sin( gamma));
         /* Vectors defining invariable plane,  expressed in B1950: */
   double x_axis[3], y_axis[3], z_axis[3];
         /* Vector defining Triton position in invariable plane space: */
   double triton[3];
         /* RA/dec of the pole: */
   double ra_dec_p[2];
   double matrix[9];
   double vect_1950[3];
   int i;

   polar3_to_cartesian( triton, lon_on_ip, lat_on_ip);

   ra_dec_p[0] = (298.72 * PI / 180.) + (2.58 * PI / 180.) * sin( n)
                     - (0.04 * PI / 180.) * sin( n + n);
   ra_dec_p[1] = (42.63 * PI / 180.) - (1.90 * PI / 180.) * cos( n)
                     + (0.01 * PI / 180.) * cos( n + n);
   setup_precession( matrix, 1950., 2000.);
   precess_ra_dec( matrix, ra_dec_p, ra_dec_p, 1);
   polar3_to_cartesian( x_axis, ra_dec_p[0] + PI / 2., 0.);
   polar3_to_cartesian( y_axis, ra_dec_p[0] + PI, PI / 2. - ra_dec_p[1]);
   polar3_to_cartesian( z_axis, ra_dec_p[0], ra_dec_p[1]);

   for( i = 0; i < 3; i++)
      vect_1950[i] = semimajor * (x_axis[i] * triton[0] +
                      y_axis[i] * triton[1] + z_axis[i] * triton[2]);
   precess_vector( matrix, vect_1950, vect);
}

