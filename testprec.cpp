/* testprec.cpp: functions for testing Earth precession code
from 'precess.cpp' and 'precess2.cpp'

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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"         /* for obliquity( ) prototype */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

/* Following are the "old" precession formulae,  in which the main */
/* precession is done in equatorial coordinates and precession in  */
/* ecliptic coordinates is derived from that.  In 'precess.cpp',   */
/* you'll see it's the other way around.  In hindsight,  that method */
/* fits the actual physics better,  and is the one I actually use.   */
/* I've preserved the "old" formulae here for testing/comparison.    */

/* setup_precession fills a 3x3 orthonormal matrix for precessing positions FROM */
/* year t1 TO year t2,  where t1 and t2 are Julian YEARS.  */

static int DLL_FUNC setup_precession_old( double DLLPTR *matrix,
                              double t1, double t2)
{
   double zeta, z, theta, czeta, cz, ctheta, szeta, sz, stheta;
   double ka, kb;
   static double t1_old = -PI, t2_old;
   static double curr_matrix[9];
   int going_backward = 0;

   if( fabs( t1 - t2) < 1.e-5)   /* dates sensibly equal;  spare the tedium */
      {                          /* of doing pointless math */
      set_identity_matrix( matrix);
      return( 0);
      }

         /* Ideally,  precessing from t1 to t2 back to t1 should get your  */
         /* original point.  To ensure that this happens,  we handle only  */
         /* the case t2 > t1;  otherwise,  we swap the times and invert    */
         /* the resulting matrix.                                          */
         /*   The reason is that the following precession formula uses     */
         /* cubic polynomials to approximate zeta,  theta,  and z.  If     */
         /* you feed it (t2,  t1),  it does _not_ create a matrix that is  */
         /* the exact inverse of (t1,  t2);  there is some accumulated     */
         /* error.  Doing it this way avoids having that show.  Also,      */
         /* there is a performance advantage:  if you _do_ call (t1, t2),  */
         /* then (t2, t1),  it's faster to invert the previous result      */
         /* than it would be to do all the math.                           */
   if( t1 < t2)
      {
      double temp = t1;

      t1 = t2;
      t2 = temp;
      going_backward = 1;
      }
         /* It's pretty common to precess a few zillion data points.  So   */
         /* it helps to cache the most recently computed precession matrix */
         /* so that repeated calls don't result in repeated computation.   */
   if( t1 == t1_old && t2 == t2_old)
      {
      FMEMCPY( matrix, curr_matrix, 9 * sizeof( double));
      if( going_backward)
         invert_orthonormal_matrix( matrix);
      return( 0);
      }
   t1_old = t1;
   t2_old = t2;

   t2 = (t2 - t1) / 100.;
   t1 = (t1 - 2000.) / 100.;
   ka = 2306.2181 + 1.39656 * t1 - .000139 * t1 * t1;
   kb = 2004.3109 - 0.85330 * t1 - .000217 * t1 * t1;
   zeta  = t2 * (ka + t2 * ( .30188 - .000345 * t1 + .017998 * t2));
   z     = t2 * (ka + t2 * (1.09468 + .000066 * t1 + .018203 * t2));
   theta = t2 * (kb + t2 * (-.42665 - .000217 * t1 - .041833 * t2));
   theta *= (PI / 180.) / 3600.;
   z     *= (PI / 180.) / 3600.;
   zeta  *= (PI / 180.) / 3600.;
   czeta = cos( zeta);
   szeta = sin( zeta);
   cz = cos( z);
   sz = sin( z);
   ctheta = cos( theta);
   stheta = sin( theta);

   *matrix++ = czeta * ctheta * cz - szeta * sz;
   *matrix++ = -szeta * ctheta * cz - czeta * sz;
   *matrix++ = -stheta * cz;

   *matrix++ = czeta * ctheta * sz + szeta * cz;
   *matrix++ = -szeta * ctheta * sz + czeta * cz;
   *matrix++ = -stheta * sz;

   *matrix++ = czeta * stheta;
   *matrix++ = -szeta * stheta;
   *matrix++ = ctheta;

   matrix -= 9;
   FMEMCPY( curr_matrix, matrix, 9 * sizeof( double));
   if( going_backward)
      invert_orthonormal_matrix( matrix);
   return( 0);
}

/* setup_ecliptic_precession fills a 3x3 orthonormal matrix for precessing */
/* positions _in ecliptic coordinates_ FROM year t1 TO year t2,  where t1  */
/* and t2 are Julian YEARS... much as setup_precession( ) does for RA/dec  */

/* 30 May 2002:  change 'obliquity#' to '-obliquity#' to fix a bug reported */
/* by Jordi Mas,  probably in place since the code was written.             */

static int DLL_FUNC setup_ecliptic_precession_old( double DLLPTR *matrix,
               const double t1, const double t2)
{
   const double obliquity1 = mean_obliquity( (t1 - 2000.) / 100.);
   const double obliquity2 = mean_obliquity( (t2 - 2000.) / 100.);

   setup_precession_old( matrix, t1, t2);
   pre_spin_matrix( matrix + 1, matrix + 2, -obliquity1);
   spin_matrix( matrix + 3, matrix + 6, -obliquity2);
   return( 0);
}

/* Some code used when I was switching over precession methods as
described above,  allowing me to compare results from the old code
(which I knew gave decent results near J2000) to the new.  It also
allowed me to try out one of Meeus' test cases.   */


static void show_matrix( const double *matrix)
{
   int i;

   for( i = 0; i < 3; i++, matrix += 3)
      printf( "%13.9f%13.9f%13.9f\n",
               matrix[0], matrix[1], matrix[2]);
}

int main( const int argc, const char **argv)
{
   double matrix_old[9], matrix_new[9];
   double vect[3], vect2[3];
   double max_diff, diff;
   int i;
   const double year_from = (argc > 1 ? atof( argv[1]) : 1950.);
   const double year_to = (argc > 2 ? atof( argv[2]) : 2000.);
   const double lon0 = 149.48194 * PI / 180.;
   const double lat0 = 1.76549 * PI / 180.;

   polar3_to_cartesian( vect, lon0, lat0);
   setup_ecliptic_precession_old( matrix_old, year_from, year_to);
   show_matrix( matrix_old);
   precess_vector( matrix_old, vect, vect2);
   printf( "%f %f\n", atan2( vect2[1], vect2[0]) * 180. / PI,
            asin( vect2[2]) * 180. / PI);

   printf( "\n  New method:\n");
   setup_ecliptic_precession( matrix_new, year_from, year_to);
   show_matrix( matrix_new);
   precess_vector( matrix_new, vect, vect2);
   printf( "%f %f\n", atan2( vect2[1], vect2[0]) * 180. / PI,
            asin( vect2[2]) * 180. / PI);

   for( i = 0, max_diff = 0.; i < 9; i++)
      if( max_diff < (diff = fabs( matrix_old[i] - matrix_new[i])))
         max_diff = diff;
   printf( "Maximum difference in matrices: %f\n", max_diff);

   printf( "\nEquatorial (old):\n");
   setup_precession_old( matrix_old, year_from, year_to);
   show_matrix( matrix_old);
   printf( "\nEquatorial (new):\n");
   setup_precession( matrix_new, year_from, year_to);
   show_matrix( matrix_new);

   for( i = 0, max_diff = 0.; i < 9; i++)
      if( max_diff < (diff = fabs( matrix_old[i] - matrix_new[i])))
         max_diff = diff;
   printf( "Maximum difference in matrices: %f\n", max_diff);
   return( 0);
}
