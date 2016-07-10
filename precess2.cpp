/* precess2.cpp: (deprecated version of) functions for computing
Earth precession;  see precess.cpp for current version,  and
'changes.txt' for info on why this is deprecated

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
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"         /* for obliquity( ) prototype */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923


/* setup_precession fills a 3x3 orthonormal matrix for precessing positions FROM */
/* year t1 TO year t2,  where t1 and t2 are Julian YEARS.  */

int DLL_FUNC setup_precession( double DLLPTR *matrix, double t1,
                               double t2)
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

static const double sin_obliq_2000 = 0.397777155931913701597179975942380896684;
static const double cos_obliq_2000 = 0.917482062069181825744000384639406458043;

void DLL_FUNC equatorial_to_ecliptic( double *vect)
{
   double temp;

   temp    = vect[2] * cos_obliq_2000 - vect[1] * sin_obliq_2000;
   vect[1] = vect[1] * cos_obliq_2000 + vect[2] * sin_obliq_2000;
   vect[2] = temp;
}

void DLL_FUNC ecliptic_to_equatorial( double *vect)
{
   double temp;

   temp    = vect[2] * cos_obliq_2000 + vect[1] * sin_obliq_2000;
   vect[1] = vect[1] * cos_obliq_2000 - vect[2] * sin_obliq_2000;
   vect[2] = temp;
}

int DLL_FUNC precess_vector( const double DLLPTR *matrix,
                                      const double DLLPTR *v1,
                                      double DLLPTR *v2)
{
   int i = 3;

   while( i--)
      {
      *v2++ = matrix[0] * v1[0] + matrix[1] * v1[1] + matrix[2] * v1[2];
      matrix += 3;
      }
   return( 0);
}

int DLL_FUNC deprecess_vector( const double DLLPTR *matrix,
                                      const double DLLPTR *v1,
                                      double DLLPTR *v2)
{
   int i = 3;

   while( i--)
      {
      *v2++ = matrix[0] * v1[0] + matrix[3] * v1[1] + matrix[6] * v1[2];
      matrix++;
      }
   return( 0);
}

int DLL_FUNC precess_ra_dec( const double DLLPTR *matrix,
                        double DLLPTR *p_out,
                        const double DLLPTR *p_in, int backward)
{
   double v1[3], v2[3];
   const double old_ra = p_in[0];

   v1[0] = cos( p_in[0]) * cos( p_in[1]);
   v1[1] = sin( p_in[0]) * cos( p_in[1]);
   v1[2] =                 sin( p_in[1]);
   if( backward)
      deprecess_vector( matrix, v1, v2);
   else
      precess_vector( matrix, v1, v2);
   if( v2[1] || v2[0])
      p_out[0] = atan2( v2[1], v2[0]);
   else
      p_out[0] = 0.;
   p_out[1] = asine( v2[2]);
   while( p_out[0] - old_ra > PI)
      p_out[0] -= PI * 2.;
   while( p_out[0] - old_ra <-PI)
      p_out[0] += PI * 2.;
   return( 0);
}

/* setup_ecliptic_precession fills a 3x3 orthonormal matrix for precessing */
/* positions _in ecliptic coordinates_ FROM year t1 TO year t2,  where t1  */
/* and t2 are Julian YEARS... much as setup_precession( ) does for RA/dec  */

/* 30 May 2002:  change 'obliquity#' to '-obliquity#' to fix a bug reported */
/* by Jordi Mas,  probably in place since the code was written.             */

int DLL_FUNC setup_ecliptic_precession( double DLLPTR *matrix, const double t1,
                              const double t2)
{
   const double obliquity1 = mean_obliquity( (t1 - 2000.) / 100.);
   const double obliquity2 = mean_obliquity( (t2 - 2000.) / 100.);

   setup_precession( matrix, t1, t2);
   pre_spin_matrix( matrix + 1, matrix + 2, -obliquity1);
   spin_matrix( matrix + 3, matrix + 6, -obliquity2);
   return( 0);
}

#ifdef TEST_MAIN
#include <stdio.h>
#include <stdlib.h>

int main( const int argc, const char **argv)
{
   double t1, t2, matrix[9];
   double p[2];
   int i;

   t1 = atof( argv[1]);
   t2 = atof( argv[2]);
   if( argc > 3)
      {
      p[0] = atof( argv[3]) * PI / 180.;
      p[1] = atof( argv[4]) * PI / 180.;
      }
   if( argc < 6)
      setup_precession( matrix, t1, t2);
   else
      setup_ecliptic_precession( matrix, t1, t2);
   for( i = 0; i < 9; i++)
      printf( "%15.11lf%s", matrix[i], (i % 3 == 2) ? "\n" : " ");
   if( argc > 3)
      {
      precess_ra_dec( matrix, p, p, 0);
      printf( "%lf %lf\n", p[0] * 180. / PI, p[1] * 180. / PI);
      precess_ra_dec( matrix, p, p, 1);
      printf( "%lf %lf\n", p[0] * 180. / PI, p[1] * 180. / PI);
      }
}
#endif
