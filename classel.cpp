/* classel.cpp: converts state vects to classical elements

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
#include <assert.h>
#ifndef __cplusplus
   #include <stdbool.h>
#endif
#include "watdefs.h"
#include "afuncs.h"
#include "comets.h"

/* MS only got around to adding 'isfinite',  asinh in VS2013 : */

#if defined( _MSC_VER) && (_MSC_VER < 1800)
#include <float.h>
#define isfinite _finite

static double asinh( const double x)
{
   return( log( x + sqrt( x * x + 1.)));
}
#endif

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define SQRT_2 1.4142135623730950488016887242096980785696718753769480731766797

/* 2011 Aug 11:  in dealing with exactly circular orbits,  and those
   nearly exactly circular,  I found several loss of precision problems
   in the angular elements (which become degenerate for e=0;  you can
   add an arbitrary amount to the mean anomaly,  as long as you subtract
   the same amount from the argument of periapsis.)  Also,  q was
   computed incorrectly for such cases when roundoff error resulted in
   the square root of a number that should be zero,  but rounded below
   zero,  was taken.  All this is fixed now.       */

/* 2009 Nov 24:  noticed a loss of precision problem in computing arg_per.
   This was done by computing the cosine of that value,  then taking the
   arc-cosine.  But if that value is close to +/-1,  precision is lost
   (you can actually end up with a domain error if the roundoff goes
   against you).  I added code so that,  if |cos_arg_per| > .7,  we
   compute the _sine_ of the argument of periapsis and use that instead.

   While doing this,  I also noticed that several variables could be made
   of type const.   */

/* calc_classical_elements( ) will take a given state vector r at a time t,
   for an object orbiting a mass gm;  and will compute the orbital elements
   and store them in the elem structure.  Normally,  ref=1.  You can set
   it to 0 if you don't care about the angular elements (inclination,
   longitude of ascending node,  argument of perihelion).         */

/* In determining the mean anomaly from the eccentricity and
eccentric anomaly,  use of the "normal" formulae

M = E - ecc * sin( E)         (elliptical case)
M = E - ecc * sinh( E)        (hyperbolic case)

   you run into nasty loss-of-precision problems with near-parabolic orbits:
E and ecc * sin(E) will be nearly equal quantities.  In such cases,  it's
better to use power series for the sin/sinh function, rearranging to get

M = E(1-ecc) - ecc( -E^3/3! + E^5/5! - E^7/7!...)  (elliptical)
M = E(1-ecc) - ecc(  E^3/3! + E^5/5! + E^7/7!...)  (hyperbolic)

   ...in which the infinite series is that for sin/sinh,  minus
the leading term E.  This can be expressed as

M = E(1-ecc) - ecc * E * remaining_terms( -E^2)    (elliptical)
M = E(1-ecc) - ecc * E * remaining_terms(  E^2)    (hyperbolic)
*/

static double remaining_terms( const double ival)
{
   double rval = 0., z = 1;
   const double tolerance = 1e-30;
   int i = 2;

   do
      {
      z *= ival / (double)( i * (i + 1));
      rval += z;
      i += 2;
      } while( fabs( z) > tolerance);
   return( rval);
}

int DLL_FUNC calc_classical_elements( ELEMENTS *elem, const double *r,
                             const double t, const int ref)
{
   const double *v = r + 3;
   const double r_dot_v = dot_product( r, v);
   const double dist = vector3_length( r);
   const double v2 = dot_product( v, v);
   double inv_major_axis = 2. / dist - v2 / elem->gm;
   double h0, n0, tval;
   double h[3], e[3], ecc2;
   double ecc;
   int i;

   assert( elem->gm != 0.);
   vector_cross_product( h, r, v);
   n0 = h[0] * h[0] + h[1] * h[1];
   h0 = n0 + h[2] * h[2];
   assert( dist > 0.);
   assert( v2 > 0.);     /* elements are undefined if the object is at rest */
   assert( h0 > 0.);     /* or if its velocity vector runs through the sun  */
   n0 = sqrt( n0);
   h0 = sqrt( h0);

                        /* See Danby,  p 204-206,  for much of this: */
   if( ref & 1)
      {
      if( !n0)    /* orbit is in xy plane;  asc node is undefined;  make */
         elem->asc_node = 0.;  /* arbitrary choice h[0] = 0, h[1] = -epsilon */
      else
         elem->asc_node = atan2( h[0], -h[1]);
      elem->incl = asine( n0 / h0);
      if( h[2] < 0.)                   /* retrograde orbit */
         elem->incl = PI - elem->incl;
      }
   vector_cross_product( e, v, h);
   for( i = 0; i < 3; i++)
      e[i] = e[i] / elem->gm - r[i] / dist;
   tval = dot_product( e, h) / h0;     /* "flatten" e vector into the rv */
   for( i = 0; i < 3; i++)             /* plane to avoid roundoff; see   */
      e[i] -= h[i] * tval;             /* above comments                 */
   ecc2 = dot_product( e, e);
   if( fabs( ecc2 - 1.) < 1.e-14)      /* avoid roundoff issues w/nearly */
      ecc2 = 1.;                       /* parabolic orbits               */
   elem->minor_to_major = sqrt( fabs( 1. - ecc2));
   ecc = elem->ecc = sqrt( ecc2);

   if( !ecc)                     /* for purely circular orbits,  e is */
      {                          /* arbitrary in the orbit plane; choose */
      for( i = 0; i < 3; i++)    /* r normalized                         */
         e[i] = r[i] / dist;
      }
   else                           /* ...and if it's not circular,  */
      for( i = 0; i < 3; i++)     /* normalize e:                  */
         e[i] /= ecc;
   if( ecc < .9)
      elem->q = (1. - ecc) / inv_major_axis;
   else        /* at eccentricities near one,  the above suffers  */
      {        /* a loss of precision problem,  and we switch to: */
      const double gm_over_h0 = elem->gm / h0;
/*    const double perihelion_speed = gm_over_h0 +
                   sqrt( fabs( gm_over_h0 * gm_over_h0 - inv_major_axis * elem->gm));
*/    const double perihelion_speed = gm_over_h0 *
               (1. + sqrt( 1. - inv_major_axis * h0 * h0 / elem->gm));

      assert( h0 != 0.);
      assert( gm_over_h0 != 0.);
      assert( isfinite( inv_major_axis));
      assert( isfinite( gm_over_h0));
      assert( isfinite( perihelion_speed));
      assert( perihelion_speed != 0.);
      elem->q = h0 / perihelion_speed;
      assert( elem->q != 0.);          /* For q=0,  nothing is defined */
      inv_major_axis = (1. - ecc) / elem->q;
      }
   assert( elem->q != 0.);          /* For q=0,  nothing is defined */
   assert( elem->q > 0.);

   if( inv_major_axis)
      {
      elem->major_axis = 1. / inv_major_axis;
      elem->t0 = elem->major_axis * sqrt( fabs( elem->major_axis) / elem->gm);
      }

   vector_cross_product( elem->sideways, h, e);
   if( ref & 1)
      {
      double cos_arg_per;

      if( n0)
         cos_arg_per = (h[0] * e[1] - h[1] * e[0]) / n0;
      else
         cos_arg_per = e[0];
      if( cos_arg_per < .7 && cos_arg_per > -.7)
         elem->arg_per = acos( cos_arg_per);
      else
         {
         double sin_arg_per;

         if( n0)
            sin_arg_per = (e[0] * h[0] * h[2] + e[1] * h[1] * h[2] - e[2] * n0 * n0)
                                            / (n0 * h0);
         else
            sin_arg_per = e[1] * h[2] / h0;

         elem->arg_per = fabs( asin( sin_arg_per));
         if( cos_arg_per < 0.)
            elem->arg_per = PI - elem->arg_per;
         }
      if( e[2] < 0.)
         elem->arg_per = PI + PI - elem->arg_per;
      }

   if( inv_major_axis && elem->minor_to_major)
      {
      const bool is_nearly_parabolic = (ecc > .99999 && ecc < 1.00001);
      const double r_cos_true_anom = dot_product( r, e);
      const double r_sin_true_anom = dot_product( r, elem->sideways) / h0;
      const double sin_E = r_sin_true_anom * inv_major_axis
                                        / elem->minor_to_major;

      assert( elem->minor_to_major);
      assert( isfinite( ecc));
      assert( isfinite( h0));
      assert( isfinite( r_cos_true_anom));
      assert( isfinite( r_sin_true_anom));
      assert( isfinite( sin_E));
      if( inv_major_axis > 0.)          /* parabolic case */
         {
         const double cos_E = r_cos_true_anom * inv_major_axis + ecc;
         const double ecc_anom = atan2( sin_E, cos_E);

         assert( isfinite( cos_E));
         assert( isfinite( ecc_anom));
         if( is_nearly_parabolic)
            elem->mean_anomaly = ecc_anom * (1 - ecc)
               - ecc * ecc_anom * remaining_terms( -ecc_anom * ecc_anom);
         else
            elem->mean_anomaly = ecc_anom - ecc * sin_E;
         assert( isfinite( elem->mean_anomaly));
         elem->perih_time = t - elem->mean_anomaly * elem->t0;
         }
      else                             /* hyperbolic case */
         {
         const double ecc_anom = asinh( sin_E);

         if( is_nearly_parabolic)
            elem->mean_anomaly = ecc_anom * (1 - ecc)
               - ecc * ecc_anom * remaining_terms( ecc_anom * ecc_anom);
         else
            elem->mean_anomaly = ecc_anom - ecc * sin_E;
         assert( isfinite( elem->mean_anomaly));
         assert( elem->t0 <= 0.);
         elem->perih_time = t - elem->mean_anomaly * fabs( elem->t0);
         h0 = -h0;
         }
      }
   else              /* parabolic case */
      {
      double tau;

      tau = sqrt( dist / elem->q - 1.);
      if( r_dot_v < 0.)
         tau = -tau;
      elem->w0 = (3. / SQRT_2) / (elem->q * sqrt( elem->q / elem->gm));
/*    elem->perih_time = t - tau * (tau * tau / 3. + 1) *                   */
/*                                elem->q * sqrt( 2. * elem->q / elem->gm); */
      elem->perih_time = t - tau * (tau * tau / 3. + 1) * 3. / elem->w0;
      }

         /* At this point,  elem->sideways has length h0.  Make it a unit vect: */
   for( i = 0; i < 3; i++)
      {
      elem->perih_vec[i] = e[i];
      elem->sideways[i] /= h0;
      }
   elem->angular_momentum = h0;
   return( 0);
}
