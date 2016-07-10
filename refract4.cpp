/* refract4.cpp: functions for very precise refraction computations

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
         All references are to the _Explanatory Supplement to the
         Astronomical Almanac_,  pages 141-143

   The following code basically follows the algorithm set forth in
the above reference.  I did modify it a little by not exactly using
Simpson's rule for integration in its "usual" form.  Instead,  the
code uses Simpson's rule over a given interval,  and compares the
result to that which you would get from the trapezoidal rule
(integral evaluated using just the endpoints.)  If the difference is
greater than a certain tolerance,  the region is split in half and
the total_refraction() function recurses.

   The benefit of this is that the region is recursively subdivided,
with most of the evaluation being done in the places where the
function is changing most rapidly.  I'm sure it's not an original
idea,  but it does seem to improve performance a bit.

   Another benefit comes into play when the iterative procedure at
(3.281-5) is done,  to determine the value of 'r' corresponding to
a given value of 'z'.  With this "subdivision" technique,  we're
looking at an intermediate value of z (between the two endpoints),
and an excellent starting point is the intermediate value of r.
This better initial value in the iteration saves us a few steps.

   NOTE that situations with observed_alt < 0. (i.e.,  looking
"below" the horizon,  zenith angle greater than 90 degrees) or
height_in_meters > 11000 (i.e.,  above the troposphere) are not
handled properly and will usually return garbage results.
'observed_alt < 0.' can happen,  for example, if you're on top of
Mauna Kea looking out over the Pacific.  If you hand the routine
something like observed_alt = -55 (i.e.,  looking downward into
the earth),  it will gamely figure out what the refraction would
be if the earth's atmosphere continued downward,  gradually
increasing in pressure.  In this mathematical model,  the index
of refraction would increase greatly,  so it would be like looking
downward into a mirrored surface;  the ray of light would bend
_up_,  eventually exiting at your level at observed_alt = +55.
Thus,  for observed_alt = -55,  the refraction is about 110 degrees.

   Also, the curvature of the earth isn't handled;  the refraction
is of the sort you'd get with a flat earth.

   All of these problems could be fixed,  but it would require some
tinkering with the _Explanatory Supplement_ algorithm. */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
            /* Constants from p 142, (3.281-3): */
#define R   8314.36
#define Md    28.966
#define Mw    18.016
#define delta 18.36
#define re 6378120.
                /* ht = height of troposphere,  in meters */
#define ht 11000.
                /* hs = height of stratosphere,  in meters */
#define hs 80000.
#define alpha .0065
#define rt (re + ht)
#define rs (re + hs)

#define REFRACT struct refract
#define LOCALS  struct locals

REFRACT
   {
   double temp_0, n0_r0_sin_z0;
   double r0, c2, gamma, c6, c7, c8, c9, nt;
   };

LOCALS
   {
   double z, r, n, dn_dr, integrand;
   };

// double minimum_refractive_altitude = 10000000;

static void compute_refractive_index( const REFRACT *ref, LOCALS *loc,
                      const int is_troposphere)
{

// if( minimum_refractive_altitude > loc->r - re)
//    minimum_refractive_altitude = loc->r - re;
   if( is_troposphere)
      {                                   /* (3.281-6) */
      const double temp_fraction = 1. - alpha * (loc->r - ref->r0) / ref->temp_0;
      const double gamma_term = pow( temp_fraction, ref->gamma - 2.);
      const double delta_term = pow( temp_fraction, delta - 2.);

//    if( temp_fraction < 0.)
//       printf( "!!! r = %lf\n", loc->r - re);
      loc->n = 1. + temp_fraction * (ref->c6 * gamma_term - ref->c7 * delta_term);
      loc->dn_dr = -ref->c8 * gamma_term + ref->c9 * delta_term;
      }
   else                                   /* (3.281-7) */
      {
      const double temp_t = ref->temp_0 - alpha * (rt - ref->r0);
      const double exp_term =
                (ref->nt - 1.) * exp( -ref->c2 * (loc->r - rt) / temp_t);

      loc->n = 1. + exp_term;
      loc->dn_dr = -(ref->c2 / temp_t) * exp_term;
      }
}

static void compute_integrand( LOCALS *loc)
{
   const double r_dn_dr = loc->r * loc->dn_dr;

   loc->integrand = r_dn_dr / (loc->n + r_dn_dr);   /* (3.281-8) */
}

static double total_refraction( const REFRACT *ref,
        const LOCALS *l1, const LOCALS *l2, const int is_troposphere,
        const int recursion_depth)
{
   LOCALS mid;
   double change;
   const double iteration_limit = 1.;     /* get 'r' within a meter */
   const double integration_tolerance = .0001;
   const int max_recursion_depth = 20;

   mid.z = (l1->z + l2->z) * .5;
   mid.r = (l1->r + l2->r) * .5;
   do
      {
      compute_refractive_index( ref, &mid, is_troposphere);
      change = -(mid.n * mid.r - ref->n0_r0_sin_z0 / sin( mid.z));
      change /= mid.n + mid.r * mid.dn_dr;           /* (3.281-5) */
      mid.r += change;
      }
      while( fabs( change) > iteration_limit);

   compute_integrand( &mid);
            /* Compute difference between a Simpson's rule integration */
            /* and a trapezoidal one... */
   change = 2. * mid.integrand - (l1->integrand + l2->integrand);
// if( recursion_depth >= max_recursion_depth)
//    printf( "!!! recursed too deep\n");
            /* ...and if it's too great,  recurse with each half: */
   if( fabs( change) > integration_tolerance && recursion_depth < max_recursion_depth)
      return( total_refraction( ref, l1, &mid, is_troposphere, recursion_depth + 1)
            + total_refraction( ref, &mid, l2, is_troposphere, recursion_depth + 1));
   else
      {                    /* Simpson's rule is good enough: */
      const double h = (l2->z - l1->z) / 6.;

      return( h * (4. * mid.integrand + l1->integrand + l2->integrand));
      }
}

double DLL_FUNC integrated_refraction( const double latitude,
                  const double observed_alt, const double wavelength_microns,
                  const double height_in_meters, const double rel_humid_pct,
                  const double temp_kelvins, const double pressure_mb)
{
   const double g_bar = 9.784 * (1. - .0026 * cos( 2. * latitude)
                        - 2.8e-7 * height_in_meters);      /* (3.281-4) */
   const double Pw0 = rel_humid_pct * pow( temp_kelvins / 247.1, delta) / 100.;
   const double l2 = wavelength_microns * wavelength_microns;
   const double A = (273.15e-6 / 1013.25)
                        * (287.607 + 1.6288 / l2 + .0136 / (l2 * l2));
   REFRACT ref;
   double c5, rval;
   LOCALS l0, lt, ls;

   ref.r0 = re + height_in_meters;
   ref.c2 = g_bar * Md / R;
   ref.gamma = ref.c2 / alpha;    /* = C3 */
   c5 = Pw0 * (1. - Mw / Md) * ref.gamma / (delta - ref.gamma);
   ref.c6 = A * (pressure_mb + c5) / temp_kelvins;
   ref.c7 = (A * c5 + 11.2684e-6 * Pw0) / temp_kelvins;
   ref.c8 = alpha * (ref.gamma - 1.) * ref.c6 / temp_kelvins;
   ref.c9 = alpha * (delta - 1.) * ref.c7 / temp_kelvins;
   ref.temp_0 = temp_kelvins;

   l0.r = ref.r0;
   compute_refractive_index( &ref, &l0, 1);
   l0.z = PI / 2. - observed_alt;
   compute_integrand( &l0);

   lt.r = rt;
   compute_refractive_index( &ref, &lt, 1);
   lt.z = asin( l0.n * ref.r0 * sin( l0.z) / (lt.n * rt));     /* (3.281-9) */
   compute_integrand( &lt);

   ref.nt = lt.n;
   ref.n0_r0_sin_z0 = l0.n * l0.r * sin( l0.z);
   rval = total_refraction( &ref, &l0, &lt, 1, 0);

         /* Now for the stratospheric portion... we need to recompute */
         /* dn/dr at the tropopause, because there's a discontinuity  */
         /* in the derivative as r crosses that point;  which also    */
         /* means we gotta recompute the integrand at r = rt.  So:    */
   compute_refractive_index( &ref, &lt, 0);
   compute_integrand( &lt);

   ls.r = rs;
   compute_refractive_index( &ref, &ls, 0);
   ls.z = asin( ls.n * ref.r0 * sin( l0.z) / (ls.n * rs));
   compute_integrand( &ls);

   return( rval + total_refraction( &ref, &lt, &ls, 0, 0));
}

/* reverse_integrated_refraction( ) attempts to find the inverse of the
above function,  using Newton-Raphson integration.  See above comments on
the limitations of the refraction code (i.e.,  can't handle heights above
the troposphere,  nor negative refracted_alt).  We start out by using the
"primitive" reverse refraction formula (see refract.cpp).  Each call to
integrated_refraction( ) is relatively expensive,  and the hope is that
use of the "primitive" function will get us close enough to the right
answer to save us an iteration or two.  */

double DLL_FUNC reverse_integrated_refraction( const double latitude,
                  const double refracted_alt, const double wavelength_microns,
                  const double height_in_meters, const double rel_humid_pct,
                  const double temp_kelvins, const double pressure_mb)
{
            /* start out with an initial "primitive" guess: */
   double x1 = refracted_alt + reverse_refraction( refracted_alt)
                     * (pressure_mb / 1010.) * (283. / temp_kelvins);
   double y1 = refracted_alt - x1 + integrated_refraction( latitude, x1,
                           wavelength_microns, height_in_meters,
                           rel_humid_pct, temp_kelvins, pressure_mb);
   double x2 = x1 + y1;
   double y2 = refracted_alt - x2 + integrated_refraction( latitude, x2,
                           wavelength_microns, height_in_meters,
                           rel_humid_pct, temp_kelvins, pressure_mb);

   const double tolerance = .1 * PI / (180. * 3600.);     /* .1 arcsec */
   int n_iterations = 0;

// printf( "x1 = %lf; y1 = %lf\n", x1 * 180. / PI, y1 * 180. / PI);
// printf( "x2 = %lf; y2 = %lf\n", x2 * 180. / PI, y2 * 180. / PI);
   do
      {
      const double x3 = x2 - (x1 - x2) * y2 / (y1 - y2);

      if( fabs( x3 - x2) < tolerance)
         return( x3 - refracted_alt);
      else
         {
         const double y3 = refracted_alt - x3
                        + integrated_refraction( latitude, x3,
                           wavelength_microns, height_in_meters,
                           rel_humid_pct, temp_kelvins, pressure_mb);

//       printf( "x3 = %lf; y3 = %lf (%d)\n", x3 * 180. / PI,
//                  y3 * 180. / PI, n_iterations);
         n_iterations++;
         x1 = x2;
         y1 = y2;
         x2 = x3;
         y2 = y3;
         }
      }
   while( n_iterations < 10);
   return( x2 - refracted_alt);
}
