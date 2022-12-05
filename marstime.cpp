/* marstime.cpp: convert Martian "standard" time to/from TT
     (Temps Terrestienne,  Earth-based atomic time)

Copyright (C) 2016, Project Pluto

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

/*  Underlying algorithms copied from :                                     */
/*     http://www.giss.nasa.gov/tools/mars24/help/algorithm.html            */
/* Routines to compute "Mars Coordinated Time" (MCT),  the Martian          */
/* equivalent of UTC,  and "Mars True Solar Time (MTST) at Airy" (Airy is   */
/* the Martian equivalent of the Greenwich meridian),  for a given JDT.     */
/* The "reverse" code to convert MTST to TT was added by me (Bill Gray).    */
/*    A test case,  from the above URL:  if run with JDT=2451549.50074,     */
/* one should get:                                                          */
/*                                                                          */
/* pbs = 0.001418; a_fms = 272.744861; v_minus_m = 4.441908                 */
/* MTC = 44795.999760 (23:59:39.281); eot = -0.014410                       */
/* LTST at Airy: 23:38:54.247                                               */
/*                                                                          */
/*    The "recovered JD" should be equal to the input JDT of 2451549.50074; */
/* i.e.,  the time transformations should all be correctly reversed.        */

#include <math.h>
#include <stdlib.h>
#ifdef TEST_PROGRAM
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "stringex.h"
#endif

const double days_per_sol = 1.0274912517;
const double zero_sol_point = 44796.0 - 0.0009626;
const double zero_jd_point = 2451549.5;

double tt_to_mtc( const double jd);
double mtc_to_tt( const double mtc);
double mars_true_solar_minus_mean_solar_time( const double jd);
double mtst_at_airy_to_tt( const double mtst);

double tt_to_mtc( const double jd)
{
                                 /* eqn C-2: */
   return( (jd - zero_jd_point) / days_per_sol + zero_sol_point);
}

double mtc_to_tt( const double mtc)
{
                                 /* C-2 equation reversed: */
   return( (mtc - zero_sol_point) * days_per_sol + zero_jd_point);
}

/* "longitude_sun" = 0 degrees at the northern hemisphere vernal equinox;
= 90 degrees at summer solstice,  = 180 at autumnal equinox,  = 270 at
winter solstice. */

double mars_true_solar_minus_mean_solar_time( const double jd)
{
   const double jd_2000 = 2451545.0;  /* JD 2451545.0 = 1.5 Jan 2000 */
   const double t = jd - jd_2000;
   const double pi =
      3.1415926535897932384626433832795028841971693993751058209749445923;
   const double D2R = pi / 180.;
            /* equations B-1 & B-2: */
   const double mars_mean_anom = 19.3871 * D2R + .52402073 * D2R * t;
        /* 'a_fms' = 'angle of fictitious mean sun' */
   const double a_fms =         270.3871 * D2R + .524038496 * D2R * t;
#ifndef IGNORE_PERTURBERS
   const double tconst = 2. * pi / 365.25;
   static const double amplit[7] = { .0071 * D2R, .0057 * D2R,
                        .0039 * D2R, .0037 * D2R, .0021 * D2R,
                        .0020 * D2R, .0018 * D2R};
   static const double freq[7] = { tconst / 2.2353, tconst / 2.7543,
               tconst / 1.1177, tconst / 15.7866, tconst / 2.1354,
               tconst / 2.4694, tconst / 32.8493 };
   static const double phase[7] = { 49.409 * D2R, 168.173 * D2R,
               191.837 * D2R, 21.736 * D2R, 15.704 * D2R,
               95.528 * D2R, 49.095 * D2R };
#endif
            /* v_minus_m = true minus mean anomaly,  a.k.a. the */
            /* equation of the center:                          */
   double v_minus_m = (10.691 * D2R + 3e-7 * D2R * t)
                             * sin( mars_mean_anom)
               + 0.623 * D2R * sin( 2. * mars_mean_anom)
               + 0.050 * D2R * sin( 3. * mars_mean_anom)
               + 0.005 * D2R * sin( 4. * mars_mean_anom)
               + 0.0005 * D2R * sin( 5. * mars_mean_anom);
   double longitude_sun, eqn_of_time;
#ifndef IGNORE_PERTURBERS
   double pbs = 0.;
   int i;

   for( i = 0; i < 7; i++)       /* eqn B-3 */
      pbs += amplit[i] * cos( freq[i] * t + phase[i]);
   v_minus_m += pbs;
#endif
   longitude_sun = a_fms + v_minus_m;       /* eqn B-5 */
#ifdef TEST_PROGRAM
   printf( "a_fms = %f; v_minus_m = %f; longitude of sun = %f\n",
            a_fms / D2R, v_minus_m / D2R, longitude_sun / D2R);
#endif
   eqn_of_time = (2.861 / 360.) * sin( 2. * longitude_sun)   /* eqn C-1 */
               - (0.071 / 360.) * sin( 4. * longitude_sun)
               + (0.002 / 360.) * sin( 6. * longitude_sun)
               - v_minus_m / (2. * pi);
            /* Above equation of time = true - mean time,  in sols. */
   return( eqn_of_time);
}

/* Reversing MTST (Mars True Solar Time) to other systems is made      */
/* slightly tricky by the fact that the first step involves computing  */
/* the equation of time,  which takes TT as an input.  So we pretend   */
/* the input MTST is actually an MTC,  and compute a TT from it using  */
/* mtc_to_tt.  This gives us an "approx_tt" which may be up to an hour */
/* off (it basically is ignoring the Martian equation of time).        */
/*    However,  if we compute the Martian EOT using approx_tt,  we'll  */
/* get a passably correct EOT and can use it to compute a better TT.   */
/* And we can then compute the EOT using this better TT to get a still */
/* better TT.  These two iterations of computing the EOT are enough to */
/* get a "real" TT that's good to machine precision.                   */

double mtst_at_airy_to_tt( const double mtst)
{
   const double approx_tt = mtc_to_tt( mtst);
   double rval = approx_tt;
   int iter;

   for( iter = 2; iter; iter--)
      {
      const double eqn_of_time =
              mars_true_solar_minus_mean_solar_time( rval);

      rval = approx_tt - eqn_of_time * days_per_sol;
      }
   return( rval);
}

#ifdef TEST_PROGRAM
static void format_time( const double day, char *buff)
{
   const double time_of_day = day - floor( day);
   const int n_millisec = (int)( time_of_day * 24. * 60. * 60. * 1000. + .5);

   snprintf_err( buff, 13, "%02d:%02d:%02d.%03d",
            n_millisec / (60 * 60 * 1000),      /* hour */
            (n_millisec / (60 * 1000)) % 60,    /* minutes */
            (n_millisec / 1000) % 60,           /* seconds */
            n_millisec % 1000);                 /* millisec */
   assert( 12 == strlen( buff));
}

int main( const int argc, const char **argv)
{
   const double jd = (argc > 1 ? atof( argv[1]) : 2451549.50074);
   const double mtc = tt_to_mtc( jd);
   const double eot = mars_true_solar_minus_mean_solar_time( jd);
   const double ltst_at_airy = mtc + eot;
   char buff[80];

            /* Above equation of time = true - mean time. */
   format_time( mtc, buff);
   printf( "MTC = %f (%s); eot = %f\n", mtc, buff, eot);
   format_time( ltst_at_airy, buff);
   printf( "LTST at Airy: %s\n", buff);
   printf( "Recovered JD: %.8f\n", mtst_at_airy_to_tt( ltst_at_airy));
   if( argc > 2)                  /* West longitudes are positive */
      {
      const double lon = atof( argv[2]);
      const double ltst = ltst_at_airy - lon / 360.;
      const double lmst = mtc - lon / 360.;

      format_time( ltst, buff);
      printf( "LTST at loc: %s\n", buff);
      format_time( lmst, buff);
      printf( "LMST at loc: %s\n", buff);
      }
   return( 0);
}
#endif
