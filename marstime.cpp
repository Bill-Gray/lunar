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
/* Routines to compute "Mars Mean Standard Time" (MST),  the Martian        */
/* equivalent of TT,  and "Mars True Solar Time (MTST) at Airy" (Airy is    */
/* the Martian equivalent of the Greenwich meridian),  for a given JDT.     */
/* The "reverse" code to convert MTST to TT was added by me (Bill Gray).    */
/*    A test case,  from the above URL:  if run with MJDT=51549.00074       */
/* (the default value),  one should get:                                    */
/*                                                                          */
/* pbs = 0.001418; a_fms = 272.74566102; v-m = 4.44192663; Ls = 277.1875876 */
/* MST = 44795.999758 (23:59:39.057); eot = -0.014410                       */
/* LTST at Airy: 23:38:53.998                                               */
/*                                                                          */
/*    The "recovered MJD" should be equal to the input MJDT of 51549.00074; */
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
const double zero_mjd_point = 51549.0;

double tt_to_mst( const double mjd);
double mst_to_tt( const double mst);
double mars_true_solar_minus_mean_solar_time( const double mjd);
double mtst_at_airy_to_tt( const double mtst);

double tt_to_mst( const double mjd)
{
                                 /* eqn C-2: */
   return( (mjd - zero_mjd_point) / days_per_sol + zero_sol_point);
}

double mst_to_tt( const double mst)
{
                                 /* C-2 equation reversed: */
   return( (mst - zero_sol_point) * days_per_sol + zero_mjd_point);
}

/* "longitude_sun" = 0 degrees at the northern hemisphere vernal equinox;
= 90 degrees at summer solstice,  = 180 at autumnal equinox,  = 270 at
winter solstice. */

double mars_true_solar_minus_mean_solar_time( const double mjd)
{
   const double mjd_2000 = 51544.5;  /* JD 2451545.0 = 1.5 Jan 2000 */
   const double t = mjd - mjd_2000;
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
   printf( "pbs = %f; ", pbs / D2R);
   printf( "a_fms = %.8f; v-m = %.8f; Ls = %.8f\n",
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
/* the input MTST is actually an MST (Mars Mean Solar Time),  and      */
/* compute a TT from it using mst_to_tt.  This gives us an "approx_tt" */
/* which may be up to an hour off (it basically is ignoring the        */
/* Martian equation of time).                                          */
/*    However,  if we compute the Martian EOT using approx_tt,  we'll  */
/* get a passably correct EOT and can use it to compute a better TT.   */
/* And we can then compute the EOT using this better TT to get a still */
/* better TT.  Continued iterations would converge linearly;  we use   */
/* just two iterations,  then compute a more 'exact' TT using Aitken's */
/* delta-squared sequence method for accelerated convergence.  The     */
/* result is good to machine precision.                                */

double mtst_at_airy_to_tt( const double mtst)
{
   const double approx_tt = mst_to_tt( mtst);
   double rval[3], d1, d2, x;
   size_t i;

   rval[0] = approx_tt;
   for( i = 1; i < 3; i++)
      {
      const double eqn_of_time =
              mars_true_solar_minus_mean_solar_time( rval[i - 1]);
      const double delta = approx_tt - eqn_of_time * days_per_sol - rval[i - 1];

      rval[i] = rval[i - 1] + delta;
      printf( "iter %d : %.10lf\n", (int)i, rval[i]);
      }
   d1 = rval[2] - rval[1];             /* Aitken delta-squared */
   d2 = rval[1] - rval[0];
   x = rval[2] - d1 * d1 / (d1 - d2);
   printf( "Aitken %.10f\n", x);
   return( x);
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
   const double mjd = (argc > 1 ? atof( argv[1]) : 51549.00074);
   const double mst = tt_to_mst( mjd);
   const double eot = mars_true_solar_minus_mean_solar_time( mjd);
   const double ltst_at_airy = mst + eot;
   char buff[80];

            /* Above equation of time = true - mean time. */
   format_time( mst, buff);
   printf( "MST = %f (%s); eot = %f\n", mst, buff, eot);
   format_time( ltst_at_airy, buff);
   printf( "LTST at Airy: %s\n", buff);
   printf( "Recovered MJD: %.10f\n", mtst_at_airy_to_tt( ltst_at_airy));
   if( argc > 2)                  /* West longitudes are positive */
      {
      const double lon = atof( argv[2]);
      const double ltst = ltst_at_airy - lon / 360.;
      const double lmst = mst - lon / 360.;

      format_time( ltst, buff);
      printf( "LTST at loc: %s\n", buff);
      format_time( lmst, buff);
      printf( "LMST at loc: %s\n", buff);
      }
   return( 0);
}
#endif
