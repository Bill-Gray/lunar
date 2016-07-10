/* relativi.cpp: shows a simple-minded way to add GR to orbital computations

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "watdefs.h"
#include "comets.h"
#include "lunar.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

int load_vsop_data( void);
int integrate_orbit( ELEMENTS *elem, const double jd_from, const double jd_to,
                              const int n_steps);

static void add_relativistic_accel( double *accel, const double *posnvel)
{
   double h[3], h_squared, r_squared, accel_mag;
   const double c = AU_PER_DAY;           /* speed of light in AU per day */
   const double alpha = 3. / (c * c);        /* Solar GM handled elsewhere */
   int i;

   h[0] = posnvel[1] * posnvel[5] - posnvel[2] * posnvel[4];
   h[1] = posnvel[2] * posnvel[3] - posnvel[0] * posnvel[5];
   h[2] = posnvel[0] * posnvel[4] - posnvel[1] * posnvel[3];
   h_squared = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
   r_squared = posnvel[0] * posnvel[0] + posnvel[1] * posnvel[1]
                                       + posnvel[2] * posnvel[2];
   accel_mag = -alpha * h_squared / (r_squared * r_squared * sqrt( r_squared));
   for( i = 0; i < 3; i++)
      accel[i] += accel_mag * posnvel[i];
}

static void set_differential_acceleration( const double *posnvel,
                      const double *posnvel_2, double *accel)
{
   double p_squared = 0., r_squared = 0.;
   double pfactor, rfactor;
   int i;

   for( i = 0; i < 3; i++)
      {
      p_squared += posnvel[i] * posnvel[i];
      r_squared += posnvel_2[i] * posnvel_2[i];
      }
               /* someday,  I'll do it right;  for the nonce,  do it quick: */
   pfactor = 1. / (p_squared * sqrt( p_squared));
   rfactor = 1. / (r_squared * sqrt( r_squared));
   for( i = 0; i < 3; i++)
      accel[i] = pfactor * posnvel[i] - rfactor * posnvel_2[i];
   add_relativistic_accel( accel, posnvel_2);
}

char *vsop_data;

static void compute_perturber( int perturber_no, double jd,
                             double *perturber_loc)
{
   const double j2000 = 2451545.;
   double loc[15];

   compute_planet( vsop_data, perturber_no, (jd - j2000) / 36525., loc);
   memcpy( perturber_loc, loc + 12, 3 * sizeof( double));
}

static int n_perturbers = 9;

#ifdef PERTURBER_2
int perturber_2_callback( const double *posnvel, const double *perturber_loc);
#endif

static int compute_derivatives( const double jd, ELEMENTS *elems,
               const double *delta, double *derivs)
{
   double accel[3], posnvel[6], posnvel_2[6];
   static double relative_mass[11] = { 1.,
         1.66e-7, 2.448e-6,            /* mercury, venus */
         2.98994e-6, 3.23e-7,          /* earth, mars */
         9.54791e-4, 2.85878e-4,       /* jupiter,  saturn */
         4.365792e-5, 5.1776e-5,       /* uranus, neptune */
         8.0e-9, 2.98994e-6 / 81.3  };       /* pluto, moon */
   int i;

   comet_posn_and_vel( elems, jd, posnvel, posnvel + 3);
   for( i = 0; i < 6; i++)
      posnvel_2[i] = posnvel[i] + delta[i];
   set_differential_acceleration( posnvel, posnvel_2, accel);
   for( i = 2; i < n_perturbers; i++)
      {
      double perturber_loc[3], diff[3], diff_squared = 0., dfactor;
      double radius_squared = 0., rfactor;
      int j;

      compute_perturber( i, jd, perturber_loc);
      for( j = 0; j < 3; j++)
         {
         diff[j] = perturber_loc[j] - posnvel_2[j];
         diff_squared += diff[j] * diff[j];
         radius_squared += perturber_loc[j] * perturber_loc[j];
         }
      dfactor = relative_mass[i] / (diff_squared * sqrt( diff_squared));
      rfactor = relative_mass[i] / (radius_squared * sqrt( radius_squared));
      for( j = 0; j < 3; j++)
         accel[j] += diff[j] * dfactor - perturber_loc[j] * rfactor;
#ifdef PERTURBER_2
      if( i == 3)          /* earth */
         perturber_2_callback( posnvel, perturber_loc);
#endif
      }
   for( i = 0; i < 3; i++)
      {
      derivs[i] = delta[i + 3];
      derivs[i + 3] = SOLAR_GM * accel[i];
      }
   return( 0);
}

static int take_step( const double jd, ELEMENTS *elems,
                const double *ival, double *ovals, double *errs,
                const int n_vals, const double step_size)
{
   double *ivals[7], *ivals_p[6];
   int i, j, k;
   const double bvals[27] = {2. / 9.,
            1. / 12., 1. / 4.,
            69. / 128., -243. / 128., 135. / 64.,
            -17. / 12., 27. / 4., -27. / 5., 16. / 15.,
            65. / 432., -5. / 16., 13 / 16., 4 / 27., 5. / 144.,
            47. / 450., 0., 12 / 25., 32. / 225., 1. / 30., 6. / 25.,
            -1. / 150., 0., .03, -16. / 75., -.05, .24};
   const double *bptr = bvals;
   const double avals[6] = { 0., 2. / 9., 1./3., .75, 1., 5./6. };

   ivals[1] = (double *)calloc( 12 * n_vals, sizeof( double));
   if( !ivals[1])
      return( -1);
   for( i = 0; i < 6; i++)
      {
      ivals[i + 1] = ivals[1] + i * n_vals;
      ivals_p[i] = ivals[1] + (i + 6) * n_vals;
      }

   compute_derivatives( jd, elems, ival, ivals_p[0]);

   for( j = 1; j < 7; j++)
      {
      for( i = 0; i < n_vals; i++)
         {
         double tval = 0.;

         for( k = 0; k < j; k++)
            tval += bptr[k] * ivals_p[k][i];
         ivals[j][i] = tval * step_size + ival[i];
         }
      bptr += j;
      if( j != 6)
         compute_derivatives( jd + step_size * avals[j], elems,
                                                  ivals[j], ivals_p[j]);
      }

   if( errs)
      for( i = 0; i < n_vals; i++)
         {
         double tval = 0.;

         for( k = 0; k < 6; k++)
            tval += bptr[k] * ivals_p[k][i];
         errs[i] = step_size * tval;
         }

   memcpy( ovals, ivals[6], n_vals * sizeof( double));
   free( ivals[1]);
   return( 0);
}

double total_change;

int integrate_orbit( ELEMENTS *elem, const double jd_from, const double jd_to,
                              const int n_steps)
{
   int i;
   double step_size = (jd_to - jd_from) / (double)n_steps;
   double delta[6], jd = jd_from, posnvel[6];

   for( i = 0; i < 6; i++)
      delta[i] = 0.;
   for( i = 0; i < n_steps; i++, jd += step_size)
      {
      double new_delta[6];

      take_step( jd, elem, delta, new_delta, NULL, 6, step_size);
      memcpy( delta, new_delta, 6 * sizeof( double));
      }
   comet_posn_and_vel( elem, jd_to, posnvel, posnvel + 3);
   total_change = sqrt( delta[0] * delta[0] +
                        delta[1] * delta[1] + delta[2] * delta[2]);
   for( i = 0; i < 6; i++)
      posnvel[i] += delta[i];
   elem->epoch = jd_to;
   elem->gm = SOLAR_GM;
   calc_classical_elements( elem, posnvel, jd_to, 1);
   return( 0);
}

int load_vsop_data( void)
{
   FILE *ifile = fopen( "vsop.bin", "rb");
   const unsigned vsop_size = 60874u;

   vsop_data = NULL;
   if( ifile)
      {
      vsop_data = (char *)calloc( vsop_size, 1);
      if( vsop_data)
         fread( vsop_data, 1, vsop_size, ifile);
      fclose( ifile);
      }
   return( ifile && vsop_data ? 0 : -1);
}

#ifdef TEST_CODE
static void show_elems( const ELEMENTS *elem)
{
   printf( "JD %.2lf         Mean anom %9.5lf\n", elem->epoch,
                                    elem->mean_anomaly * 180. / PI);
   printf( "Asc node: %9.5lf     ", elem->asc_node * 180. / PI);
   printf( "Arg per: %9.5lf     ", elem->arg_per * 180. / PI);
   printf( "Incl: %9.5lf\n", elem->incl * 180. / PI);
   printf( "Ecc: %.8lf     ", elem->ecc);
   printf( "Major axis: %11.8lf    ", elem->major_axis);
   printf( "Perihelion dist: %11.8lf\n", elem->q);
   printf( "Perihelion date: %.6lf\n", elem->perih_time);
}

#ifndef __linux
#include <conio.h>
#else
   #define _getch getchar
#endif
#include <time.h>

int main( const int argc, const char **argv)
{
   ELEMENTS elem;
   int c, n_repeats = 1;
   time_t t0 = time( NULL), t1;
   double step0 = 10.;

   if( argc == 2)                /* shut off perturbing objects */
      n_perturbers = 2;
   elem.epoch = 2450000.;
   elem.mean_anomaly = 10. * PI / 180.;
   elem.arg_per = 20. * PI / 180.;
   elem.asc_node = 30. * PI / 180.;
   elem.incl = 7. * PI / 180.;
   elem.ecc = .20563175;
   elem.major_axis = .38710353;
   elem.q = elem.major_axis * (1. - elem.ecc);
   derive_quantities( &elem, SOLAR_GM);
   elem.perih_time = elem.epoch - elem.mean_anomaly * elem.t0;
   elem.angular_momentum = sqrt( SOLAR_GM * elem.major_axis);
   elem.angular_momentum *= sqrt( 1. - elem.ecc * elem.ecc);

   if( load_vsop_data( ))
      {
      printf( "VSOP.BIN not loaded!\n");
      return( -1);
      }
   show_elems( &elem);
   while( (c = _getch( )) != 27)
      {
      if( c >= 'a')
         {
         int i, n_steps = 1 << (c - 'a');
         double delta_t = step0 * n_steps;

         for( i = 0; i < n_repeats; i++)
            integrate_orbit( &elem, elem.epoch, elem.epoch + delta_t, n_steps);
         t1 = time( NULL);
         if( t1 != t0)
            {
            double posnvel[4];

            show_elems( &elem);
            comet_posn_and_vel( &elem, elem.epoch, posnvel, NULL);
            printf( "   %11.8lf %11.8lf %11.8lf\n", posnvel[0],
                                        posnvel[1], posnvel[2]);
            printf( "Precession: %lf arcsec/century\n",
                  (elem.arg_per * 180. / PI - 20.) * 3600. * 36525.
                        / (elem.epoch - 2450000.));
            printf( "Net Precession: %lf arcsec/century\n",
                  ((elem.asc_node + elem.arg_per) * 180. / PI - 50.) * 3600. * 36525.
                        / (elem.epoch - 2450000.));
            printf( "Total change (magnitude of delta): %lf\n", total_change);
            t0 = t1;
            }
         }
      else if( c >= '0' && c <= '9')
         {
         n_repeats = 1 << (c - '0');
         printf( "n_repeats reset to %d\n", n_repeats);
         }
      else if( c >= 'A' && c <= 'Z')
         {
         step0 = (double)( 1L << (c - 'A'));
         printf( "step0 reset to %f\n", step0);
         }
      }
   return( 0);
}
#endif
