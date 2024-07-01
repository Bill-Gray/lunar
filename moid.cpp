/* moid.cpp: computes MOID (Minimum Orbital Intersection Distance)
between two orbits.

Copyright (C) 2018, Project Pluto

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

/* The MOID search works as follows.  We take two orbits as input,
and rotate both so that one of them (usually the one with the
smaller eccentricity) has its major axis running through the x-axis
and minor axis parallel to the y-axis.  (The MOID is,  of course,
unchanged by rotation.)

   We take advantage of the fact that computing the point on an
ellipse that is closest to a given point is a solved problem. There
are iterative techniques and an analytical method that involves
solving a quartic polynomial.  The latter method is used here;  see
the point_to_ellipse() function in mpc_code.cpp.

   Next,  we hunt around in the second orbit,  computing the minimum
distance between various points in that orbit and the first orbit.
This proceeds along a grid in true anomaly with spacing of,  at
most,  'moid_step' radians,  with various tricks used to figure out
what range(s) of true anomaly can be eliminated from consideration.
As we step along the second orbit,  we may find that we have
bracketed a "minimum minimum".  If so,  we can refine that minimum
and find a true MOID.  (Though there may be up to four minima,  so
our search isn't over yet,  and it is theoretically possible that
two minima may exist between points in the search grid.)

   After writing this,  I found a very similar approach discussed at

https://www.researchgate.net/publication/325922027_On_the_minimum_orbital_intersection_distance_computation_A_new_effective_method
*/

#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "watdefs.h"
#include "brentmin.h"
#include "comets.h"
#include "afuncs.h"
#include "mpc_func.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define J2000 2451545.0

static void fill_matrix( double mat[3][3], const ELEMENTS *elem)
{
   memcpy( mat[0], elem->perih_vec, 3 * sizeof( double));
   memcpy( mat[1], elem->sideways, 3 * sizeof( double));
         /* mat[2] is the cross-product of mat[0] & mat[1]: */
   vector_cross_product( mat[2], mat[0], mat[1]);
}

static double compute_posn_and_vel( const ELEMENTS *elem,
            const double true_anom, const double matrix[3][3],
            double *posn, double *vel)
{
   const double cos_true_anom = cos( true_anom);
   const double sin_true_anom = sin( true_anom);
   const double denom = 1. + elem->ecc * cos_true_anom;
   double true_r = elem->q * (1. + elem->ecc) / denom;
   double x, y;
   int i;

   if( true_r > 1000. || true_r < 0.)
      true_r = 1000.;
   x = true_r * cos_true_anom;
   y = true_r * sin_true_anom;
   for( i = 0; i < 3; i++)
      posn[i] = x * matrix[0][i] + y * matrix[1][i];
   if( vel)
      {
      const double dx_dtheta = -y / denom;
      const double dy_dtheta = (x + elem->ecc * true_r) / denom;
      const double vel_at_perihelion =    /* in AU/day */
                   GAUSS_K * sqrt( (1. + elem->ecc) / elem->q);
      const double dtheta_dt = vel_at_perihelion * elem->q / (true_r * true_r);
      const double dx_dt = dx_dtheta * dtheta_dt;
      const double dy_dt = dy_dtheta * dtheta_dt;

      for( i = 0; i < 3; i++)
         vel[i] = dx_dt * matrix[0][i] + dy_dt * matrix[1][i];
      }
   return( true_r);
}

static double true_anomaly_to_eccentric( const double true_anom,
                                         const double ecc)
{
   const double r = (1. - ecc * ecc) / (1. + ecc * cos( true_anom));
   const double x = r * cos( true_anom) + ecc;
   const double y = r * sin( true_anom) / sqrt( 1. - ecc * ecc);
   const double ecc_anom = PI + atan2( -y, -x);

   return( ecc_anom);
}

typedef struct
{
   double xform_matrix[3][3];
   double elem1_b, lat, r2;
   const ELEMENTS *elem1, *elem2;
   bool compute_obj_1_data;
   moid_data_t *mdata;
} internal_moid_t;

static double true_anom_from_r( const ELEMENTS *elem, const double r)
{
   double rval = 0.;

   if( elem->ecc)
      {
      const double cos_true_anom =
                      (elem->q * (1. + elem->ecc) - r) / (r * elem->ecc);

      if( cos_true_anom >= 1.)
         rval = 0.;
      else if( cos_true_anom <= -1.)
         rval = PI;
      else
         rval = acos( cos_true_anom);
      }
   return( rval);
}

static double find_point_moid_2( internal_moid_t *iptr, const double true_anomaly2)
{
   double vect2[3], dist, dist_squared;
   double x, y;

   iptr->r2 = compute_posn_and_vel( iptr->elem2, true_anomaly2,
                 iptr->xform_matrix, vect2, NULL);
   x = vect2[0] + iptr->elem1->q * iptr->elem1->ecc / (1. - iptr->elem1->ecc);
   y = vect2[1];   /* above shifts origin from focus to center of ellipse */
   iptr->lat = point_to_ellipse( iptr->elem1->major_axis, iptr->elem1_b,
                x, y, &dist);
   dist_squared = dist * dist + vect2[2] * vect2[2];
   if( iptr->compute_obj_1_data)
      {
      x = vect2[0] - dist * cos( iptr->lat);
      y = vect2[1] - dist * sin( iptr->lat);
      iptr->mdata->obj1_true_anom = atan2( y, x);
//    printf( "Dist %f; loc %f, %f; true anom %f/%f\n", dist, x, y,
//             iptr->mdata->obj1_true_anom * 180. / PI,
//             true_anomaly2 * 180. / PI);
      }
   return( dist_squared);
}

/* Ensure that we check for minima at least once per five degrees in
true anomaly of elem2's orbit : */

static double moid_step = 5. * PI / 180.;

double DLL_FUNC find_moid_full( const ELEMENTS *elem1, const ELEMENTS *elem2, moid_data_t *mdata)
{
   double mat1[3][3], mat2[3][3];
   internal_moid_t idata;
   moid_data_t mdata2;
   double least_dist_squared = 10000.;
   int i, j, n_steps, pass;
   double true_anomaly2, dist_squared, min_true2 = 0.;
   double min_true_anom = 0., max_true_anom = PI, dist;
   const double Q1 = elem1->major_axis * 2. - elem1->q;
   const double Q2 = elem2->major_axis * 2. - elem2->q;

   assert( elem1->ecc < 1.);
   fill_matrix( mat1, elem1);
   fill_matrix( mat2, elem2);
   for( i = 0; i < 3; i++)
      for( j = 0; j < 3; j++)
         idata.xform_matrix[i][j] = dot_product( mat1[j], mat2[i]);
   idata.mdata = (mdata ? mdata : &mdata2);
   idata.compute_obj_1_data = false;
   idata.elem1_b = elem1->major_axis * elem1->minor_to_major;
   idata.elem1 = elem1;
   idata.elem2 = elem2;
   true_anomaly2 = true_anom_from_r( elem2, elem1->major_axis);
// printf( "Nominal MOID true anom %f\n", true_anomaly2 * 180. / PI);
   least_dist_squared = find_point_moid_2( &idata, true_anomaly2);
   if( true_anomaly2 && true_anomaly2 != PI)
      {
      dist_squared = find_point_moid_2( &idata, -true_anomaly2);
      if( least_dist_squared > dist_squared)
         least_dist_squared = dist_squared;
      }
   if( elem2->ecc < 1.)
   {
      double ivect[3], arg_per;

//    printf( "starting distance constraint %f\n", sqrt( least_dist_squared));
      ivect[0] = -idata.xform_matrix[2][1];   /* a vector perpendicular to zaxis, */
      ivect[1] =  idata.xform_matrix[2][0];    /* lying in the x/y plane */
      ivect[2] = 0.;
      arg_per = -atan2( dot_product( ivect, idata.xform_matrix[1]),
                        dot_product( ivect, idata.xform_matrix[0]));
//    printf( " Arg per %f\n", arg_per * 180. / PI);
      for( pass = 0; pass < 2; pass++)
         {                     /* check ascending & descending nodes */
         dist_squared = find_point_moid_2( &idata, pass * PI - arg_per);
//       printf( "dist at node : %f\n", sqrt( dist_squared));
         if( dist_squared < least_dist_squared)
            least_dist_squared = dist_squared;
         }
   }
   for( pass = 0; pass < 2; pass++)
      {
      double x[3], y[3];

      dist = sqrt( least_dist_squared);
//    printf( "Curr distance constraint %f\n", dist);
      if( elem2->ecc >= 1. || Q2 > Q1 + dist)
         {
         max_true_anom = true_anom_from_r( elem2, Q1 + dist);
//       printf( "New max true anom %f\n", max_true_anom * 180. / PI);
         }
      if( elem2->q < elem1->q - dist)
         {
         min_true_anom = true_anom_from_r( elem2, elem1->q - dist);
//       printf( "New min true anom %f\n", min_true_anom * 180. / PI);
         }
      n_steps = (int)( (max_true_anom - min_true_anom) / moid_step) + 1;
      for( i = 0; i < 3; i++)
         x[i] = y[i] = 0.;
      for( i = -1; i <= n_steps + 1; i++)
         {
         true_anomaly2 = (max_true_anom - min_true_anom) * (double)i / (double)n_steps
                              + min_true_anom;
         if( pass)
            true_anomaly2 = -true_anomaly2;
         dist_squared = find_point_moid_2( &idata, true_anomaly2);
//       printf( "%f: %f %c\n", true_anomaly2 * 180. / PI, sqrt( dist_squared),
//                (least_dist_squared > dist_squared) ? '*' : ' ');
         if( least_dist_squared > dist_squared)
            {
            min_true2 = true_anomaly2;
            least_dist_squared = dist_squared;
            }
         x[0] = x[1];   y[0] = y[1];
         x[1] = x[2];   y[1] = y[2];
         x[2] = true_anomaly2;
         y[2] = dist_squared;
         if( i > 0 && y[1] < y[0] && y[1] < y[2])
            {                 /* we've got a MOID bracketed;  use Brent  */
            brent_min_t b;    /* algorithm to find the 'true' minimum */

            brent_min_init( &b, x[0], y[0], x[1], y[1], x[2], y[2]);
            b.tolerance  = .000001 * PI / 180.;
            b.ytolerance  = .0000001;
            while( b.step_type)
               {
               const double new_true = brent_min_next( &b);

               assert( b.n_iterations < 200);
               dist_squared = find_point_moid_2( &idata, new_true);
               if( least_dist_squared > dist_squared)
                  {
                  min_true2 = new_true;
                  least_dist_squared = dist_squared;
                  }
               if( b.step_type)
                  brent_min_add( &b, dist_squared);
               }
            }
         }
      }
   {
   double vdiff[3], vel1[3], vel2[3], unused_posn[3], ecc_anom, mean_anom;

   idata.compute_obj_1_data = true;
   find_point_moid_2( &idata, min_true2);
   compute_posn_and_vel( elem1, idata.mdata->obj1_true_anom, mat1, unused_posn, vel1);
// printf( "Obj1 posn: %f %f %f\n", unused_posn[0], unused_posn[1], unused_posn[2]);
// printf( "Obj1 vel : %f %f %f\n", vel1[0] * AU_IN_KM / seconds_per_day,
//                                  vel1[1] * AU_IN_KM / seconds_per_day,
//                                  vel1[2] * AU_IN_KM / seconds_per_day);
   compute_posn_and_vel( elem2, min_true2,            mat2, unused_posn, vel2);
// printf( "Obj2 posn: %f %f %f\n", unused_posn[0], unused_posn[1], unused_posn[2]);
// printf( "Obj2 vel : %f %f %f\n", vel2[0] * AU_IN_KM / seconds_per_day,
//                                  vel2[1] * AU_IN_KM / seconds_per_day,
//                                  vel2[2] * AU_IN_KM / seconds_per_day);
   for( i = 0; i < 3; i++)
      vdiff[i] = vel2[i] - vel1[i];
   idata.mdata->barbee_speed = vector3_length( vdiff);
   ecc_anom = true_anomaly_to_eccentric( idata.mdata->obj1_true_anom, idata.elem1->ecc);
   mean_anom = ecc_anom - idata.elem1->ecc * sin( ecc_anom);
   idata.mdata->jd1 = idata.elem1->perih_time + mean_anom * idata.elem1->t0;
   }

   return( sqrt( least_dist_squared));
}

static inline double centralize_angle( double ang)
{
   ang = fmod( ang, PI + PI);
   if( ang < -PI)
      ang += PI + PI;
   else if( ang > PI)
      ang -= PI + PI;
   return( ang);
}

#define N_PLANET_ELEMS 15
#define N_PLANET_RATES 9

int setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                             const double t_cen)
{
/* Planet elems taken straight from https://ssd.jpl.nasa.gov/txt/p_elem_t1.txt
See https://ssd.jpl.nasa.gov/?planet_pos for a discussion.  For MOID finding,
the longitude is irrelevant.  I left it in on the assumption that we might
want it someday. Asteroid elements are from BC-405,  epoch 2451535.0.  The
"longitudes" are actually mean anomalies,  and the "LonPer"s are actually
arguments of perihelion;  the code corrects for that last,  and the mean
anomaly/longitude isn't used (yet).   */

static const double planet_elem[N_PLANET_ELEMS * 6] = {
         /*     a          eccent    inclin    AscNode   LonPer      Longit */
/* Merc */   0.38709927,  0.20563593,  7.00497902,   48.33076593,  77.45779628,   252.25032350,
/* Venu */   0.72333566,  0.00677672,  3.39467605,   76.67984255, 131.60246718,   181.97909950,
/* EMB  */   1.00000261,  0.01671123, -0.00001531,    0.0       , 102.93768193,   100.46457166,
/* Mars */   1.52371034,  0.09339410,  1.84969142,   49.55953891, -23.94362959,    -4.55343205,
/* Jupi */   5.20288700,  0.04838624,  1.30439695,  100.47390909,  14.72847983,    34.39644051,
/* Satu */   9.53667594,  0.05386179,  2.48599187,  113.66242448,  92.59887831,    49.95424423,
/* Uran */  19.18916464,  0.04725744,  0.77263783,   74.01692503, 170.95427630,   313.23810451,
/* Nept */  30.06992276,  0.00859048,  1.77004347,  131.78422574,  44.96476227,   -55.12002969,
/* Plut */  39.48211675,  0.24882730, 17.14001206,  110.30393684, 224.06891629,   238.92903833,
/*  (1) */  2.7664603,  0.0783638, 10.583360,  80.494464,  73.921341,   4.036019,
/*  (2) */  2.7723257,  0.2296435, 34.846130, 173.197757, 310.264059, 350.826074,
/*  (4) */  2.3615363,  0.0900245,  7.133918, 103.951631, 149.589094, 338.305822,
/* (29) */  2.5543838,  0.0722511,  6.102741, 356.567840,  62.015715,  20.150301,
/* (16) */  2.9204983,  0.1382234,  3.093382, 150.465894, 229.122381, 333.613957,
/* (15) */  2.6437135,  0.1862108, 11.747399, 293.516504,  96.956836, 104.873024 };

         /* Rates are in same units,  but per century. */
static const double planet_elem_rate[N_PLANET_RATES * 6] = {
/* Merc */   0.00000037,  0.00001906, -0.00594749,   -0.12534081,   0.16047689,  149472.67411175,
/* Venu */   0.00000390, -0.00004107, -0.00078890,   -0.27769418,   0.00268329, 58517.81538729,
/* EMB  */   0.00000562, -0.00004392, -0.01294668,    0.0       ,   0.32327364, 35999.37244981,
/* Mars */   0.00001847,  0.00007882, -0.00813131,   -0.29257343,   0.44441088, 19140.30268499,
/* Jupi */  -0.00011607, -0.00013253, -0.00183714,    0.20469106,   0.21252668,  3034.74612775,
/* Satu */  -0.00125060, -0.00050991,  0.00193609,   -0.28867794,  -0.41897216,  1222.49362201,
/* Uran */  -0.00196176, -0.00004397, -0.00242939,    0.04240589,   0.40805281,   428.48202785,
/* Nept */   0.00026291,  0.00005105,  0.00035372,   -0.00508664,  -0.32241464,   218.45945325,
/* Plut */  -0.00031596,  0.00005170,  0.00004818,   -0.01183482,  -0.04062942,   145.20780515 };

   const double *pdata = planet_elem + (planet_idx - 1) * 6;
   const double *rate_data = planet_elem_rate + (planet_idx - 1) * 6;
   double elem_array[6];
   int i;

   if( planet_idx > N_PLANET_ELEMS || planet_idx < 1)
      return( -1);
   for( i = 0; i < 6; i++)
      {
      elem_array[i] = pdata[i];
      if( planet_idx <= N_PLANET_RATES)
         elem_array[i] += rate_data[i] * t_cen;
      }
   for( i = 2; i < 6; i++)
      elem_array[i] = centralize_angle( elem_array[i] * PI / 180.);
   memset( elem, 0, sizeof( ELEMENTS));
   elem->ecc = elem_array[1];
   elem->q = (1. - elem->ecc) * elem_array[0];
   elem->incl = elem_array[2];
   elem->asc_node = elem_array[3];
               /* For planets,  the longitude of perihelion is given. */
               /* For asteroids,  it's the argument of perihelion.    */
   if( planet_idx < 9)
      elem->arg_per = elem_array[4] - elem_array[3];
   else
      elem->arg_per = elem_array[4];
   elem->mean_anomaly = centralize_angle( elem_array[5] - elem->asc_node - elem->arg_per);
      /* l = (100.46435 + (129597740.63 / 3600.) * t_cen) * PI / 180.; */
   derive_quantities( elem, SOLAR_GM);
   elem->epoch = J2000 + t_cen * 36525.;
   elem->perih_time = elem->epoch - elem->mean_anomaly * elem->t0;
// printf( "epoch %f; t0 %f; mean anom %f\n", elem->perih_time, elem->t0,
//          elem->mean_anomaly * 180. / PI);
   return( 0);
}

