/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* 2012 Jul 22:  (BJG) Debugging statements are now shown only if
one defines DEBUGGING_STATEMENTS.  The Vincenty routine will compute
a back-azimuth only if you pass it a non-NULL pointer to 'back_azimuth'.

   The following code implements three methods for computing the
distance between two lat/lon points on the earth.  The first formula
uses a method due to H. Andoyer,  from _Annuaire du Bureau des
Longitudes pour 1950_, p. 145.  It takes the earth's flattening into
account,  resulting in a relative error of the order of the square
of the earth's flattening (i.e.,  about one part in 90000).  I got
the formula "second hand" from Jean Meeus' _Astronomical
Algorithms_,  pp 80-82.

   The second method just uses plain spherical trigonometry,  as if
the earth really was round.  This can cause an error of about .5
percent.

   The third method is due to Thaddeus Vincenty,  and is documented at
http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf .  An implementation in
JavaScript is available at
http://www.movable-type.co.uk/scripts/latlong-vincenty.html , and at
http://www.ga.gov.au/geodesy/datums/vincenty_inverse.jsp .  The method
is supposed to be accurate to within .5 millimeters.  I've made some
modifications to the basic method,  as described below.

   The actual error of Andoyer's method,  compared to the presumably
near-error-free method of Vincenty,  can be as much as 15 km in some
_very_ unusual near-antipodal cases.  (Try,  for example, N 30, E 0
to S 29.99, E or W 179.99.)  However,  except for places within a
degree or so of the antipodes,  the error seems to be about a kilometer
at most (one part in 20000),  and therefore not usually a problem.
*/

double earth_dist( const double lat1, const double lon1,
                   const double lat2, const double lon2,
                   const double flattening);
double spherical_earth_dist( const double lat1, const double lon1,
                             const double lat2, const double lon2);
double vincenty_earth_dist( const double lat1, const double lon1,
                            const double lat2, const double lon2,
                            const double flat, double *azimuth,
                            double *back_azimuth);
void vincenty_direct( const double lat1, const double lon1,
                            double *lat2, double *lon2,
                            const double flat, const double azimuth,
                            const double dist);

double earth_dist( const double lat1, const double lon1,
                   const double lat2, const double lon2,
                   const double flattening)
{
   const double f = (lat1 + lat2) / 2.;
   const double g = (lat1 - lat2) / 2.;
   const double lambda = (lon1 - lon2) / 2.;
   const double sin_g = sin( g), cos_g = cos( g);
   const double sin_f = sin( f), cos_f = cos( f);
   const double sin_lambda = sin( lambda), cos_lambda = cos( lambda);
   const double sin_g2 = sin_g * sin_g, cos_g2 = cos_g * cos_g;
   const double sin_f2 = sin_f * sin_f, cos_f2 = cos_f * cos_f;
   const double sin_lambda2 = sin_lambda * sin_lambda;
   const double cos_lambda2 = cos_lambda * cos_lambda;
   const double s = sin_g2 * cos_lambda2 + cos_f2 * sin_lambda2;
   const double c = cos_g2 * cos_lambda2 + sin_f2 * sin_lambda2;
   const double omega = atan( sqrt( s / c));
   const double r = sqrt( s * c) / omega;
   const double d = 2. * omega;
   const double h1 = (3. * r - 1.) / (2. * c);
   const double h2 = (3. * r + 1.) / (2. * s);
   const double rval = d * (1. + flattening *
                       (h1 * sin_f2 * cos_g2 - h2 * cos_f2 * sin_g2));

   return( rval);
}

/* This uses the 'haversine' method for solving the spherical triangle.
The commonly-used method of computing distance over a spherical earth is

(1) cos( dist) = sin(lat1) sin(lat2) + cos( lat1) cos(lat2) cos( delta_lon)

   The problem is that if the points are close to each other,  cos( dist)
is very close to 1.  If they're near the antipodes,  it's very close to
-1.  So any roundoff may push you outside the domain of the cosine function,
or at the very least produce a horrendously magnified error when you take
the arccosine.

   The haversine formula doesn't lose precision for the two-point-near-
each-other case.  It still would have trouble near the antipodes,
where you'd get hav(dist) close to 1.  I say "would" because if hav(dist)
is greater than 0.9,  you'll see that we turn it into 180 degrees minus
a short-distance problem.

   Note that hav(x) = sin^2(x / 2) = (1-cos(x)) / 2,  and that
x = 2 * asin( sqrt( hav(x))).

The usual 'law of haversines' for a spherical triangle with sides
a, b, c and angles A, B, C is

(2) hav(a) = hav(b - c) + sin( b) sin( c) hav( A)

   (this is actually a rearrangement of (1).)  However, sin( b) sin( c)
= 1/2 (cos( b-c) - cos( b+c)) = hav( b+c) - hav( b-c).  So (2) can
be rewritten as

(3) hav(a) = hav(b - c) + (hav(b + c) - hav( b - c)) hav( A)

   which is the form used below.  Note that while I've used it here to
determine a side given two sides and an included angle,  it could be
trivially rearranged if you wanted to get an angle from three sides.
You'd have to do similar "flip to the antipodes" tricks to handle the
case where hav( A) is close to 1,  i.e.,  A is nearly 180 degrees. */

const double pi = 3.14159265358979323846264338327950288419716939937510582;

double spherical_earth_dist( const double lat1, const double lon1,
                             const double lat2, const double lon2)
{
   double x, y, z;
   double hav_rval;

   x = sin( (lon1 - lon2) / 2);
   x *= x;                          /* x = hav( lon1 - lon2) */
   y = sin ((lat1 - lat2) / 2);
   y *= y;                          /* y = hav( lat1 - lat2) */

   /* Seem to be more precise : */
   z = cos ((lat1 + lat2) / 2);
   z*=z;                            /* z = hav( lat1 + lat2) */

   hav_rval = x * (z - y) + y;
   if( hav_rval > 0.9)              /* close to antipodes */
      return( pi - spherical_earth_dist( lat1, lon1,
                     -lat2, lon2 + pi));
   return 2 * asin( sqrt( hav_rval));
}

/* Some modifications were made in this implementation of the
Vincenty method:

  -- The "reduced latitudes",  u1 and u2,  are calculated slightly
differently to avoid domain errors at the poles:  with the original
formulation,  if a latitude happens to be _exactly_ +/-90 degrees,
you got a domain error in the tangent function.  If it rounded off to
a hair greater than +/-90 degrees, you got a reduced latitude at the
opposite pole.

  -- The core of the code is basically a root-finding algorithm,
trying to find the value of lambda for which delta_lambda (in the
original formulation,  lambda - lambdaprime) is zero.  If we are
far from the antipodes,  the root-finding proceeds briskly with
the original method (and slightly faster with the secant method
implemented below).  Closer to the antipodes,  though,  delta_lambda
becomes a much more complicated function,  and finding the actual
root requires some more care.  The original method could fail to
converge,  or enter a loop between two values.  Also, near the
antipodes,  there can be two solutions,  one leading to a
maximum-distance geodesic and the other to a minimum-distance
geodesic.  We want the former (corresponding to lambda < 180),
but the original code could get the latter (lambda > 180).

   To fix this,  I added some bracketing code:  we know that
delta_lambda is positive for lambda=0 and negative for lambda=180,
so that's our initial bracket,  narrowed down each time we compute
a new delta_lambda.  If our next step would take us outside that
bracket (or if too many steps have occurred, indicating slow
convergence),  we do a bisection step.

   Near the antipodes,  delta_lambda becomes ill-behaved,  and we
almost always go to the bisection code.  Elsewhere,  the newly-added
secant method makes for slightly faster convergence (though
convergence away from the antipodes was already pretty good).

   Also,  _really_ near the antipodes,  delta_lambda becomes
so poorly behaved that with mere 64-bit floats,  we can't get
delta < tolerance. That's why we break out of the bisection
code if it's clear that lambda isn't actually changing.  (Which
isn't a problem,  because we've determined lambda to sufficient
accuracy.)

   -- All this narrowed the "trouble zone" to the area where the
latitudes are of opposite sides,  and the longitudes are within
(flattening*180) degrees of being 180 degrees apart.  To evade this,
if we're within antipode_tol of this line segment,  the code just
bumps the second latitude by antipode_tol,  resulting in a possible
two-millimeter error.  Even much of this could be eliminated,  by
"adjusting" the result by the latitude correction multiplied by the
cosine of the azimuth... something for another day.

*/

static double big_a_poly( const double usquared)
{
   return( 1. + usquared *
       (4096 + usquared * (-768 + usquared * (320 - 175 * usquared))) / 16384.);
}

static double big_b_poly( const double usquared)
{
   return( usquared *
       (256. + usquared * (-128. + usquared * (74 - 47 * usquared))) / 1024.);
}

static double cvt_lat( const double b, const double lat)
{
   return( atan2( b * sin( lat), cos( lat)));;
}

double vincenty_earth_dist( const double lat1, const double lon1,
                            const double lat2, const double lon2,
                            const double flat, double *azimuth,
                            double *back_azimuth)
{
   const double b = 1 - flat;       /* if a = 1 */
            /* u1, u2 are "reduced latitudes" */
   const double u1 = cvt_lat( b, lat1);
   const double u2 = cvt_lat( b, lat2);
   const double cos_u1 = cos( u1), sin_u1 = sin( u1);
   const double cos_u2 = cos( u2), sin_u2 = sin( u2);
   double lambda, delta_lambda = 2.;
   double prev_lambda[2], prev_delta[2];
   double low_lambda = 0., high_lambda = pi;
   double dlon = fmod( lon2 - lon1, 2. * pi);
   const double antipode_tol = 3e-10;  /* corresponds to 2 mm */
   const double tolerance = 1e-12;     /* corresponds to .006 mm */
   double temp1, temp2, temp3;
   double sin_lambda, cos_lambda;
   double sin_o, cos_o, o, cos2_a, cos_2om, big_c;
   int iter = 0, flipped_azimuth = 0;
   const int max_iter = 100;

   if( dlon > pi)          /* Things are a hair simpler if 0 < dlon < 180. */
      dlon -= 2. * pi;     /* We note if the azimuths have to be flipped.  */
   else if( dlon < -pi)
      dlon += 2. * pi;
   if( dlon < 0.)
      {
      dlon = -dlon;
      flipped_azimuth = 1;
      }
                       /* If we're too close to a line near the antipodes, */
                       /* we try again slightly offset from that line:     */
   if( fabs( u1 + u2) < .9 * antipode_tol && dlon > pi * b)
      {
#ifdef DEBUGGING_STATEMENTS
      printf( "Using lat %.10f\n", (antipode_tol - lat1) * 180. / pi);
#endif
      return( vincenty_earth_dist( lat1, lon1, antipode_tol - lat1,
                                 lon2, flat, azimuth, back_azimuth));
      }
   if( lat1 == lat2 && dlon == 0.)      /* points are identical */
      {
      if( azimuth)
         *azimuth = 0.;
      if( back_azimuth)
         *back_azimuth = 0.;
      return( 0.);
      }
   lambda = dlon;
   prev_lambda[0] = prev_delta[0] = 0.;
       /* above line added solely to remove a compiler warning */
   do
      {
      double sin_a;

      sin_lambda = sin( lambda);
      cos_lambda = cos( lambda);
      temp1 = cos_u2 * sin_lambda;
      temp2 = cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda;
      cos_o = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda;
      sin_o = sqrt( temp1 * temp1 + temp2 * temp2);
      o = atan2( sin_o, cos_o);
      sin_a = cos_u1 * cos_u2 * sin_lambda / sin_o;
//    printf( "temp1 %.14f temp2 %.14f  sin_o %.14f cos_o %.14f sin_a %.14f\n",
//             temp1, temp2, sin_o, cos_o, sin_a);
      cos2_a = 1 - sin_a * sin_a;  /* trig identity */
      if( cos2_a == 0.)
         cos_2om = 0.;
      else
         cos_2om = cos_o - 2 * sin_u1 * sin_u2 / cos2_a;
      big_c = flat * cos2_a * (4. + flat * (4 - 3 * cos2_a)) / 16.;
      temp3 = 2. * cos_2om * cos_2om - 1.;
      delta_lambda = dlon - lambda + (1.-big_c) * flat * sin_a
                   * (o + big_c * sin_o * (cos_2om + big_c * cos_o * temp3));
//    printf( "Iter %d: lambda %.14f; delta %.14f\n", iter,
//                                  lambda * 180. / pi,
//                                  delta_lambda * 180. / pi);
      if( delta_lambda > 0.)     /* update our root brackets */
         low_lambda = lambda;
      else
         high_lambda = lambda;
      prev_lambda[1] = prev_lambda[0];
      prev_lambda[0] = lambda;       /* store previous lambda & delta_lambda */
      prev_delta[1] = prev_delta[0];   /* values for use in secant method */
      prev_delta[0] = delta_lambda;
      if( !iter)
         lambda += delta_lambda;
      else
         {
         const double dx = prev_lambda[1] - prev_lambda[0];
         const double dy = prev_delta[1] - prev_delta[0];

         if( dy != 0.)
            lambda -= delta_lambda * dx / dy;
#ifdef DEBUGGING_STATEMENTS
         else
            printf( "???? dy = 0\n");
#endif
         }
                  /* following basically means "if the proposed lambda    */
                  /* lies outside the bracket,  or if we've done too many */
                  /* steps and obviously aren't converging briskly,  then*/
                  /* skip the secant steps and do a bisection step."    */
      if( lambda < low_lambda || lambda > high_lambda
                     || (iter > 10 && !(iter % 3)))
         {
#ifdef DEBUGGING_STATEMENTS
         printf( "Gone outside brackets (%.14f)!  Doing a bisection\n",
                        lambda * 180. / pi);
         printf( "%.14f to %.14f (%.7lg)\n",
               low_lambda * 180. / pi, high_lambda * 180. / pi,
               (high_lambda - low_lambda) * 180. / pi);
#endif
         delta_lambda = (high_lambda - low_lambda) / 2.;
         lambda = low_lambda + delta_lambda;
         if( high_lambda == pi)
            lambda += delta_lambda * .8;
         if( lambda == prev_lambda[0])     /* if no change occurs, */
            {                              /* we've  "maxed out" and */
#ifdef DEBUGGING_STATEMENTS
            printf( "NO ACTUAL CHANGE\n"); /* should leave the loop  */
#endif
            delta_lambda = 0.;
            }
         }
      iter++;
      }
   while( fabs( delta_lambda) > tolerance && iter < max_iter);

   const double usquared = cos2_a * (1.-b*b) / (b*b);
   const double big_b = big_b_poly( usquared);
   const double delta_o = big_b * sin_o *
              (cos_2om + (big_b / 4.) * (cos_o * temp3
            - (big_b / 6.) * cos_2om * (-3 + 4 * sin_o * sin_o)
            * (-3. + 4. * cos_o * cos_o)));
   const double rval = b * big_a_poly( usquared) * (o - delta_o);
   if( azimuth)
      {
      *azimuth = atan2( temp1, temp2);
      if( *azimuth < 0.)
         *azimuth += 2. * pi;
      if( flipped_azimuth)
         *azimuth = 2. * pi - *azimuth;
      }
   if( back_azimuth)
      {
      *back_azimuth = atan2( cos_u1 * sin_lambda,
             cos_u1 * sin_u2 * cos_lambda - sin_u1 * cos_u2);
      if( *back_azimuth < 0.)
         *back_azimuth += 2. * pi;
      if( flipped_azimuth)
         *back_azimuth = 2. * pi - *back_azimuth;
      }
   return( rval);
}

void vincenty_direct( const double lat1, const double lon1,
                            double *lat2, double *lon2,
                            const double flat, const double azimuth,
                            const double dist)
{
   const double pi = 3.14159265358979323846264338327950288419716939937510582;
   const double b = 1 - flat;       /* if a = 1 */
            /* u1, u2 are "reduced latitudes" */
   const double u1 = cvt_lat( b, lat1);
   const double cos_u1 = cos( u1), sin_u1 = sin( u1);
   const double sin_az = sin( azimuth), cos_az = cos( azimuth);
   const double sigma1 = atan2( sin_u1, cos_az * cos_u1);
   const double sin_alpha = cos_u1 * sin_az;
   const double cos2_alpha = 1. - sin_alpha * sin_alpha;
   const double u_squared = cos2_alpha * (1 - b * b) / (b * b);
   const double big_a = big_a_poly( u_squared);
   const double big_b = big_b_poly( u_squared);
   double sigma = dist / (b * big_a);
   double sigma_prime = 2. * pi;
   double cos_2sigma_m, sin_sigma, cos_sigma;
   const double tolerance = 1e-12;     /* corresponds to .006 mm */

   while( fabs( sigma - sigma_prime) > tolerance)
      {
      double delta_sigma, cos_2sigma_m2;

      cos_2sigma_m = cos( 2. * sigma1 + sigma);
      cos_2sigma_m2 = cos_2sigma_m * cos_2sigma_m;
      sin_sigma = sin( sigma);
      delta_sigma = big_b * sin_sigma *
               (cos_2sigma_m + (big_b / 4.) * (cos( sigma) *
               (2. * cos_2sigma_m2 - 1.) - (big_b / 6.)
               * cos_2sigma_m * (4. * sin_sigma * sin_sigma - 3.) *
               (4. * cos_2sigma_m2 - 3.)));

      sigma_prime = sigma;
      sigma = dist / (b * big_a) + delta_sigma;
      }

   sin_sigma = sin( sigma);
   cos_sigma = cos( sigma);
   cos_2sigma_m = cos( 2. * sigma1 + sigma);
   const double tval = sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_az;
   *lat2 = atan2( sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_az,
            b * sqrt( sin_alpha * sin_alpha + tval * tval));
   double lambda = atan2( sin_sigma * sin_az,
            cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_az);
// double big_c = flat / (16. * cos2_alpha * (4. + flat * (4. - 3 * cos2_alpha)));
   double big_c = (flat / 16) * cos2_alpha * (4. + flat * (4. - 3 * cos2_alpha));
   double big_l = lambda - (1. - big_c) * flat * sin_alpha *
         (sigma + big_c * sin_sigma * (cos_2sigma_m + big_c * cos_sigma *
         (2. * cos_2sigma_m * cos_2sigma_m - 1.)));

   *lon2 = lon1 + big_l;
}

static void reset_max_diff( double *max_diff, const double a, const double b)
{
   if( *max_diff < fabs( a))
      *max_diff = fabs( a);
   if( *max_diff < fabs( b))
      *max_diff = fabs( b);
}

/* A few test cases from Vincenty's paper
(http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf).   'data' gives,  in
order:  lat1, lat2, delta_lon, az1, distance, az2.  The first test
case,  (a),  is for the Bessel ellipsoid,  so I skipped it.    */

static void run_test_cases( void)
{
   const double data[] =  {
        37. + 19. / 60. + 54.95367 / 3600., 26. + 7/60. + 42.83946 / 3600.,
                     41 + 28./60. + 35.50729 / 3600.,
        95. + 27. / 60. + 59.63089 / 3600., 4085966.703,   /* (b) */
        118 + 5/60. + 58.96161 / 3600.,

        35. + 16/60. + 11.24862 / 3600., 67. + 22/60. + 14.77638 / 3600.,
                     137. + 47./60. + 28.31435 / 3600.,
        15 + 44./60. + 23.74850 / 3600., 8084823.839,
        144.+55./60.+39.92147 / 3600.,

        1., -59. / 60. - 53.83076 / 3600.,
                  179. + 17. / 60. + 48.02997 / 3600.,
        89., 19960000., 91. + 6.11733 / 3600.,

        1., 1. + 1./60. + 15.18952/3600.,
               179. + 46./60. + 17.84244 / 3600.,
        4. + 59./60. + 59.99995/3600., 19780006.558,
              174. + 59./60. + 59.88481/3600.,
        0. };
   const double *tptr = data;
   const double pi = 3.14159265358979323846264338327950288419716939937510582;
   double max_diff = 0.;
   const double semimajor = 6378388.0;   /* Int'l spheroid */

   while( *tptr)
      {
      double lat, lon, dist, az1, az2;
      const double flattening = 1. / 297.;

      vincenty_direct( tptr[0] * pi / 180., 0., &lat, &lon, flattening,
                       tptr[3] * pi / 180., tptr[4] / semimajor);
      lat *= 180. / pi;
      lon *= 180. / pi;
      printf( "Result : %.9f %.9f (%.9f %.9f)\n", lat, lon,
                     lat - tptr[1], lon - tptr[2]);
      reset_max_diff( &max_diff, lat - tptr[1], lon - tptr[2]);
      vincenty_direct( tptr[1] * pi / 180., 0., &lat, &lon, flattening,
                  pi + tptr[5] * pi / 180., tptr[4] / semimajor);
      lat *= 180. / pi;
      lon *= 180. / pi;
      printf( "Reverse: %.9f %.9f (%.9f %.9f)\n", lat, lon,
                     lat - tptr[0], lon + tptr[2]);
      reset_max_diff( &max_diff, lat - tptr[0], lon + tptr[2]);
      dist = vincenty_earth_dist( tptr[0] * pi / 180., 0., tptr[1] * pi / 180.,
                                  tptr[2] * pi / 180., flattening, &az1, &az2);
      dist *= semimajor;
      az1 *= 180. / pi;
      az2 *= 180. / pi;
      printf( "Dist: %.4f (%.4f)\n", dist, dist - tptr[4]);
      reset_max_diff( &max_diff, (dist - tptr[4]) / semimajor, 0.);
      printf( "Az1, 2: %.9f (%.9f)  %.9f (%.9f)\n",
                  az1, az1 - tptr[3],
                  az2, az2 - tptr[5]);
      reset_max_diff( &max_diff, az1 - tptr[3], az2 - tptr[5]);
      tptr += 6;
      }
   printf( "Max difference : %.9f degrees = %.5f meters\n",
                 max_diff, max_diff * semimajor);
}

/* For the 'spherical earth' formula,  we use a radius of 6371 km (close to
the average of the polar and equatorial radii.)  For the 'flattened earth',
the equatorial radius of 6378.140 km is used. */

int main( const int argc, const char **argv)
{
   const double pi = 3.14159265358979323846264338327950288419716939937510582;
   const double flattening = 1. / 298.257223563;
   const double semimajor = 6378.137;
// const double flattening = 1. / 298.257;
// const double semimajor = 6378.14;
   double lat1, lon1;

   if( argc < 5)
      {
      printf( "usage: dist lat1 lon1 lat2 lon2\n");
      printf( "or:    dist lat1 lon1 dist(km) azim -d\n");
      run_test_cases( );
      return( 0);
      }
   lat1 = atof( argv[1]) * pi / 180.;
   lon1 = atof( argv[2]) * pi / 180.;
   if( argc == 6)
      {
      double lat2, lon2, back_azimuth;
      const double dist = atof( argv[3]) / semimajor;
      double azimuth = atof( argv[4]) * pi / 180.;

      vincenty_direct( lat1, lon1, &lat2, &lon2, flattening, azimuth, dist);
      printf( "Vincenty direct: %.10f %.10f\n",
               lat2 * 180. / pi, lon2 * 180. / pi);
      vincenty_earth_dist( lat1, lon1, lat2, lon2, flattening, &azimuth,
                                    &back_azimuth);
      printf( "Reverse az: %.10f\n", back_azimuth * 180. / pi);
      }
   else if( argc == 5)
      {
      double lat2 = atof( argv[3]) * pi / 180.;
      double lon2 = atof( argv[4]) * pi / 180.;
      double azimuth, dist, back_azimuth;

      const double d1 = spherical_earth_dist( lat1, lon1, lat2, lon2);
      const double d2 =           earth_dist( lat1, lon1, lat2, lon2, flattening);

      printf( "From lat=%f, lon=%f to lat=%f, lon=%f\n",
                           lat1 * 180. / pi, lon1 * 180. / pi,
                           lat2 * 180. / pi, lon2 * 180. / pi);
      printf( "Distance from 'round earth': %f\n", d1 * 6371.);
      printf( "Distance from H. Andoyer: %f\n", d2 * semimajor);
      dist = vincenty_earth_dist( lat1, lon1, lat2, lon2, flattening, &azimuth,
                                    &back_azimuth);
      printf( "Distance from Vincenty: %.6f\n", semimajor * dist);
      printf( "Azimuth %.10f (back az: %.10f)\n", azimuth * 180. / pi,
                                              back_azimuth * 180. / pi);

      lat2 = lon2 = -1.;
      vincenty_direct( lat1, lon1, &lat2, &lon2, flattening, azimuth, dist);
      printf( "Vincenty direct: %.10f %.10f\n",
               lat2 * 180. / pi, lon2 * 180. / pi);
      }
   return( 0);
}
