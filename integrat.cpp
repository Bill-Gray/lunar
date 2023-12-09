/* integrat.cpp: numerically integrates 'mpcorb.dat' to arbitrary epochs

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
#ifdef _MSC_VER
#include <conio.h>
#endif
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "watdefs.h"
#include "stringex.h"
#include "comets.h"
#include "lunar.h"
#include "date.h"
#include "afuncs.h"        /* for rotate_vector( ) proto */

#if defined( __has_include)
   #if !__has_include(<jpleph.h>)
   #error   \
        'jpleph.h' not found.  This project depends on the 'jpl_eph'\
        library.  See www.github.com/Bill-Gray/jpl_eph .\
        Clone that repository,  'make'  and 'make install' it.
#ifdef __GNUC__
   #include <stop_compiling_here>
         /* Above line suppresses cascading errors. */
#endif
#endif
#endif

#include "jpleph.h"

/* On some (non-Windows) system,  spreading the integration out to
multiple processes is possible.  It's handled a little oddly.  The number
of processes can be specified on the command line with the -z switch.
Things proceed unchanged at first,  until data for the three perturbing
asteroids (Ceres,  Pallas,  Vesta) have been read and their ephemerides
computed. _Then_ we fork,  with each process reading every n_processes
line,  and "chunk" files being created.  Thus,  if n_processes is 7,  the
first process will create a chunk file with asteroids 5, 12, 19, 26, ...
The second will process asteroids 6, 13, 20, ...

   Once all the processes complete,  the original process zippers the
results from the chunk files together and unlinks them.     */

#if defined( __linux) || defined( __unix__) || defined( __APPLE__)
   #define FORKING

   #include <unistd.h>     /* Symbolic Constants */
   #include <sys/types.h>  /* Primitive System Data Types */
   #include <errno.h>      /* Errors */
   #include <stdio.h>      /* Input/Output */
   #include <sys/wait.h>   /* Wait for Process Termination  */
            /* above basically allows for forking so we can */
            /* run different objects on different cores     */
   #include <sys/time.h>         /* these allow resource limiting */
   #include <sys/resource.h>     /* see '-r' command switch below */
#endif


#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

#define PERTURBERS_MERCURY_TO_NEPTUNE 0xff
#define PERTURBERS_PLUTO 0x100
#define PERTURBERS_MOON  0x200
#define PERTURBERS_PLUTO_AND_MOON (PERTURBERS_PLUTO | PERTURBERS_MOON)
/* Following macro isn't actually used...
   #define PERTURBERS_CERES_PALLAS_VESTA 0x1c00   */
#define N_PERTURBERS 13
         /* hash table sizes should be prime numbers: */
#define HASH_TABLE_SIZE 3000017

static int verbose = 0, n_steps_taken = 0, resync_freq = 50;
static int asteroid_perturber_number = -1;
static double *position_cache;
static int position_cache_size = 0;
static unsigned long perturber_mask = PERTURBERS_MERCURY_TO_NEPTUNE;
      /*  PERTURBERS_MERCURY_TO_NEPTUNE | PERTURBERS_CERES_PALLAS_VESTA; */

int integrate_orbit( ELEMENTS *elem, double jd, const double jd_end,
                              const double max_err, const double stepsize);

/* 28 Feb 2003:  modified heavily after getting an e-mail from Werner
Huget. See his e-mail and page 281 of the _Explanatory Supplement to the
Astronomical Almanac_ for details.  Basically,  computing the relativistic
acceleration in the simple manner I previously used led to significant
errors for Mercury,  and presumably for other objects orbiting close
to the Sun. */

static void add_relativistic_accel( double *accel, const double *posnvel)
{
   int i;
   const double c = AU_PER_DAY;           /* speed of light in AU per day */
   const double r_squared = posnvel[0] * posnvel[0] + posnvel[1] * posnvel[1]
                                                    + posnvel[2] * posnvel[2];
   const double v_squared = posnvel[3] * posnvel[3] + posnvel[4] * posnvel[4]
                                                    + posnvel[5] * posnvel[5];
   const double v_dot_r   = posnvel[0] * posnvel[3] + posnvel[1] * posnvel[4]
                                                    + posnvel[2] * posnvel[5];
   const double r = sqrt( r_squared), r_cubed_c_squared = r_squared * r * c * c;
   const double r_component =
                  (4. * SOLAR_GM / r - v_squared) / r_cubed_c_squared;
   const double v_component = 4. * v_dot_r / r_cubed_c_squared;

   for( i = 0; i < 3; i++)
      accel[i] += r_component * posnvel[i] + v_component * posnvel[i + 3];
}

static FILE *err_fopen( const char *filename, const char *permits)
{
   FILE *rval = fopen( filename, permits);

   if( !rval)
      {
      printf( "Couldn't open file '%s'\n", filename);
#ifdef _MSC_VER
      printf( "Hit any key:\n");
      getch( );
#endif
      exit( -1);
      }
   return( rval);
}

#ifdef FORKING
static char *chunk_filename( char *filename, const int chunk_number, const int type)
{
   const char *format = (type ? "ephem%d.ugh" : "chunk%d.ugh");

   assert( chunk_number >= 0 && chunk_number <= 99);
   snprintf_err( filename, 12, format, chunk_number);
   return( filename);
}
#endif

static void set_differential_acceleration( const double *posnvel,
                      const double *delta, double *accel)
{
   double p_squared = 0., r_squared = 0.;
   double pfactor, rfactor, posnvel_2[6];
   int i;

   for( i = 0; i < 6; i++)
      posnvel_2[i] = posnvel[i] + delta[i];
   for( i = 0; i < 3; i++)
      {
      p_squared += posnvel[i] * posnvel[i];
      r_squared += posnvel_2[i] * posnvel_2[i];
      }
               /* someday,  I'll do it right;  for the nonce,  do it quick: */
               /* SEE: \useless\smalldif.cpp */
   pfactor = 1. / (p_squared * sqrt( p_squared));
   rfactor = 1. / (r_squared * sqrt( r_squared));
   for( i = 0; i < 3; i++)
      accel[i] = pfactor * posnvel[i] - rfactor * posnvel_2[i];
   add_relativistic_accel( accel, posnvel_2);
}

int load_vsop_data( void);

static char *vsop_data;
static void *jpl_ephemeris;

/* perturber_loc[0, 1, 2] = heliocentric ecliptic coords */

static void compute_perturber( int perturber_no, double jd,
                             double *perturber_loc)
{
   if( vsop_data)
      {
      const double j2000 = 2451545.;
      double loc[15];

      compute_planet( vsop_data, perturber_no, (jd - j2000) / 36525., loc);
      memcpy( perturber_loc, loc + 12, 3 * sizeof( double));
      }
   else
      {
      static double jd0 = -1, posns[11][6];
      int i;

      if( jd0 != jd)
         {
         int list[14];
         const double ratio = 1. + jpl_get_double( jpl_ephemeris,
                                          JPL_EPHEM_EARTH_MOON_RATIO);

         for( i = 0; i < 14; i++)
            list[i] = (i < 10);
         jpl_state( jpl_ephemeris, jd, list, posns, NULL, 0);
         for( i = 0; i < 3; ++i)
            {
            posns[2][i] -= posns[9][i] / ratio;
            posns[9][i] += posns[2][i];
            }
         jd0 = jd;
         for( i = 0; i < 10; i++)
            {
            const double sin_obliq_2000 = 0.397777155931913701597179975942380896684;
            const double cos_obliq_2000 = 0.917482062069181825744000384639406458043;
            double temp = posns[i][1] * cos_obliq_2000 + posns[i][2] * sin_obliq_2000;

            posns[i][2] = posns[i][2] * cos_obliq_2000 - posns[i][1] * sin_obliq_2000;
            posns[i][1] = temp;
            }
         }

      for( i = 0; i < 3; ++i)
         perturber_loc[i] = posns[perturber_no - 1][i];  /* - posns[10][i]; */
               /* rotate equatorial J2000.0 into ecliptical J2000: */
      }
}

static double *make_position_cache( double jd0, double jd_end,
                                 const double stepsize)
{
   double *rval, *tptr;
   int i, j, n_steps, step;
   const double max_jd = 2451545.0;    /* 2000 jan 1.5 */

   if( jd0 > jd_end)
      {
      const double tval = jd0;

      jd0 = jd_end;
      jd_end = tval;
      }
   if( jd0 > max_jd)    /* make sure cache goes back to at least 2000 */
      jd0 = max_jd;
   jd0 = floor( (jd0 - .5) / stepsize) * stepsize + 0.5;
         /* Make a few extra steps before & after the planned time range */
   n_steps = floor( (jd_end - jd0) / stepsize) + 5;
   jd0 -= stepsize * 2.;
   rval = (double *)calloc( 2 + (size_t)n_steps * N_PERTURBERS * 6 * 3,
                                           sizeof( double));
   tptr = rval + 2;
   if( !rval)
      {
      printf( "Ran out of memory!\n");
      exit( -1);
      }
   rval[0] = jd0;
   rval[1] = stepsize;
   position_cache_size = n_steps;
   for( step = 0; step < n_steps; step++)
      {
      for( j = 0; j < 6; j++)
         {
         const double avals[6] = { 0., 2. / 9., 1./3., .75, 1., 5./6. };

         for( i = 0; i < N_PERTURBERS; i++)
            {
            if( (i < 10) && ((perturber_mask >> i) & 1ul))
               compute_perturber( i + 1, jd0 + avals[j] * stepsize, tptr);
            else        /* put it far,  far away where it won't do anything: */
               tptr[0] = tptr[1] = tptr[2] = 1.e+8;
            tptr += 3;
            }
         }
      jd0 += stepsize;
#ifdef OBSOLETE_DEBUGGING_CODE
      while( step * 70 / n_steps > counter)
         {
         printf( "%d", counter % 10);
         counter++;
         }
#endif
      }
   printf( "\n");
   return( rval);
}

#define EARTH_MOON_RATIO 81.30056

static double relative_mass[14] = { 1.,
         1.660136795271931e-007,                /* mercury */
         2.447838339664545e-006,                /* venus */
         3.003489596331057e-006,                /* Earth */
         3.227151445053866e-007,                /* Mars */
         0.0009547919384243268,                 /* Jupiter */
         0.0002858859806661309,                 /* saturn */
         4.366244043351564e-005,                /* Uranus */
         5.151389020466116e-005,                /* Neptune */
         7.396449704142013e-009,                /* Pluto */
         3.003489596331057e-006 / EARTH_MOON_RATIO, /* Moon */
         4.7622e-10, 1.0775e-10, 1.3412e-10 };    /* Ceres,  Pallas, Vesta */

#define MERCURY_R   (2439.4 / AU_IN_KM)
#define VENUS_R     (6051. / AU_IN_KM)
#define EARTH_R     (6378.140 / AU_IN_KM)
#define MARS_R      (3397.0 / AU_IN_KM)
#define JUPITER_R   (71492. / AU_IN_KM)
#define SATURN_R    (60330. / AU_IN_KM)
#define URANUS_R    (25559. / AU_IN_KM)
#define NEPTUNE_R   (25225. / AU_IN_KM)
#define PLUTO_R     (1500. / AU_IN_KM)
#define MOON_R      (1748.2 / AU_IN_KM)

/* See 'runge.cpp' in Find_Orb for an explanation of this.  Basically,
it keeps accelerations from reaching infinity as an object passes through
a planet.  Integrate backward,  and 2018 LA,  2008 TC3,  and 2014 AA
will do exactly that,  and the integration step size can drop to zero.
The following code ramps acceleration _down_ as you approach the center
of a planet.  */

#define FUDGE_FACTOR   0.9

static double compute_accel_multiplier( double fraction)
{
   const double r0 = .8;  /* acceleration drops to zero at 80% of planet radius */
   double rval;

   assert( fraction >= 0. && fraction <= 1.);
   if( fraction < r0)
      rval = 0.;
   else
      {
      fraction = (fraction - r0) / (1. - r0);
      assert( fraction >= 0. && fraction <= 1.);
      rval =  fraction * fraction * (3. - 2. * fraction);
      }
   return( rval);
}

static int compute_derivatives( const double jd, ELEMENTS *elems,
               double *delta, double *derivs, double *posn_data)
{
   double accel[3], posnvel[6];
   int i;

   comet_posn_and_vel( elems, jd, posnvel, posnvel + 3);
   set_differential_acceleration( posnvel, delta, accel);
   for( i = 0; i < N_PERTURBERS; i++)       /* include perturbers */
      if( (perturber_mask >> i) & 1ul)
         {
         double perturber_loc[3], diff[3], diff_squared = 0., dfactor;
         double radius_squared = 0., rfactor, d, r;
         int j;
         static const double planet_radius[10] = {
                         MERCURY_R * FUDGE_FACTOR,
                        VENUS_R * FUDGE_FACTOR, EARTH_R * FUDGE_FACTOR,
                        MARS_R * FUDGE_FACTOR, JUPITER_R * FUDGE_FACTOR,
                        SATURN_R * FUDGE_FACTOR, URANUS_R * FUDGE_FACTOR,
                        NEPTUNE_R * FUDGE_FACTOR, PLUTO_R * FUDGE_FACTOR,
                        MOON_R * FUDGE_FACTOR };

         if( posn_data)
            memcpy( perturber_loc, posn_data + i * 3, 3 * sizeof( double));
         else
            if( i < 10)
               compute_perturber( i + 1, jd, perturber_loc);
            else
               perturber_loc[0] = perturber_loc[1] = perturber_loc[2] = 1.e+8;
         for( j = 0; j < 3; j++)
            {
            diff[j] = perturber_loc[j] - (posnvel[j] + delta[j]);
            diff_squared += diff[j] * diff[j];
            radius_squared += perturber_loc[j] * perturber_loc[j];
            }
         d = sqrt( diff_squared);
         r = sqrt( radius_squared);
         dfactor = relative_mass[i + 1] / (diff_squared * d);
         rfactor = relative_mass[i + 1] / (radius_squared * r);
         if( i < 10)
            {
            if( d < planet_radius[i])
               dfactor *= compute_accel_multiplier( d / planet_radius[i]);
            if( r < planet_radius[i])
               rfactor *= compute_accel_multiplier( r / planet_radius[i]);
            }
         for( j = 0; j < 3; j++)
            accel[j] += diff[j] * dfactor - perturber_loc[j] * rfactor;
         }

                      /* copy in Ceres,  Pallas, Vesta loc if needed: */
   if( posn_data && asteroid_perturber_number >= 0)
      memcpy( posn_data + asteroid_perturber_number * 3, posnvel,
                        3 * sizeof( double));
   for( i = 0; i < 3; i++)
      {
      derivs[i] = delta[i + 3];
      derivs[i + 3] = SOLAR_GM * accel[i];
      }
   return( 0);
}

#define N_VALUES 6
      /* i.e.,  a state vector consumes six values: x, y, z, vx, vy, vz */

static int take_step( const double jd, ELEMENTS *elems,
                double *ival, double *ovals, double *errs,
                double step_size)
{
   double *ivals[7], *ivals_p[6];
   double ivals_1_buff[12 * N_VALUES];
   double *posn_data = NULL;
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

   ivals[1] = ivals_1_buff;
   for( i = 0; i < 6; i++)
      {
      ivals[i + 1] = ivals[1] + i * N_VALUES;
      ivals_p[i] = ivals[1] + (i + 6) * N_VALUES;
      }

   if( fabs( step_size - position_cache[1]) < .000001)
      {
      int cache_loc = (int)floor( (jd - position_cache[0]) / step_size + .5);

      if( cache_loc >= 0 && cache_loc < position_cache_size)
         posn_data = position_cache + 2 + cache_loc * 6 * N_PERTURBERS * 3;
      }

   compute_derivatives( jd, elems, ival, ivals_p[0], posn_data);

   for( j = 1; j < 7; j++)
      {
      for( i = 0; i < N_VALUES; i++)
         {
         double tval = 0.;

         for( k = 0; k < j; k++)
            tval += bptr[k] * ivals_p[k][i];
         ivals[j][i] = tval * step_size + ival[i];
         }
      bptr += j;
      if( j != 6)
         compute_derivatives( jd + step_size * avals[j], elems,
                     ivals[j], ivals_p[j], posn_data ?
                     posn_data + j * N_PERTURBERS * 3 : NULL);
      }

   if( errs)
      for( i = 0; i < N_VALUES; i++)
         {
         double tval = 0.;

         for( k = 0; k < 6; k++)
            tval += bptr[k] * ivals_p[k][i];
         errs[i] = step_size * tval;
         }

   memcpy( ovals, ivals[6], N_VALUES * sizeof( double));
   n_steps_taken++;
   return( 0);
}

/* The following 'full_rk_step' integrates using the Runge-Kutta-Fehlberg
fifth-order integrator with automatic stepsize,  as described in J M A
Danby's _Fundamentals of Celestial Mechanics_,  second edition,  pages
297-299.  Basically,  the integration is done both to fourth and fifth
order.  The difference gives us an idea of the error for that step.  If
it is greater than some desired amount,  we can try again with a smaller
step size.

   After each step,  we recompute the step size:  if the previous step
resulted in a really low error,  we need to raise the step size,  and
if it caused a lot of error,  we decrease the step size.  */

static int full_rk_step( ELEMENTS *elems, double *ivals, double *ovals,
                double t0, double t1, double max_err)
{
   double step = t1 - t0;
   double errs[N_VALUES], new_vals[N_VALUES];
   int n_chickens = 0;

   memcpy( ovals, ivals, N_VALUES * sizeof( double));
   max_err *= max_err;
   while( t0 != t1)
      {
      double err_val = 0.;
      const double chicken_factor = .9;
      int i;

      take_step( t0, elems, ovals, new_vals, errs, step);
      for( i = 0; i < N_VALUES; i++)
         err_val += errs[i] * errs[i];
      if( err_val < max_err)   /* yeah,  it was a good step */
         {
         memcpy( ovals, new_vals, N_VALUES * sizeof( double));
         t0 += step;
         }
      else
         n_chickens++;
      step *= chicken_factor * exp( log( max_err / err_val) / 5.);
      if( t0 < t1)
         if( t0 + step > t1)
            step = t1 - t0;
      if( t1 < t0)
         if( t0 + step < t1)
            step = t1 - t0;
/*    if( err_val >= max_err)                             */
/*       printf( "Chickened out: new step %lf\n", step);  */
      }
   return( n_chickens);
}

/* Used for caching integrated positions (if n_xyzs is set to zero).
The cached positions are written to a file;  the file can then be
interpolated within to look for mutual close approaches.  This can
be done for asteroid mass measurements ("asteroid A came close to
more massive asteroid B;  measurements of how much A was perturbed
can determine B's mass") and to see what asteroids a spacecraft
might approach.         */

static int n_xyzs = -1;
static double *xyzs = NULL, ephem_start = 0., ephem_stepsize;

/* 'integrate_orbit' integrates the elements over the desired time span to
   the desired maximum error,  using the number of steps requested.  The
   orbit is broken up into that many steps,  and 'full_rk_step' is then
   called for each step.  The individual steps will probably be taken in
   one RKF step,  but if their errors prove to be too great,  they'll
   be broken into sub-steps.  See comments for the above code.

   The reason for this is speed.  Much of Integrat's time is spent in
   computing planetary positions.  If the steps fall on an evenly spaced
   grid,  the positions can be drawn from a precomputed array.  For the
   cases that break up into sub-steps,  planetary positions have to be
   computed "from scratch".  But with a suitably short step size,  you
   can keep that from happening too often.

   The down side to all of this is complexity and (often) taking some
   unnecessary steps for main-belt objects,  where a larger step size
   would work just fine.  I _do_ have a better scheme in mind,  and it's
   implemented in my Find_Orb software... but not here (yet).   */

int integrate_orbit( ELEMENTS *elem, double jd, const double jd_end,
                              const double max_err, const double stepsize)
{
   double delta[6],  posnvel[6];
   int i, j, n_steps = 0;

   for( i = 0; i < 6; i++)
      delta[i] = 0.;
   if( !n_xyzs)
      {
      ephem_stepsize = ( jd_end > jd ? stepsize : -stepsize);
      xyzs = (double *)calloc( 3 * ((size_t)( fabs( jd - jd_end) / stepsize) + 2),
                                 sizeof( double));
      }
   while( jd != jd_end)
      {
      double new_delta[6], jd2;

      jd2 = floor( (jd - 0.5) / stepsize + 0.5) * stepsize + 0.5;
      if( jd < jd_end)    /* integrating forward */
         {
         jd2 += stepsize;
         if( jd2 > jd_end)       /* going past the end;  truncate step */
            jd2 = jd_end;
         }
      else                /* integrating backward */
         {
         jd2 -= stepsize;
         if( jd2 < jd_end)
            jd2 = jd_end;
         }
      assert( jd != jd2);
      assert( fabs( jd - jd2) < stepsize * 2.);
      full_rk_step( elem, delta, new_delta, jd, jd2, max_err);
      memcpy( delta, new_delta, 6 * sizeof( double));
      jd = jd2;
      comet_posn_and_vel( elem, jd, posnvel, posnvel + 3);
      for( j = 0; j < 6; j++)
         {
         posnvel[j] += delta[j];
         delta[j] = 0.;
         }
      elem->epoch = jd;
      elem->gm = SOLAR_GM;
      calc_classical_elements( elem, posnvel, jd, 1);
      if( !ephem_start)
         {
         ephem_start = jd;    /* we start on an even 'stepsize' boundary */
         printf( "ephem start %f\n", ephem_start);
         }
      if( n_xyzs >= 0)
         memcpy( xyzs + n_steps * 3, posnvel, 3 * sizeof( double));
      n_steps++;
      }
   if( !n_xyzs)
      n_xyzs = n_steps;
   comet_posn_and_vel( elem, jd_end, posnvel, posnvel + 3);
   for( i = 0; i < 6; i++)
      posnvel[i] += delta[i];
   elem->epoch = jd_end;
   elem->gm = SOLAR_GM;
   calc_classical_elements( elem, posnvel, jd_end, 1);
   return( 0);
}

int load_vsop_data( void)
{
   FILE *ifile = err_fopen( "vsop.bin", "rb");
   const unsigned vsop_size = 60874u;

   vsop_data = NULL;
   if( ifile)
      {
      vsop_data = (char *)calloc( vsop_size, 1);
      if( vsop_data)
         {
         const size_t bytes_read = fread( vsop_data, 1, vsop_size, ifile);

         assert( bytes_read == vsop_size);
         }
      fclose( ifile);
      }
   return( ifile && vsop_data ? 0 : -1);
}

static double centralize( double ang)
{
   while( ang < 0.)
      ang += PI + PI;
   while( ang > PI + PI)
      ang -= PI + PI;
   return( ang);
}

static int put_elem_into_sof( const char *header, char *buff, const ELEMENTS *elem)
{
   const char *tptr;
   char tbuff[40];

   tptr = strstr( header, "|Tp ");
   assert( tptr);
   full_ctime( tbuff, elem->perih_time, FULL_CTIME_YMD | FULL_CTIME_NO_SPACES
               | FULL_CTIME_LEADING_ZEROES | FULL_CTIME_MONTHS_AS_DIGITS
               | FULL_CTIME_FORMAT_DAY | FULL_CTIME_7_PLACES);
   memcpy( buff + (tptr - header) + 1, tbuff, strlen( tbuff));

   tptr = strstr( header, "|Te ");
   assert( tptr);
   full_ctime( tbuff, elem->epoch,  FULL_CTIME_YMD | FULL_CTIME_NO_SPACES
               | FULL_CTIME_LEADING_ZEROES | FULL_CTIME_MONTHS_AS_DIGITS
               | FULL_CTIME_FORMAT_DAY);
   memcpy( buff + (tptr - header) + 1, tbuff, strlen( tbuff));

   tptr = strstr( header, "|q ");
   assert( tptr);
   snprintf_err( tbuff, 13, " %11.8f", elem->q);
   memcpy( buff + (tptr - header), tbuff, 12);

   tptr = strstr( header, "|i ");
   assert( tptr);
   snprintf_err( tbuff, 12, " %10.6f", elem->incl * 180. / PI);
   memcpy( buff + (tptr - header), tbuff, 11);

   tptr = strstr( header, "|Om ");
   assert( tptr);
   snprintf_err( tbuff, 12, " %10.6f", centralize( elem->asc_node) * 180. / PI);
   memcpy( buff + (tptr - header), tbuff, 11);

   tptr = strstr( header, "|om ");
   assert( tptr);
   snprintf_err( tbuff, 12, " %10.6f", centralize( elem->arg_per) * 180. / PI);
   memcpy( buff + (tptr - header), tbuff, 11);

   tptr = strstr( header, "|e ");
   assert( tptr);
   snprintf_err( tbuff, 12, " %10.8f", elem->ecc);
   memcpy( buff + (tptr - header), tbuff, 11);

   return( 0);
}

static int integrate_unperturbed = 0;

static double try_to_integrate( const char *header, char *buff, const double dest_jd,
                         const double max_err, const double stepsize)
{
   ELEMENTS elem;
   int got_it = 0, pluto_removed = 0;

   if( !integrate_unperturbed)
      {
      const char *tptr = strstr( header, "|Perts");

      if( tptr && !memcmp( buff + (tptr - header), "      ", 6))
         return( 0.);
      }
   got_it = !extract_sof_data_ex( &elem, buff, header, NULL);
   assert( got_it);
   if( !memcmp( buff, "      134340 ", 13))   /* don't let (134340) Pluto */
      if( perturber_mask & PERTURBERS_PLUTO)    /* perturb itself! */
         {
         pluto_removed = 1;
         perturber_mask ^= PERTURBERS_PLUTO;
         }

   if( got_it && dest_jd != 0. && elem.epoch != 0.)
      {
      if( !position_cache)       /* gotta initialize it: */
         position_cache = make_position_cache( elem.epoch, dest_jd, stepsize);
      integrate_orbit( &elem, elem.epoch, dest_jd, max_err, stepsize);
      put_elem_into_sof( header, buff, &elem);
      }

   if( pluto_removed)
      perturber_mask ^= PERTURBERS_PLUTO;
   return( elem.epoch);
}

static void get_sof_element( char *obuff, const size_t obuff_size,
            const char *buff, const char *header, const char *tag)
{
   const char *tptr = strstr( header, tag);

   if( tptr && (tptr == header || tptr[-1] == '|'))
      {
      size_t i;

      buff += tptr - header;
      for( i = 0; i < obuff_size - 1 && tptr[i] != '|' && tptr[i] != '^'; i++)
         obuff[i] = buff[i];
      obuff[i] = '\0';
      }
   else
      obuff[0] = '\0';
}

/* If we're updating a previous result,  we check to see if the designation,
H, G,  reference,  number of observations,  etc.  have changed.  If they
have, the data underlying the orbit have presumably changed,  and we need to
re-integrate that object's orbit.  But if our previous result does contain an
object with the same name and other details,  we don't have to do all the
math to integrate it all over again just to get the same result as before. */

static long compute_hash( const char *buff, const char *header)
{
   long rval = 0;
   const long big_prime = 2141592701L;
   size_t i, j;
   char obuff[30];
   const char *tags[] =
             { "H .", "G .", "Tlast", "Tfirst", "n_o", "rms", "Name" };

   for( i = 0; i < sizeof( tags) / sizeof( tags[0]); i++)
      {
      get_sof_element( obuff, sizeof( obuff), buff, header, tags[i]);
      for( j = 0; obuff[j]; j++)
         rval = rval * big_prime + (long)obuff[j];
      }
   return( rval);
}

static unsigned find_in_table( const long *hashes, const long hash_val)
{
   unsigned i, loc = (unsigned)hash_val % HASH_TABLE_SIZE;

   for( i = 1; hashes[loc] && hashes[loc] != hash_val; i += 2)
      loc = (loc + i) % HASH_TABLE_SIZE;
   return( loc);
}

static void error_exit( void)
{
#ifdef _MSC_VER
   printf( "Hit any key:\n");
   getch( );
#endif
}

#define JAN_1970 2440587.5
#define LINE_SIZE 300

int main( int argc, const char **argv)
{
   FILE *ifile, *ofile, *update_file = NULL, *ephem_file = NULL;
   const char *temp_file_name = "ickywax.ugh";
   const char *output_filename = argv[2];
   long *hashes, *file_offsets, hash_val;
   const char *ephem_filename = NULL;
   double dest_jd, max_err = 1.e-12, stepsize = 2., t_last_printout = 0.;
   double starting_jd = 0.;
   char buff[LINE_SIZE], time_buff[60], header[LINE_SIZE];
   int i, n_integrated = 0, total_asteroids_in_file;
   int max_asteroids = (1 << 30);
#ifdef FORKING
   int n_processes = 0, process_number = 0, child_status;
   bool forking_has_happened = false;
#endif
   int quit = 0, n_found_from_update = 0;
   clock_t t0;
   bool update_existing_file = true;

   if( argc < 4)
      {
      printf( "INTEGRAT takes as command-line arguments the name of an input\n");
      printf( "file in .sof format;  the name of the output\n");
      printf( "file that is to be created;  and the epoch (JD or YYYYMMDD)\n");
      printf( "of that file.  For example: either\n\n");
      printf( "integrat mpcorb.sof 2452600.sof 2452600.5\n\n");
      printf( "integrat mpcorb.sof 2452600.sof 20021122\n\n");
      printf( "would read in the 'mpcorb.sof' file,  and create a new file\n");
      printf( "updated to the epoch JD 2452600.5 = 22 Nov 2002.\n");
      error_exit( );
      return( -1);
      }
   for( i = 1; i < argc; i++)
      if( !strcmp( argv[i], "-u"))
         update_existing_file = false;
   setvbuf( stdout, NULL, _IONBF, 0);
   ifile = err_fopen( argv[1], "rb");
   if( !fgets( header, sizeof( header), ifile))
      {
      fprintf( stderr, "Couldn't read header from original file\n");
      error_exit( );
      return( -1);
      }
   if( update_existing_file && !rename( output_filename, temp_file_name))
      {
      int n_hashes = 0;

      printf( "Using an update\n");
      update_file = err_fopen( temp_file_name, "rb");
      hashes = (long *)calloc( HASH_TABLE_SIZE * 2, sizeof( long));
      file_offsets = hashes + HASH_TABLE_SIZE;
      while( fgets( buff, sizeof( buff), update_file))
         if( (hash_val = compute_hash( header, buff)) != 0L)
            {
            const unsigned hash_loc = find_in_table( hashes, hash_val);

            hashes[hash_loc] = hash_val;
            file_offsets[hash_loc] = ftell( update_file) - strlen( buff);
            n_hashes++;
                      /* We shouldn't try to fill the table more than 80%. */
                      /* If that happens,  raise HASH_TABLE_SIZE.          */
            assert( n_hashes < HASH_TABLE_SIZE * 4 / 5);
            }
      printf( "Got %d hashes\n", n_hashes);
      }
   else
      hashes = file_offsets = NULL;

   dest_jd = get_time_from_string( 0., argv[3], FULL_CTIME_YMD, NULL);
   full_ctime( time_buff, dest_jd, 0);
   snprintf( buff, sizeof( buff),
                   "Integrat version %s %s\nIntegrating to %s = JD %.5f\n",
                    __DATE__, __TIME__, time_buff, dest_jd);
   printf( "%s", buff);
   ofile = err_fopen( output_filename, "wb");
   setvbuf( ofile, NULL, _IONBF, 0);
   if( dest_jd != floor( dest_jd) + .5)
      {
      printf( "WARNING: the MPCORB format can only handle 'standard' 0h TD epochs.\n");
      printf( "Integrat will create elements that give the correct position and velocity\n");
      printf( "at the epoch you've requested;  but the epoch stored in MPCORB format\n");
      printf( "will be rounded to the nearest day (and the mean anomaly suitably\n");
      printf( "corrected.)\n\nHit any key:\n");
#ifdef _MSC_VER
      getch( );
#endif
      }
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'c':
               n_xyzs = 0;
               printf( "Creating 'ephem.dat' file\n");
               break;
            case 'f':
               ephem_filename = argv[i] + 2;
               if( !*ephem_filename && i < argc - 1)
                  ephem_filename = argv[i + 1];
               break;
            case 'n':
               max_asteroids = atoi( argv[i] + 2);
               printf( "Only integrating up to %d objects\n", max_asteroids);
               break;
            case 'p':
               integrate_unperturbed = 1;
               printf( "Integrating unperturbed objects,  too\n");
               break;
            case 'r':
               resync_freq = atoi( argv[i] + 2);
               break;
            case 's':
               stepsize = atof( argv[i] + 2);
               printf( "Step size set at %.2f days\n", stepsize);
               break;
            case 't':
               max_err = atof( argv[i] + 2);
               break;
            case 'u':      /* handled above */
               break;
            case 'v':
               verbose = 1 + atoi( argv[i] + 2);
               printf( "Setting verbose output\n");
               break;
#ifdef FORKING
            case 'z':
               n_processes = atoi( argv[i] + 2);
               break;
#endif
            default:
               printf( "Command-line option '%s' ignored\n", argv[i]);
               break;
            }
   if( ephem_filename)
      {
      jpl_ephemeris = jpl_init_ephemeris( ephem_filename, NULL, NULL);
      if( !jpl_ephemeris)
         {
         printf( "JPL ephemeris file '%s' not found\n", ephem_filename);
         error_exit( );
         return( -3);
         }
      perturber_mask |= PERTURBERS_PLUTO_AND_MOON;
      if( verbose)
         printf( "Using JPL ephemeris file '%s'\n", ephem_filename);
      }

   if( !jpl_ephemeris && load_vsop_data( ))
      {
      printf( "VSOP.BIN not loaded!\n");
      error_exit( );
      return( -4);
      }

   if( !ephem_filename)          /* gotta lump the Moon in with the earth: */
      relative_mass[3] += relative_mass[10];

   /* first,  go through the file to figure out how many asteroids  */
   /* we'll have integrate: */

   total_asteroids_in_file = 0;
   while( fgets( buff, sizeof( buff), ifile)
                     && total_asteroids_in_file < max_asteroids)
      {
      const double tval = try_to_integrate( header, buff, 0., max_err, stepsize);

      if( tval != 0. && starting_jd == 0.)
         {
         starting_jd = tval;
         full_ctime( time_buff, starting_jd, FULL_CTIME_DATE_ONLY | 0x30);
         snprintf_err( buff, sizeof( buff),
                      "'%s' has elements for %s = JD %.1f (and possibly other epochs)\n",
                      argv[1], time_buff, starting_jd);
         printf( "%s", buff);
         }
      if( tval != 0.)
         total_asteroids_in_file++;
      }

   snprintf_err( buff, sizeof( buff),
                 "%d asteroids to be integrated\n", total_asteroids_in_file);
   printf( "%s", buff);
   fputs( header, ofile);

   fseek( ifile, strlen( header), SEEK_SET);

   t0 = clock( );
   while( !quit && fgets( buff, sizeof( buff), ifile)
                                 && n_integrated < max_asteroids)
      {
      bool got_it_from_update = false;

      asteroid_perturber_number = -1;
      if( n_integrated < 4 && !memcmp( buff, "           ", 11))
         switch( atoi( buff))
            {
            case 1:              /* Ceres */
               asteroid_perturber_number = 10;
               break;
            case 2:              /* Pallas */
               asteroid_perturber_number = 11;
               break;
            case 4:              /* Vesta */
               asteroid_perturber_number = 12;
               break;
            default:
               break;
            }
#ifdef FORKING
      if( (perturber_mask & 0x1c00) == 0x1c00 && n_processes
               && !forking_has_happened)
         {
         char outfile_name[50];
         int j;
         const long offset = ftell( ifile);

         forking_has_happened = true;
         fclose( ofile);
         fclose( ifile);
         if( jpl_ephemeris)
            jpl_close_ephemeris( jpl_ephemeris);
         if( update_file)
            fclose( update_file);
         while( process_number < n_processes - 1)
            {
            const pid_t childpid = fork( );

            if( childpid == -1)      /* fork( ) returns -1 on failure */
               {
               perror( "fork"); /* display error message */
               exit(0);
               }
            else if( childpid == 0)     /* we're a child process */
               {
//             printf( "Hi!  I'm child %d.  My PID is %d; parent's is %d\n",
//                      process_number, getpid( ), getppid( ));
               }
            else
               break;       /* break out of loop,  signalling we're a parent */
            process_number++;
            }
         printf( "Hi!  I've got process number %d,  PID %d,  parent's is %d\n",
                        process_number, getpid( ), getppid( ));
         chunk_filename( outfile_name, process_number, 0);
         ofile = err_fopen( outfile_name, "wb");
         if( n_xyzs > 0)
            {
            chunk_filename( outfile_name, process_number, 1);
            ephem_file = err_fopen( outfile_name, "wb");
            }
         ifile = err_fopen( argv[1], "rb");
         fseek( ifile, offset, SEEK_SET);
         if( jpl_ephemeris)
            jpl_ephemeris = jpl_init_ephemeris( ephem_filename, NULL, NULL);
         if( update_file)
            update_file = err_fopen( temp_file_name, "rb");
         j = 0;
         while( j < process_number && fgets( buff, sizeof( buff), ifile))
            j++;
         }
#endif
      if( update_file && asteroid_perturber_number == -1
                          && (hash_val = compute_hash( header, buff)) != 0L)
         {
         char buff2[220];
         const unsigned hash_loc = find_in_table( hashes, hash_val);

         if( hashes[hash_loc])
            {
            assert( hashes[hash_loc] == hash_val);
            fseek( update_file, file_offsets[hash_loc], SEEK_SET);
            if( fgets( buff2, sizeof( buff2), update_file)
                         && !memcmp( buff2, buff, 20)
                         && !memcmp( buff2 + 105, buff + 105, 97))
               {
               strcpy( buff, buff2);
               got_it_from_update = true;
               n_found_from_update++;
               }
            }
         }

      if( !got_it_from_update &&
                  try_to_integrate( header, buff, dest_jd, max_err, stepsize) != 0.)
         {
         clock_t t = clock( );
         const double elapsed_time = (double)(t - t0) / (double)CLOCKS_PER_SEC;

         if( asteroid_perturber_number > 0)
            {
            const char *pert_text[3] = { "(1) Ceres", "(2) Pallas", "(4) Vesta" };

            assert( asteroid_perturber_number >= 10);
            assert( asteroid_perturber_number < 13);
            printf( "Perturber %s calculated\n", pert_text[asteroid_perturber_number - 10]);
            perturber_mask |= (1L << asteroid_perturber_number);
            }
         n_integrated++;
         if( verbose > 1)
            {
            char tbuff[30];

            memcpy( tbuff, buff, 29);
            tbuff[29] = '\0';
            printf( "%s: %.2f seconds;  %5d steps: %5d integrated\n",
                            tbuff, elapsed_time, n_steps_taken, n_integrated);
            t0 = t;        /* restart the clock */
            n_steps_taken = 0;
            }
         else if( elapsed_time > t_last_printout + 1.)
            {
            int divisor = 1;
            double t_remains;

            t_last_printout = elapsed_time;
#ifdef FORKING
            if( forking_has_happened)
               divisor = n_processes;
#endif
            t_remains = (double)(total_asteroids_in_file / divisor - n_integrated)
                           * elapsed_time / (double)n_integrated;
            printf( "%.0f seconds elapsed;  %.0f seconds remain; %d done %d    \r",
                        elapsed_time, t_remains,
                       n_integrated * divisor,
                       n_found_from_update);
            }
#ifdef _MSC_VER
         if( kbhit( ))
            if( getch( ) == 27)
               quit = 1;
#endif
         }
      fputs( buff, ofile);
      if( ephem_file)
         {
         const size_t n_written = fwrite( xyzs, 3 * sizeof( double), n_xyzs, ephem_file);

         assert( n_written == (size_t)n_xyzs);
         }
#ifdef FORKING
      if( forking_has_happened)
         {
         int j = 1;

         while( j < n_processes && fgets( buff, sizeof( buff), ifile))
            j++;
         }
#endif
      }
   if( jpl_ephemeris)
      jpl_close_ephemeris( jpl_ephemeris);
   fclose( ifile);
   fclose( ofile);
   if( ephem_file)
      fclose( ephem_file);
#ifdef FORKING
   if( forking_has_happened)
      {
      printf( "Process %d is done\n", process_number);
      wait( &child_status); /* wait for child to exit, and store its status */
      printf( "Waiting is over for process %d\n", process_number);
      if( !process_number)
         {
         FILE **ifiles = (FILE **)calloc( n_processes, sizeof( FILE *));

         for( i = 0; i < n_processes; i++)
            ifiles[i] = err_fopen( chunk_filename( buff, i, 0), "rb");
         ofile = err_fopen( output_filename, "ab");
         i = 0;
         while( fgets( buff, sizeof( buff), ifiles[i]))
            {
            fputs( buff, ofile);
            i = (i + 1) % n_processes;
            }
         for( i = 0; i < n_processes; i++)
            {
            fclose( ifiles[i]);
            unlink( chunk_filename( buff, i, 0));
            }
         fclose( ofile);

         if( ephem_file)
            {
            for( i = 0; i < n_processes; i++)
               ifiles[i] = err_fopen( chunk_filename( buff, i, 1), "rb");
            ofile = err_fopen( "ephem.dat", "wb");
            fprintf( ofile, "%f %f %d\n", ephem_start, ephem_stepsize, n_xyzs);
            i = 0;
            while( fread( xyzs, 3 * sizeof( double), n_xyzs, ifiles[i]) == (size_t)n_xyzs)
               {
               const size_t n_written = fwrite( xyzs, 3 * sizeof( double), n_xyzs, ofile);

               assert( n_written == (size_t)n_xyzs);
               i = (i + 1) % n_processes;
               }
            for( i = 0; i < n_processes; i++)
               {
               fclose( ifiles[i]);
               unlink( chunk_filename( buff, i, 1));
               }
            fclose( ofile);
            }
         free( ifiles);
         }
      }
#endif
   assert( position_cache);
   free( position_cache);
   return( 0);
}
