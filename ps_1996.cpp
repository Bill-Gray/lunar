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
#include <string.h>
#include "watdefs.h"
#include "mpc_func.h"
#include "lunar.h"
#include "afuncs.h"
#include "date.h"

/* Examples for computing topocentric positions of the moon
and planets,  using PS1996 and ELP-82.  For JD 2456789.0,
one should get the following output:

Example for ELP82 and PS1996
Lunar distance from geocenter:    389666.44589 km
Topocentric Lunar position: 12h 31m 29.5284s   -05 56' 05.657"
EMB position: 15h 12m 21.0066s   -17 52' 47.858"
Earth position: 15h 12m 21.7098s   -17 52' 43.784"
Mercury    04h 19m 45.1664s   +23 29' 30.443"
Venus      00h 37m 22.6957s   +02 10' 06.901"
Mars       12h 35m 57.4619s   -02 38' 06.749"
Jupiter    07h 11m 27.6660s   +22 44' 07.445"
Saturn     15h 12m 00.3509s   -15 16' 54.541"
Uranus     00h 54m 01.5742s   +05 04' 33.671"
Neptune    22h 36m 32.9155s   -09 31' 54.182"
Pluto      18h 56m 15.2128s   -20 08' 01.579"

   In computing planetary/lunar/asteroid positions,  there are
several odd things that need to be done.  First,  you have to
account for the fact that three different time scales are in
use.  "UT" (Universal Time) corresponds to the "usual" concept
of time as being linked to the earth's rotation;  noon UT is
expected to match up with when the sun (on average) is highest
in the sky as seen from Greenwich,  England.  "TD" (Temps
Dynamic,  or Dynamical Time) corresponds to a physicist's idea
of time:  time as measured by an atomic clock.  Since the
earth's rotation isn't perfectly uniform,  there is a varying
difference between TD and UT,  known as "Delta-T".  (It's
currently about 70 seconds,  and increasing.)

   "UTC" (Coordinated Universal Time) is something of an hybrid.
The lengths of seconds are determined by atomic clocks.  When
the clocks don't quite line up with the earth's rotation (which
is inevitable because the clocks represent a "uniform" flow of
time,  but the earth's rotation is irregular),  a "leap second"
is inserted.  Strange as it may seem,  the last day of December
or June is sometimes one second longer than other days,  in order
to ensure that UTC stays within .9 seconds of UT (i.e.,  within
.9 seconds of the earth's rotation).  Because the rotation really
is irregular,  you can't tell very far in advance when a leap
second will be inserted;  you usually get about six months of
warning.

   Anyway.  The upshot of all this is that we need to handle
UTC,  because that's what the user will usually "think" in;
it's assumed that a UTC time will be given to this program.
We also need to know UT,  because that corresponds to the earth's
rotation;  and we need to know how much the earth has rotated
in order to compute where the observer,  on the surface of the
earth,  is positioned relative to the geocenter.  And finally,
we need to know TD,  because that is the "uniformly flowing"
time scale best suited to planetary and lunar ephemerides.

   (Note that a purist would turn pale at the brief history
of time systems given above.  I've omitted lots of nuances for
the sake of clarity.)

   Anyway.  The next odd thing we need to consider is that PS1996
(the theory used to compute planetary positions) gives the
coordinates of the Earth/Moon barycenter (EMB),  rather than those
of the Earth.  (JPL ephemerides do the same thing.  VSOP comes with
"Earth center" series and "EMB" series.)  Using EMB is actually a
pretty reasonable thing to do.  The EMB moves along in a fairly
smooth and continuous manner.  The earth's center wobbles around
the EMB once a month,  resulting in more complicated motion and a
longer series of coefficients.

   To get around this,  the code computes both the EMB position
and the earth-moon vector.  (The latter isn't part of PS1996;
it's computed using the separate ELP-82 theory.)  That lets us
"adjust" the EMB position,  after which we've a good handle on
where the center of the earth is.

   However,  one is often more interested in where an observer, on the
center of the earth,  happens to be. 'compute_topocentric_offset()'
figures out the vector between the earth's center and an observer at a
specified lat/lon and altitude above the geoid.

   Next,  we compute the position of the target planet.  Then
we hit a small snag:  the speed of light is finite,  so what
we really need is to know where it was back when the light left
the planet and headed our way.  So if we find that Jupiter,  for
example,  was 4.5 AU away,  we divide the distance by the speed
of light and then say,  "OK... where was Jupiter that much
earlier (about 37 minutes)?"  After a couple of iterations
(this code does three,  which is probably one more than is
really needed),  we have a good position for the target planet.

   Subtract the (heliocentric) Earth vector from the (heliocentric)
planet vector,  and you've a good Earth-to-planet vector. You can
then take that vector,  convert to polar coordinates,  and show an
RA/dec (as is done in the 'make_ra_dec_string' function.)

   At this point,  you may be done.  You have a topocentric position
in the J2000,  or ICRF,  frame of reference. (Purists will point out
that J2000 isn't exactly the same as ICRF, but it's close enough for
current purposes.)  This is the position in a system that is lined
up with distant quasars and corresponds to a physicist's concept of
an inertial frame of reference:  one whose axes are fixed.  If you
looked up a position on a star chart,  this is the coordinate system
you would use.

   However,  you may need "apparent coordinates",  in which the
axes are defined by the earth's rotation.  Since the earth wobbles
due to precession and nutation,  the axes used for apparent coordinates
wobble as well;  and one needs to account for the effects of the
aberration of light as well.  All that is outside the scope of this
little demo program (though routines for precession,  nutation,  and
aberration are part of this library).   */

#define J2000 2451545.0
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define EARTH_MOON_BALANCE 81.300569404492336
#define J2000_OBLIQUITY  (23.4392911 * PI / 180)

static double make_ra_dec_string( double *vect, char *buff)
{
   double ra = atan2( vect[1], vect[0]) * 12. / PI;
   const double dist =
              sqrt( vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
   double dec = asin( vect[2] / dist) * 180. / PI;
   long units;

   if( ra < 0.)
      ra += 24.;

   units = (long)( ra * 3600. * 10000.);
   sprintf( buff, "%02ldh %02ldm %02ld.%04lds   ",
         units / 36000000L, (units / 600000L) % 60L,
         (units / 10000L) % 60L, units % 10000L);
   buff += strlen( buff);

   if( dec < 0.)
      {
      *buff++ = '-';
      dec = -dec;
      }
   else
      *buff++ = '+';

   units = (long)( dec * 3600. * 1000.);
   sprintf( buff, "%02ld %02ld' %02ld.%03ld\"   ",
         units / 3600000L, (units / 60000L) % 60L,
         (units / 1000L) % 60L, units % 1000L);
   return( dist);
}

#define EARTH_MAJOR_AXIS    6378137.
#define EARTH_MINOR_AXIS    6356752.

static void compute_topocentric_offset( double *earth_vect, const double lat,
               const double lon, const double alt_in_meters, const double jd_ut)
{
   double precess_matrix[9];
   double rho_cos_phi, rho_sin_phi;
   int i;

   lat_alt_to_parallax( lat, alt_in_meters, &rho_cos_phi, &rho_sin_phi,
            EARTH_MAJOR_AXIS, EARTH_MINOR_AXIS);
   rho_cos_phi *= EARTH_MAJOR_AXIS / AU_IN_METERS;
   rho_sin_phi *= EARTH_MAJOR_AXIS / AU_IN_METERS;
   calc_planet_orientation( 3, 0, jd_ut, precess_matrix);
   spin_matrix( precess_matrix, precess_matrix + 3, lon);
   for( i = 0; i < 3; i++)
      earth_vect[i] = (rho_cos_phi * precess_matrix[i]
                     + rho_sin_phi * precess_matrix[i + 6]);
}

#define N_PASSES 3

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( "elp82.dat", "rb");
   double xyz[4], xyz_topo[3], topo_offset[3];
   double jd_utc = get_time_from_string( 0., argv[1],
                  CALENDAR_JULIAN_GREGORIAN, NULL);
   double jd_td = jd_utc + td_minus_utc( jd_utc) / seconds_per_day;
   double jd_ut = jd_td - td_minus_ut( jd_utc) / seconds_per_day;
   double earth_vect[6];
   double lon = -69.9 * PI / 180.;  /* Project Pluto corporate headquarters is */
   double lat = 44.01 * PI / 180.;  /* roughly at W 69.9, N 44.01, 100 meters alt */
   double alt_in_meters = 100.;
   int i, planet_no;
   char buff[80];
   void *p;

   printf( "Example for ELP82 and PS1996\n");
   if( !ifile)
      {
      printf( "Couldn't find 'elp82.dat'\n");
      return( -1);
      }
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'l':
               sscanf( argv[i] + 2, "%lf,%lf,%lf",
                        &lon, &lat, &alt_in_meters);
               printf( "Setting lon = %f, lat = %f, alt = %f meters\n",
                        lon, lat, alt_in_meters);
               lon *= PI / 180.;
               lat *= PI / 180.;
               break;
            default:
               printf( "Command line option '%s' ignored\n", argv[i]);
            }

   compute_elp_xyz( ifile, (jd_td - J2000) / 36525., 0., xyz);
   printf( "Lunar distance from geocenter: %15.5f km\n", xyz[3]);
                  /* Rotate from ecliptic J2000 to equatorial J2000: */
   rotate_vector( xyz, J2000_OBLIQUITY, 0);
   for( i = 0; i < 3; i++)          /* cvt earth-moon vector from km */
      xyz[i] /= AU_IN_KM;           /* into AU */
   make_ra_dec_string( xyz, buff);
   printf( "Geocentric lunar position: %s\n", buff);
   compute_topocentric_offset( topo_offset, lat, lon, alt_in_meters, jd_ut);
   for( i = 0; i < 3; i++)           /* make a topocentric version of the */
      xyz_topo[i] = xyz[i] - topo_offset[i];         /* earth-moon vector */
   make_ra_dec_string( xyz_topo, buff);
   printf( "Topocentric Lunar position: %s\n", buff);
   printf( "Topocentric lunar distance: %f km\n",
               vector3_length( xyz_topo) * AU_IN_KM);
   fclose( ifile);

               /* We first earth_vect with the position of the         */
               /* Earth-Moon barycenter (EMB).  Yes,  it would be nice */
               /* if the PS1996 gave us the earth's position instead.  */
               /* But it does not.   */

   ifile = fopen( "ps_1996.dat", "rb");
   if( !ifile)
      {
      printf( "Couldn't find 'ps_1996.dat'\n");
      return( -1);
      }
   p = load_ps1996_series( ifile, jd_td, 3);
   if( !p)
      {
      printf( "No planetary data before 1900 or after 2100\n");
      return( -2);
      }
   get_ps1996_position( jd_td, p, earth_vect, 0);
   free( p);

   make_ra_dec_string( earth_vect, buff);
   printf( "EMB position: %s\n", buff);
   for( i = 0; i < 3; i++)
      earth_vect[i] -= xyz[i] / (EARTH_MOON_BALANCE + 1.);
   for( i = 0; i < 3; i++)
      earth_vect[i] += topo_offset[i];
   make_ra_dec_string( earth_vect, buff);
   printf( "Earth position: %s\n", buff);

   for( planet_no = 1; planet_no <= 9; planet_no++)
      if( planet_no != 3)
         {
         static const char *planet_names[10] = {
               "Sun", "Mercury", "Venus", "Earth",
               "Mars", "Jupiter", "Saturn", "Uranus",
               "Neptune", "Pluto" };

         p = load_ps1996_series( ifile, jd_td, planet_no);
         if( !p)
            strcpy( buff, "No data for that time");
         else
            {
            double dist = 0., state_vect[6], delta[3];
            int pass;

            for( pass = 0; pass < N_PASSES; pass++)
               {
               get_ps1996_position( jd_td - dist * AU_IN_KM / (SPEED_OF_LIGHT * seconds_per_day),
                               p, state_vect, 0);
               dist = 0.;
               for( i = 0; i < 3; i++)
                  {
                  delta[i] = state_vect[i] - earth_vect[i];
                  dist += delta[i] * delta[i];
                  }
               dist = sqrt( dist);
               }
            free( p);
            make_ra_dec_string( delta, buff);
            }
         printf( "%-10s %s\n", planet_names[planet_no], buff);
         }

   fclose( ifile);
   return( 0);
}
