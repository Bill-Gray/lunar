/* getplane.cpp: functions for computing planetary ephems

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
#include "watdefs.h"
#include "lunar.h"
#include "afuncs.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

/* The compute_planet( ) function is supposed to provide a "standardized"
   way to get positions for planets and the moon,  in any of four systems.
   You give it a pointer to the VSOP data,  a planet number (0=sun,
   1=mercury, ... 9 = Pluto,  10 = the Moon),  and a time in centuries
   from J2000.  The returned array of fifteen doubles gives you the object
   position in five different systems:

      ovals[0, 1, 2] = lon, lat, r,  heliocentric ecliptic of date;
      ovals[3, 4, 5] = x, y, z,  heliocentric ecliptic of date;
      ovals[6, 7, 8] = x, y, z,  heliocentric equatorial of date;
      ovals[9, 10, 11] = x, y, z,  heliocentric equatorial, J2000.0
      ovals[12, 13, 14] = x, y, z,  heliocentric ecliptical, J2000.0

   I intended to rig this up so you could keep on going to 11=Io, 12=Europa,
etc.  This would make all sorts of sense.  But I've not done it yet. */

int DLL_FUNC compute_planet( const char FAR *vsop_data, const int planet_no,
                         const double t_c, double DLLPTR *ovals)
{
   double lat, lon, r;
   const double obliquit = mean_obliquity( t_c);
   double matrix[9];
   const double obliq_2000 = 23.4392911 * PI / 180.;

            /* first,  compute polar heliocentric in eclip of date */
   if( planet_no != 10)
      {
      lon = calc_vsop_loc( vsop_data, planet_no, 0, t_c, 0.);
      lat = calc_vsop_loc( vsop_data, planet_no, 1, t_c, 0.);
      r   = calc_vsop_loc( vsop_data, planet_no, 2, t_c, 0.);
      }
   else
      {
      double fund[N_FUND];

      lunar_fundamentals( vsop_data, t_c, fund);
      lat = lunar_lat( vsop_data, fund, 0L);
      lunar_lon_and_dist( vsop_data, fund, &lon, &r, 0L);
      lon *= PI / 180.;
      lat *= PI / 180.;
      r /= AU_IN_KM;      /* from km to AU */
      }
   ovals[0] = lon;
   ovals[1] = lat;
   ovals[2] = r;
            /* next, compute polar cartesian in eclip of date */
   ovals[3] = cos( lon) * cos( lat) * r;
   ovals[4] = sin( lon) * cos( lat) * r;
   ovals[5] =             sin( lat) * r;
            /* next, compute polar cartesian in eclip of date, */
            /* but in equatorial coords */
   FMEMCPY( ovals + 6, ovals + 3, 3 * sizeof( double));
   rotate_vector( ovals + 6, obliquit, 0);
            /* next, precess to get J2000.0 equatorial values */
   setup_precession( matrix, 2000. + t_c * 100., 2000.);
   precess_vector( matrix, ovals + 6, ovals + 9);
            /* Finally,  rotate equatorial J2000.0 into ecliptical J2000 */
   FMEMCPY( ovals + 12, ovals + 9, 3 * sizeof( double));
   rotate_vector( ovals + 12, -obliq_2000, 0);
   return( 0);
}
