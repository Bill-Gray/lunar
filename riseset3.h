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

#define PLANET_DATA struct planet_data

PLANET_DATA
   {
   double ecliptic_loc[3], equatorial_loc[3], altaz_loc[3];
   double r, ecliptic_lon, ecliptic_lat, jd;
   double hour_angle;
   };

int fill_planet_data( PLANET_DATA *pdata, const int planet_no, const double jd,
                  const double observer_lat, const double observer_lon,
                  const char *vsop_data);
double look_for_rise_set( const int planet_no,
                  const double jd0, const double jd1,
                  const double observer_lat, const double observer_lon,
                  const char *vsop_data, int *is_setting);
char *load_file_into_memory( const char *filename, size_t *filesize);
