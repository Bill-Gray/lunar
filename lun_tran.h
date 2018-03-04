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

double get_lunar_transit_time( const int year, const int month,
   const int day, const double latitude, const double longitude,
   const int time_zone, const int dst, const int real_transit);
void format_hh_mm( char *buff, const double time);
int get_zip_code_data( const int zip_code, double *latitude,
      double *longitude, int *time_zone, int *use_dst,
      char *place_name);
