/* date.h: header file for color conversions
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

double v_minus_i_to_b_minus_v( const double v_minus_i);         /* colors.c */
double b_minus_v_to_v_minus_i( const double b_minus_v);         /* colors.c */
double v_minus_i_to_v_minus_r( const double v_minus_i);         /* colors.c */
double v_minus_r_to_v_minus_i( const double v_minus_r);         /* colors.c */
double v_minus_r_to_b_minus_v( const double v_minus_r);         /* colors.c */
double b_minus_v_to_v_minus_r( const double b_minus_v);         /* colors.c */
double johnson_b_minus_v_from_tycho_b_minus_v( const double b_v_t);
double johnson_v_from_tycho_b_minus_v( const double b_v_t, const double tycho_v);
int tycho_to_johnson_colors( double bt_minus_vt, double *results);
