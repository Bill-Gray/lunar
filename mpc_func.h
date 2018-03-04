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

#ifndef MPC_FUNC_H_INCLUDED
#define MPC_FUNC_H_INCLUDED

bool is_valid_mpc_code( const char *mpc_code);        /* mpc_fmt.cpp */
double extract_date_from_mpc_report( const char *buff, unsigned *format);
int get_ra_dec_from_mpc_report( const char *ibuff,    /* mpc_fmt.cpp */
                       int *ra_format, double *ra, double *ra_precision,
                       int *dec_format, double *dec, double *dec_precision);

char net_name_to_byte_code( const char *net_name);
const char *byte_code_to_net_name( const char byte_code);

typedef struct
{
   double lat, lon;     /* in radians */
   double alt;          /* in metres */
   double rho_cos_phi, rho_sin_phi;    /* in planet radii */
   int planet, prec1, prec2, format;
   const char *name;
   char code[5];
} mpc_code_t;

#define MPC_CODE_PARALLAXES         1
#define MPC_CODE_LAT_LON_ALT        2
#define MPC_CODE_SATELLITE          3

int get_mpc_code_info( mpc_code_t *cinfo, const char *buff);
double point_to_ellipse( const double a, const double b,
                         const double x, const double y, double *dist);

#endif
