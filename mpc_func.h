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

#ifndef DLL_FUNC
   #include "watdefs.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

bool is_valid_mpc_code( const char *mpc_code);        /* mpc_fmt.cpp */
double extract_date_from_mpc_report( const char *buff, unsigned *format);
int get_ra_dec_from_mpc_report( const char *ibuff,    /* mpc_fmt.cpp */
                       int *ra_format, double *ra, double *ra_precision,
                       int *dec_format, double *dec, double *dec_precision);

char net_name_to_byte_code( const char *net_name);
const char *byte_code_to_net_name( const char byte_code);
int extract_region_data_for_lat_lon( FILE *ifile, char *buff,
            const double lat_in_degrees, const double lon_in_degrees);

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

int DLL_FUNC text_search_and_replace( char *str, const char *oldstr,
                                     const char *newstr);
int get_mpc_code_info( mpc_code_t *cinfo, const char *buff);
int get_xxx_location_info( mpc_code_t *cinfo, const char *buff);
int get_lat_lon_info( mpc_code_t *cinfo, const char *buff);
double get_ra_from_string( const char *buff, int *bytes_read);
double get_dec_from_string( const char *buff, int *bytes_read);
void output_angle_to_buff( char *obuff, double angle, int precision);
void output_signed_angle_to_buff( char *obuff, const double angle,
                               const int precision);

double point_to_ellipse( const double a, const double b,
                         const double x, const double y, double *dist);
int lat_alt_to_parallax( const double lat, const double ht_in_meters,
            double *rho_cos_phi, double *rho_sin_phi,
            const double major_axis_in_meters,
            const double minor_axis_in_meters);    /* mpc_code.cpp */
int create_mpc_packed_desig( char *packed_desig, const char *obj_name);

void *init_ades2mpc( void);
int xlate_ades2mpc( void *context, char *obuff, const char *buff);
int xlate_ades2mpc_in_place( void *context, char *buff);
int free_ades2mpc_context( void *context);
int fgets_with_ades_xlation( char *buff, const size_t len,
                                      void *ades_context, FILE *ifile);
int mutant_hex_char_to_int( const char c);
char int_to_mutant_hex_char( const int ival);
int get_mutant_hex_value( const char *buff, size_t n_digits);
int encode_value_in_mutant_hex( char *buff, size_t n_digits, int value);
int unpack_mpc_desig( char *obuff, const char *packed);
int unpack_unaligned_mpc_desig( char *obuff, const char *packed);

#define OBJ_DESIG_ASTEROID_PROVISIONAL   0
#define OBJ_DESIG_ASTEROID_NUMBERED      1
#define OBJ_DESIG_COMET_PROVISIONAL      2
#define OBJ_DESIG_COMET_NUMBERED         3
#define OBJ_DESIG_NATSAT_PROVISIONAL     4
#define OBJ_DESIG_NATSAT_NUMBERED        5
#define OBJ_DESIG_ARTSAT                 6
#define OBJ_DESIG_OTHER                 -1

/*   The offset between a satellite observation and the earth or sun
is stored in a second line,  as described at

https://www.minorplanetcenter.net/iau/info/SatelliteObs.html

   These are sometimes -- nay,  frequently -- mangled (decimal points
in odd places,  etc.)  get_satellite_offset() will do its best to
recover from such things,  but may fail.        */

#define SATELL_COORD_ERR_NO_ERROR            0
#define SATELL_COORD_ERR_BAD_SIGN           -1
#define SATELL_COORD_ERR_BAD_NUMBER         -2
#define SATELL_COORD_ERR_NO_DECIMAL         -3
#define SATELL_COORD_ERR_DECIMAL_MISPLACED  -4
#define SATELL_COORD_ERR_UNKNOWN_OFFSET     -5
#define SATELL_COORD_ERR_EXACTLY_ZERO       -6
#define SATELL_COORD_ERR_INSIDE_EARTH       -7

#define N_SATELL_COORD_ERRORS                8

int DLL_FUNC get_satellite_offset( const char *iline, double xyz[3]);

#ifdef __cplusplus
}
#endif

#endif
