/* comets.h: header file for comet/asteroid functions
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

#define ELEMENTS struct elements

#pragma pack(4)
ELEMENTS
   {
   double perih_time, q, ecc, incl, arg_per, asc_node;
   double epoch,  mean_anomaly;
            /* derived quantities: */
   double lon_per, minor_to_major;
   double perih_vec[3], sideways[3];
   double angular_momentum, major_axis, t0, w0;
   double abs_mag, slope_param, gm;
   int is_asteroid, central_obj;
   };
#pragma pack( )

/* Note that in the above structure,  t0 = 1/n = time to move */
/* one radian in mean anomaly = orbital period / (2*pi).      */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

// void calc_vectors( ELEMENTS *elem, const double sqrt_gm);
int DLL_FUNC calc_classical_elements( ELEMENTS *elem, const double *r,
                             const double t, const int ref);
int DLL_FUNC comet_posn_and_vel( ELEMENTS DLLPTR *elem, double t,
                  double DLLPTR *loc, double DLLPTR *vel);
int DLL_FUNC comet_posn( ELEMENTS DLLPTR *elem, double t, double DLLPTR *loc);       /* astfuncs.c */
void DLL_FUNC derive_quantities( ELEMENTS DLLPTR *e, const double gm);
int DLL_FUNC setup_elems_from_ast_file( ELEMENTS DLLPTR *class_elem,
              const uint32_t DLLPTR *elem, const double t_epoch);
double DLL_FUNC phase_angle_correction_to_magnitude( const double phase_angle,
                                 const double slope_param);
int setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                             const double t_cen);
int extract_sof_data( ELEMENTS *elem, const char *buff, const char *header);
int extract_sof_data_ex( ELEMENTS *elem, const char *buff, const char *header,
                        double *extra_info);                /* sof.cpp */
double extract_yyyymmdd_to_jd( const char *buff);           /* sof.cpp */

typedef struct
{
   double obj1_true_anom, jd1;       /* these are set in find_moid_full */
   double obj2_true_anom, jd2;       /* NOT SET YET */
   double barbee_speed;              /* in AU/day */
} moid_data_t;

double DLL_FUNC find_moid_full( const ELEMENTS *elem1, const ELEMENTS *elem2, moid_data_t *mdata);

#ifdef __cplusplus
}
#endif
