/* lunar.h: header file for basic astrometric functions
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

#ifndef LUNAR_H_INCLUDED
#define LUNAR_H_INCLUDED

#define N_FUND 9

#ifdef __cplusplus
extern "C" {
#endif

#ifndef AU_IN_KM
#define AU_IN_KM 1.495978707e+8
#endif

#ifndef AU_IN_METERS
#define AU_IN_METERS (AU_IN_KM * 1000.)
#endif

#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 299792.458
#endif

#ifndef AU_PER_DAY
#define AU_PER_DAY (86400. * SPEED_OF_LIGHT / AU_IN_KM)
#endif

double DLL_FUNC lunar_lat( const void FAR *data, const double DLLPTR *fund,
                                           const long precision);
int DLL_FUNC lunar_lon_and_dist( const void FAR *data, const double DLLPTR *fund,
                 double DLLPTR *lon, double DLLPTR *r, const long precision);

int DLL_FUNC unload_ps1996_series( void *p);
int DLL_FUNC get_ps1996_position( const double jd, const void *iptr,
                        double *state_vect, const int compute_velocity);
#ifdef SEEK_CUR
void * DLL_FUNC load_ps1996_series( FILE *ifile, double jd, int planet_no);
int DLL_FUNC compute_elp_xyz( FILE *ifile, const double t_cen, const double prec,
                     double *ecliptic_xyz_2000);
int DLL_FUNC calc_big_vsop_loc( FILE *ifile, const int planet,
                      double *ovals, double t, const double prec0);
#endif

int DLL_FUNC lunar_fundamentals( const void FAR *data, const double t,
                                        double DLLPTR *fund);
long double DLL_FUNC find_nearest_lunar_phase_time(
                         const int phase_idx, const long double t2k);
double DLL_FUNC mean_obliquity( const double t_cen);
int DLL_FUNC calc_pluto_loc( const void FAR *data, double DLLPTR *loc,
                                const double t, const long precision);
int DLL_FUNC calc_jsat_loc( const double jd, double DLLPTR *jsats,
                         const int sats_wanted, const long precision);
int DLL_FUNC calc_ssat_loc( const double t, double DLLPTR *ssat,
                                const int sat_wanted, const long precision);
void DLL_FUNC calc_triton_loc( const double jd, double *vect);
double DLL_FUNC calc_vsop_loc( const void FAR *data, const int planet,
                          const int value, double t, double prec);
int DLL_FUNC nutation( const double t, double DLLPTR *d_lon,
                                       double DLLPTR *d_obliq);
int DLL_FUNC compute_planet( const char FAR *vsop_data, const int planet_no,
            const double t_c, double DLLPTR *ovals);
int DLL_FUNC calc_planet_orientation( const int planet_no, const int system_no,
               const double jd, double *matrix);
int DLL_FUNC planet_radii( const int planet_no, double *radii_in_km);
double DLL_FUNC planet_rotation_rate( const int planet_no, const int system_no);
int DLL_FUNC load_cospar_file( const char *filename);
int DLL_FUNC evaluate_rock( const double jd, const int jpl_id,
                                                  double *output_vect);
double planet_radius_in_meters( const int planet_idx);   /* mpc_code.cpp */
double planet_axis_ratio( const int planet_idx);         /* mpc_code.cpp */

#ifdef __cplusplus
}
#endif

#endif /* #ifndef LUNAR_H_INCLUDED */
