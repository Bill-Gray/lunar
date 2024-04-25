/* afuncs.h: header files for basic astronomical functions
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
#ifndef AFUNCS_H_INCLUDED
#define AFUNCS_H_INCLUDED

#include <stdint.h>              /* required for int64_t #define */

#ifndef DPT
#define DPT struct dpt

DPT
   {
   double x, y;
   };
#endif /* #ifndef DPT */

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

#ifdef FRENCH_REPUBLICAN_CLOCK
   #define  hours_per_day      10
   #define minutes_per_hour   100
   #define seconds_per_minute 100
#else
   #define hours_per_day       24
   #define minutes_per_hour    60
   #define seconds_per_minute  60
#endif

#define seconds_per_hour   (seconds_per_minute * minutes_per_hour)
#define seconds_per_day    (seconds_per_hour * hours_per_day)
#define minutes_per_day    (minutes_per_hour * hours_per_day)

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

int DLL_FUNC constell_from_ra_dec( const double ra_degrees_1875,
                                   const double dec_degrees_1875,
                                   char DLLPTR *constell_name);
void DLL_FUNC make_var_desig( char DLLPTR *buff, int var_no);
int DLL_FUNC decipher_var_desig( const char DLLPTR *desig);
int DLL_FUNC setup_precession( double DLLPTR *matrix, const double year_from,
                               const double year_to);   /* precess.c */
int DLL_FUNC setup_ecliptic_precession( double DLLPTR *matrix,
                    const double year_from, const double year_to);
int DLL_FUNC setup_precession_with_nutation( double DLLPTR *matrix,
                    const double year);         /* precess.c */
int DLL_FUNC setup_precession_with_nutation_delta( double DLLPTR *matrix,
                    const double year,          /* precess.c */
             const double delta_nutation_lon, const double delta_nutation_obliq);
int DLL_FUNC precess_vector( const double DLLPTR *matrix,
                                      const double DLLPTR *v1,
                                      double DLLPTR *v2);    /* precess.c */
int DLL_FUNC deprecess_vector( const double DLLPTR *matrix,
                                      const double DLLPTR *v1,
                                      double DLLPTR *v2);    /* precess.c */
int DLL_FUNC precess_ra_dec( const double DLLPTR *matrix,
                        double DLLPTR *p_out,
                        const double DLLPTR *p_in, int backward);
void DLL_FUNC rotate_vector( double DLLPTR *v, const double angle,
                                          const int axis);
void DLL_FUNC polar3_to_cartesian( double *vect, const double lon,
                                          const double lat);
double DLL_FUNC vector3_length( const double *vect);
double DLL_FUNC dot_product( const double *a, const double *b);
double DLL_FUNC normalize_vect3( double *vect);
void DLL_FUNC vector_cross_product( double *xprod, const double *a,
                                const double *b);
         /* Following two functions convert 3-D vectors to/from */
         /* ecliptic & equatorial _J2000_ systems:              */
void DLL_FUNC equatorial_to_ecliptic( double *vect);
void DLL_FUNC ecliptic_to_equatorial( double *vect);
void DLL_FUNC set_identity_matrix( double DLLPTR *matrix);
void DLL_FUNC invert_orthonormal_matrix( double DLLPTR *matrix);
void DLL_FUNC spin_matrix( double *v1, double *v2, const double angle);
void DLL_FUNC pre_spin_matrix( double *v1, double *v2, const double angle);
double DLL_FUNC td_minus_ut( const double jd);                /* delta_t.c */
double DLL_FUNC td_minus_utc( const double jd_utc);           /* delta_t.c */
double DLL_FUNC tdb_minus_utc( const double jd_utc);
long double DLL_FUNC tdb_minus_tdt( const long double t_centuries);
void DLL_FUNC reset_td_minus_dt_string( const char *string);  /* delta_t.c */

double DLL_FUNC green_sidereal_time( double jd_ut);          /* alt_az.c */
double DLL_FUNC acose( const double arg);                    /* alt_az.c */
double DLL_FUNC asine( const double arg);                    /* alt_az.c */
void DLL_FUNC full_ra_dec_to_alt_az( const DPT DLLPTR *ra_dec,
                DPT DLLPTR *alt_az,
                DPT DLLPTR *loc_epoch, const DPT DLLPTR *latlon,
                const double jd_utc, double DLLPTR *hr_ang);
void DLL_FUNC full_alt_az_to_ra_dec( DPT DLLPTR *ra_dec,
                              const DPT DLLPTR *alt_az,
                              const double jd_utc, const DPT DLLPTR *latlon);
void DLL_FUNC ra_dec_to_galactic( const double ra, const double dec,
                   double DLLPTR *glat, double DLLPTR *glon);
const double * DLL_FUNC j2000_to_galactic_matrix( void);
void DLL_FUNC galactic_to_ra_dec( const double glat, double glon,
                    double DLLPTR *ra, double DLLPTR *dec);
void DLL_FUNC ra_dec_to_supergalactic( const double ra, const double dec,
                   double DLLPTR *glat, double DLLPTR *glon);
const double * DLL_FUNC j2000_to_supergalactic_matrix( void);
void DLL_FUNC supergalactic_to_ra_dec( const double glat, double glon,
                    double DLLPTR *ra, double DLLPTR *dec);
void DLL_FUNC precess_pt( DPT DLLPTR *opt, const DPT DLLPTR *ipt,
                       double from_time, double to_time);
int DLL_FUNC get_comet_file( const char *cd_path,
                 const double year, const double mag_limit);
int DLL_FUNC extract_periodic_name( const char *istr, char *ostr);
double DLL_FUNC refraction( const double observed_alt);
double DLL_FUNC reverse_refraction( const double true_alt);
double DLL_FUNC saasta_refraction( const double observed_alt,
         const double pressure_mb, const double temp_kelvin,
         const double relative_humidity);
double DLL_FUNC reverse_saasta_refraction( const double true_alt,
         const double pressure_mb, const double temp_kelvin,
         const double relative_humidity);
double DLL_FUNC integrated_refraction( const double latitude,
                  const double observed_alt, const double wavelength_microns,
                  const double height_in_meters, const double rel_humid_pct,
                  const double temp_kelvins, const double pressure_mb);
double DLL_FUNC reverse_integrated_refraction( const double latitude,
                  const double refracted_alt, const double wavelength_microns,
                  const double height_in_meters, const double rel_humid_pct,
                  const double temp_kelvins, const double pressure_mb);

int DLL_FUNC calc_dist_and_posn_ang( const double *p1, const double *p2,
                                       double *dist, double *posn_ang);
void DLL_FUNC reverse_dist_and_posn_ang( double *to, const double *from,
                                 const double dist, const double posn_ang);

typedef struct
{
   double dX, dY, dPsi, dEps;    /* all in radians */
   double tdt_minus_ut1;         /* in seconds */
} earth_orientation_params;

#define EOP_FILE_NOT_FOUND     -1
#define EOP_FILE_WRONG_FORMAT  -2
#define EOP_ALLOC_FAILED       -3

extern int mjd_end_of_predictive_leap_seconds;       /* see delta_t.cpp */

int DLL_FUNC load_earth_orientation_params( const char *filename,
                                             int *file_date);
int DLL_FUNC get_earth_orientation_params( const double jd,  /* eop_prec.c */
                              earth_orientation_params *params,
                              int desired_params_mask);
int DLL_FUNC setup_precession_with_nutation_eops( double DLLPTR *matrix,
                    const double year);            /* eop_prec.c */
int64_t DLL_FUNC nanoseconds_since_1970( void);    /* nanosecs.c */
double DLL_FUNC current_jd( void);                 /* nanosecs.c */

#ifdef __cplusplus
}
#endif  /* #ifdef __cplusplus */
#endif  /* #ifndef AFUNCS_H_INCLUDED */
