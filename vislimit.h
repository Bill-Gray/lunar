/* vislimith: header file for visibility/sky brightness computations
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

#define BRIGHTNESS_DATA struct brightness_data

#pragma pack(4)
BRIGHTNESS_DATA
   {
               /* constants for a given time: */
   double zenith_ang_moon, zenith_ang_sun, moon_elongation;
   double ht_above_sea_in_meters, latitude;
   double temperature_in_c, relative_humidity;
   double year, month;
               /* values varying across the sky: */
   double zenith_angle;
   double dist_moon, dist_sun;         /* angular,  not real linear */
   int mask;   /* indicates which of the 5 photometric bands we want */
               /* Items computed in set_brightness_params: */
   double air_mass_sun, air_mass_moon, lunar_mag;
   double k[5], c3[5], c4[5];
   double ka[5];        /* Aerosol extinction coeffs in UBVRI */
   double kr[5];        /* Rayleigh extinction coeffs */
   double ko[5];        /* Ozone extinction */
   double kw[5];        /* Water extinction */
   double year_term;    /* variation due to 11-year solar cycle */
               /* Items computed in compute_limiting_mag: */
   double air_mass_gas, air_mass_aerosol, air_mass_ozone;
   double extinction[5];
               /* Internal parameters from compute_sky_brightness: */
   double air_mass;
   double brightness[5];
   };
#pragma pack( )

#ifdef _WIN32
#define DLL_FUNC __stdcall
#else
#define DLL_FUNC
#endif

#ifdef __cplusplus
extern "C" {
#endif

int DLL_FUNC set_brightness_params( BRIGHTNESS_DATA *b);
int DLL_FUNC compute_sky_brightness( BRIGHTNESS_DATA *b);
double DLL_FUNC compute_limiting_mag( BRIGHTNESS_DATA *b);
int DLL_FUNC compute_extinction( BRIGHTNESS_DATA *b);

#ifdef __cplusplus
}
#endif
