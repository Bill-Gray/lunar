/* alt_az.cpp: functions for coordinate conversions

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
#include <stdio.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define CONVERT (1000000. * 180. / PI)
#define TWO_PI (2. * PI)
#define J2000 2451545.0

static void ra_dec_to_alt_az( const double hr_ang, const double dec,
        double DLLPTR *alt, double DLLPTR *az, const double lat)
{
   double temp, cos_lat = cos( lat);

   *alt = asine( sin( lat) * sin( dec) + cos_lat * cos( dec) * cos( hr_ang));
   if( cos_lat < .00001)         /* polar case */
      *az = hr_ang;
   else
      {
      temp = (sin( dec) - sin( *alt) * sin( lat)) / (cos( *alt) * cos_lat);
      temp = PI - acose( temp);
      *az = ((sin( hr_ang) < 0.) ? temp : -temp);
      }
}

/* Normally,  the following will take the J2000 RA/dec and compute
the mean location at the epoch given by 'jd_utc'.  The result is
then stored in *loc_epoch.  However,  one can pass a NULL RA/dec;
in that case,  it's assumed that the location at epoch is already
stored in *loc_epoch.  The hour angle and nutation are then computed,
and the alt/azimuth.  You can optionally pass in NULLs for 'hr_ang'
and 'loc_epoch' in the (rather common) situations where you're
really just interested in the alt/az.   */

void DLL_FUNC full_ra_dec_to_alt_az( const DPT DLLPTR *ra_dec,
                DPT DLLPTR *alt_az,
                DPT DLLPTR *loc_epoch, const DPT DLLPTR *latlon,
                const double jd_utc, double DLLPTR *hr_ang)
{
   double ha, nutation_lon;
   const double t_centuries = (jd_utc - J2000) / 36525.;
   const double cos_obliq_2000 = 0.917482062069181825744000384639406458043;
   DPT loc_at_epoch;

   if( ra_dec)
      precess_pt( &loc_at_epoch, ra_dec, 2000., 2000. + t_centuries * 100.);
   else     /* no RA/dec at J2000 supplied */
      loc_at_epoch = *loc_epoch;
   ha = -loc_at_epoch.x - (green_sidereal_time( jd_utc) + latlon->x);
   nutation( t_centuries, &nutation_lon, NULL);
   ha -= cos_obliq_2000 * nutation_lon * (PI / 180.) / 3600.;
   ha = fmod( ha, TWO_PI);
   if( ha > PI) ha -= TWO_PI;
   if( ha <-PI) ha += TWO_PI;
   ra_dec_to_alt_az( ha, loc_at_epoch.y, &alt_az->y, &alt_az->x, latlon->y);
   if( hr_ang)
      *hr_ang = ha;
   if( loc_epoch)
      *loc_epoch = loc_at_epoch;
}

static void alt_az_to_ra_dec( double alt, double az,
                         double DLLPTR *hr_ang,
                         double DLLPTR *dec, const double lat)
{
   double temp, sin_dec, cos_lat = cos( lat);

   if( alt > PI / 2.)
      {
      alt = PI - alt;
      az += PI;
      }
   if( alt < -PI / 2.)
      {
      alt = -PI - alt;
      az -= PI;
      }
   sin_dec = sin( lat) * sin( alt) + cos_lat * cos( alt) * cos( az);
   *dec = asine( sin_dec);
   if( cos_lat < .00001)         /* polar case */
      *hr_ang = az + PI;
   else
      {
      temp = cos_lat * cos( *dec);
      temp = (sin( alt) - sin( lat) * sin_dec) / temp;
      temp = acose( -temp);
      if( sin( az) > 0.)
         *hr_ang = PI - temp;
      else
         *hr_ang = PI + temp;
      }
}

void DLL_FUNC full_alt_az_to_ra_dec( DPT DLLPTR *ra_dec,
                              const DPT DLLPTR *alt_az,
                              const double jd_utc, const DPT DLLPTR *latlon)
{
   double hr_ang, ra;
   DPT tmp;

   alt_az_to_ra_dec( alt_az->y, alt_az->x, &hr_ang,
                                         &tmp.y, latlon->y);
   ra = hr_ang + green_sidereal_time( jd_utc) + latlon->x;
   tmp.x = fmod( -ra, TWO_PI);
   precess_pt( ra_dec, &tmp, 2000. + (jd_utc - J2000) / 365.25, 2000.);
}

/* The following matrix was derived from the code in 'superga2.cpp'.  */

const double * DLL_FUNC j2000_to_supergalactic_matrix( void)
{
   static const double rval[9] = {
           0.37501548,  0.34135896,  0.86188018,
          -0.89832046, -0.09572714,  0.42878511,
           0.22887497, -0.93504565,  0.27075058 };

   return( rval);
}

void DLL_FUNC ra_dec_to_supergalactic( const double ra, const double dec,
                   double DLLPTR *glat, double DLLPTR *glon)
{
   double ipt[2], opt[2];

   ipt[0] = ra;
   ipt[1] = dec;
   precess_ra_dec( j2000_to_supergalactic_matrix( ), opt, ipt, 0);
   *glon = opt[0];
   *glat = opt[1];
}

void DLL_FUNC supergalactic_to_ra_dec( const double glat, double glon,
                    double DLLPTR *ra, double DLLPTR *dec)
{
   double ipt[2], opt[2];

   ipt[0] = glon;
   ipt[1] = glat;
   precess_ra_dec( j2000_to_supergalactic_matrix( ), opt, ipt, 1);
   *ra  = opt[0];
   *dec = opt[1];
}

const double * DLL_FUNC j2000_to_galactic_matrix( void)
{
               /* The following matrix comes from _The Hipparcos & Tycho */
               /* Catalogues:  Introduction & Guide to the Data_, p 92:  */
   static const double rval[9] = {
      -.0548755604, -.8734370902, -.4838350155,
       .4941094279, -.4448296300,  .7469822445,
      -.8676661490, -.1980763734,  .4559837762 };

   return( rval);
}

void DLL_FUNC ra_dec_to_galactic( const double ra, const double dec,
                   double DLLPTR *glat, double DLLPTR *glon)
{
   double ipt[2], opt[2];

   ipt[0] = ra;
   ipt[1] = dec;
   precess_ra_dec( j2000_to_galactic_matrix( ), opt, ipt, 0);
   *glon = opt[0];
   *glat = opt[1];
}

void DLL_FUNC galactic_to_ra_dec( const double glat, double glon,
                    double DLLPTR *ra, double DLLPTR *dec)
{
   double ipt[2], opt[2];

   ipt[0] = glon;
   ipt[1] = glat;
   precess_ra_dec( j2000_to_galactic_matrix( ), opt, ipt, 1);
   *ra  = opt[0];
   *dec = opt[1];
}

void DLL_FUNC precess_pt( DPT DLLPTR *opt, const DPT DLLPTR *ipt,
                        double from_time, double to_time)
{
   double precess[9];
   double temp_opt[2], temp_ipt[2];
   int dir = 0;

   if( from_time == to_time)
      {
      *opt = *ipt;
      return;
      }
   if( from_time == 2000.)
      {
      from_time = to_time;
      to_time = 2000.;
      dir = 1;
      }
   setup_precession( precess, from_time, to_time);
   temp_ipt[0] = -ipt->x;
   temp_ipt[1] = ipt->y;
   precess_ra_dec( precess, temp_opt, temp_ipt, dir);
   opt->x = -temp_opt[0];
   opt->y = temp_opt[1];
}
