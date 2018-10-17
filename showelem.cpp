/* showelem.cpp: "pretty formatting" of orbital elements

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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "watdefs.h"
#include "comets.h"
#include "date.h"
#include "afuncs.h"
#include "showelem.h"

/* 23 Sep 2006:  revised to show 'M1' and 'K1' for comets instead of
the (very wrong) 'H' and 'G',  and made sure those parameters were shown
for hyperbolic/parabolic orbits. */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

double DLL_FUNC decimal_day_to_dmy( double jd, long *year, int *month,
                              const int calendar);
int DLL_FUNC elements_in_mpc_format( char *obuff, const ELEMENTS *elem,
                const char *obj_id, const int is_cometary, const int format);

static double zero_to_two_pi( double ang)
{
   ang = fmod( ang, PI + PI);
   if( ang < 0.)
      ang += PI + PI;
   return( ang);
}

static void lop_digits( char *buff, int precision)
{
   while( precision && *buff >= '0' && *buff <= '9')
      {
      buff++;
      precision--;
      }
   while( *buff >= '0' && *buff <= '9')
      *buff++ = ' ';
}

double DLL_FUNC decimal_day_to_dmy( double jd, long *year, int *month,
                              const int calendar)
{
   int day;
   long floor_jd;

   jd += .5;
   floor_jd = (long)floor( jd);
   day_to_dmy( (long)floor_jd, &day, month, year, calendar);
   return( (double)day + jd - floor_jd);
}

      /* "normal_vector" takes a vector in ecliptic J2000 coordinates and */
      /* returns a _normalized_ vector rotated to _equatorial_ J2000. */
static double normal_vect( double *norm, const double *ival)
{
   double rval = 0.;
   int i;
   static const double sin_obliq_2000 = 0.397777155931913701597179975942380896684;
   static const double cos_obliq_2000 = 0.917482062069181825744000384639406458043;

   norm[0] = ival[0];
   norm[1] = ival[1] * cos_obliq_2000 - ival[2] * sin_obliq_2000;
   norm[2] = ival[2] * cos_obliq_2000 + ival[1] * sin_obliq_2000;
   for( i = 0; i < 3; i++)
      rval += ival[i] * ival[i];
   rval = sqrt( rval);
   if( rval != 0.)
      for( i = 0; i < 3; i++)
         norm[i] /= rval;
   return( rval);
}

static int show_formatted_dist( char *obuff, const double dist_in_au,
                     const int precision)
{
   int in_km = 0, i, n_digits = 0;

   if( dist_in_au < 0.)
      {
      *obuff++ = '-';
      return( show_formatted_dist( obuff, -dist_in_au, precision));
      }
   if( dist_in_au > 999999.)
      {
      strcpy( obuff, " <HUGE>");
      return( 0);
      }
   else if( dist_in_au > 9999.)
      sprintf( obuff, "%23.15f", dist_in_au);
   else if( dist_in_au > 999999. / AU_IN_KM)      /* within a million km */
      sprintf( obuff, "%23.18f", dist_in_au);
   else
      {
      sprintf( obuff, "%23.16f", dist_in_au * AU_IN_KM);
      in_km = 1;
      }
   for( i = 0; obuff[i] == ' '; i++)
      ;
   while( obuff[i++] != '.')
      n_digits++;
   if( n_digits < precision + 4)
      obuff[i + precision + 4 - n_digits] = '\0';
   else
      obuff[i - 1] = '\0';
   if( in_km)
      strcat( obuff, "km");
   return( 0);
}

static void add_pq_data( char *obuff, const double p, const double q,
                           const int precision)
{
   const int n_digits_to_show = (precision > 13 ? precision : 13);

   sprintf( obuff + strlen( obuff), "%+*.*f%+*.*f",
            n_digits_to_show + 3, n_digits_to_show, p,
            n_digits_to_show + 3, n_digits_to_show, q);
   lop_digits( obuff + 11, precision);
   lop_digits( obuff + 14 + n_digits_to_show,     precision);
   lop_digits( obuff + 17 + n_digits_to_show * 2, precision);
}


#define N_PERICENTERS      12

/* REMEMBER:  set 'central_obj', 'epoch', 'abs_mag', 'slope_param', 'gm' fields */

int DLL_FUNC elements_in_mpc_format( char *obuff, const ELEMENTS *elem,
                  const char *obj_id, const int is_cometary, const int format)
{
   const char *nineteen_blank_spaces = "                   ";
   double p_vect[3], q_vect[3];
   double dday;
   int month, i;
   long year;
   double perihelion_dist = elem->q;
   const double asc_node = zero_to_two_pi( elem->asc_node);
   const double arg_per = zero_to_two_pi( elem->arg_per);
   const int precision = (format & SHOWELEM_PRECISION_MASK);
   int n_lines = 0, n_digits_to_show;
   char *tptr;

   FSTRCPY( obuff, obj_id);
   obuff += strlen( obuff) + 1;
   n_lines++;
   if( elem->perih_time > 2. && elem->perih_time < 3000000.)
      {
      static const char *perinames[N_PERICENTERS] = {   "barion",
            "helion", "mercury", "venus", "gee", "mars", "jove",
            "saturn", "uranus", "neptune", "pluto", "lune" };
      const char *pericenter_name =
            (elem->central_obj >= -1 && elem->central_obj < N_PERICENTERS - 1) ?
            perinames[elem->central_obj + 1] : "----";

      dday = decimal_day_to_dmy( elem->perih_time, &year, &month,
                                             CALENDAR_JULIAN_GREGORIAN);
      sprintf( obuff, "   Peri%s %ld %s %.6f TT", pericenter_name, year,
              set_month_name( month, NULL), dday);
      if( format & SHOWELEM_PERIH_TIME_MASK)
         {
         char hhmmss[20];

         full_ctime( hhmmss, elem->perih_time,
                        FULL_CTIME_TIME_ONLY | CALENDAR_JULIAN_GREGORIAN);
         sprintf( obuff + strlen( obuff), " = %s (JD %.6f)",
                               hhmmss, elem->perih_time);
         }
      obuff += strlen( obuff) + 1;
      n_lines++;
      }

   dday = decimal_day_to_dmy( elem->epoch, &year, &month,
                                             CALENDAR_JULIAN_GREGORIAN);
   sprintf( obuff, "Epoch %4ld %s %9.6f TT = JDT %.6f", year,
              set_month_name( month, NULL), dday + 1.e-7, elem->epoch + 1.e-7);
                     /* lop off trailing zeroes after JD...: */
   for( i = 0; i < 5 && obuff[strlen( obuff) - (size_t)i - 1] == '0'; i++)
      ;
   obuff[strlen( obuff) - (size_t)i] = '\0';
                     /* ...and similar zeroes after the day: */
   tptr = strstr( obuff, " TT =") - i;
   memmove( tptr, tptr + i, strlen( tptr) + 1);

   obuff += strlen( obuff) + 1;
   n_lines++;
   n_digits_to_show = (precision > 10 ? precision : 10);
   if( is_cometary || elem->ecc >= 1.)
      {
      *obuff = 'q';
      show_formatted_dist( obuff + 1, perihelion_dist, precision);
      for( i = (int)strlen( obuff); i < n_digits_to_show + 9; i++)
         obuff[i] = ' ';
      obuff[n_digits_to_show + 9] = '\0';
      }
   else
      {
      const double mean_anomaly = zero_to_two_pi( elem->mean_anomaly);

      sprintf( obuff, "M%20.15f", mean_anomaly * 180. / PI);
      lop_digits( obuff + 9, precision);
      }
   if( !(format & SHOWELEM_OMIT_PQ_MASK))
      strcat( obuff, " (2000.0)            P               Q");
   obuff += strlen( obuff) + 1;
   n_lines++;

   normal_vect( p_vect, elem->perih_vec);
   normal_vect( q_vect, elem->sideways);
   if( is_cometary || elem->ecc >= 1.)
      {
      if( elem->abs_mag != 0.)
         {
         sprintf( obuff,
             (elem->is_asteroid ? "H%7.1f  G %4.2f   " :
                                  "M(T)%5.1f  K%5.1f  " ),
                                   elem->abs_mag + .05,
                                   elem->slope_param);
         if( !elem->is_asteroid)
            if( format & SHOWELEM_COMET_MAGS_NUCLEAR)
               obuff[2] = 'N';
         }
      else
         strcpy( obuff, nineteen_blank_spaces);
      for( i = 19; i < n_digits_to_show + 9; i++)
         obuff[i] = ' ';
      obuff[i] = '\0';
      }
   else
      {
      sprintf( obuff, "n%*.*f", n_digits_to_show + 8, n_digits_to_show + 3,
                              (180 / PI) / elem->t0);
      lop_digits( obuff + 9, precision);
      }
   obuff += strlen( obuff);

   sprintf( obuff, "Peri.%*.*f",
               n_digits_to_show + 6, n_digits_to_show, arg_per * 180. / PI);
   if( format & SHOWELEM_OMIT_PQ_MASK)
      obuff[11 + precision] = '\0';
   else
      add_pq_data( obuff, p_vect[0], q_vect[0], precision);
   obuff += strlen( obuff) + 1;
   n_lines++;

   if( is_cometary || elem->ecc >= 1.)
      strcpy( obuff, nineteen_blank_spaces);
   else
      {
      *obuff = 'a';
      show_formatted_dist( obuff + 1, elem->major_axis, precision);
      }
   for( i = (int)strlen( obuff); i < n_digits_to_show + 9; i++)
      obuff[i] = ' ';
   obuff[i] = '\0';
   obuff += strlen( obuff);

   sprintf( obuff, "Node %*.*f", n_digits_to_show + 6, n_digits_to_show,
               asc_node * 180. / PI);
   if( format & SHOWELEM_OMIT_PQ_MASK)
      obuff[11 + precision] = '\0';
   else
      add_pq_data( obuff, p_vect[1], q_vect[1], precision);
   obuff += strlen( obuff) + 1;
   n_lines++;

   if( is_cometary)
      sprintf( obuff, "e   1.0            ");
   else
      {
      sprintf( obuff, "e%*.*f", n_digits_to_show + 8,
                       n_digits_to_show + 3, elem->ecc);
      lop_digits( obuff + 8, precision);
      }
   obuff += strlen( obuff);

   sprintf( obuff, "Incl.%*.*f", n_digits_to_show + 6, n_digits_to_show,
                 elem->incl * 180. / PI);
   if( format & SHOWELEM_OMIT_PQ_MASK)
      obuff[11 + precision] = '\0';
   else
      add_pq_data( obuff, p_vect[2], q_vect[2], precision);
   obuff += strlen( obuff) + 1;
   n_lines++;

   if( !is_cometary && elem->ecc < 1.)
      {
      const double GAUSS_K = .01720209895; /* Gauss' gravitational constant */
      const double SOLAR_GM = (GAUSS_K * GAUSS_K);
      const double t0 = elem->major_axis *
                      sqrt( elem->major_axis * SOLAR_GM / elem->gm);
      const double t0_in_days = t0 * 365.25;
      const double apoapsis_dist =
             perihelion_dist * (1. + elem->ecc) / (1. - elem->ecc);
      char tbuff[40];

      if( elem->central_obj <= 0 || t0 > 1.)         /* helio or barycentric */
         {
         if( t0 > 1e+8 - 1.)        /* too big to fit in buffer */
            obuff += sprintf( obuff, "P!!!!!!!");
         else if( t0 > 9999.)
            obuff += sprintf( obuff, "P%7ld", (long)t0);
         else
            obuff += sprintf( obuff, "P%7.2f", t0);
         if( t0_in_days < 999.9)
            obuff += sprintf( obuff, "/%6.2fd ", t0_in_days);
         else
            obuff += sprintf( obuff, "         ");
         }
      else
         {
         if( t0_in_days * minutes_per_day < 9999.)  /* about 6.944 days */
            obuff += sprintf( obuff, "P%7.2fm/%5.3fd ",
                              t0_in_days * minutes_per_day,
                              t0_in_days);
         else
            obuff += sprintf( obuff, "P%7.2fd        ", t0_in_days);
         }
      if( elem->abs_mag != 0.)
         {
         sprintf( obuff, (elem->is_asteroid ? "  H%7.1f     G   %4.2f" :
                                              "  M(T)%5.1f    K %5.1f"),
                                              elem->abs_mag + .05,
                                              elem->slope_param);
         if( !elem->is_asteroid)
            if( format & SHOWELEM_COMET_MAGS_NUCLEAR)
               obuff[4] = 'N';
         obuff += strlen( obuff);
         }

      strcat( obuff, "   q ");
      show_formatted_dist( tbuff, perihelion_dist, precision);
      for( i = 0; tbuff[i] == ' '; i++)   /* skip leading spaces */
         ;
      strcat( obuff, tbuff + i);

      strcat( obuff, "  Q ");
      show_formatted_dist( tbuff, apoapsis_dist, precision);
      for( i = 0; tbuff[i] == ' '; i++)   /* skip leading spaces */
         ;
      strcat( obuff, tbuff + i);
      n_lines++;
      }
   return( n_lines);
}
