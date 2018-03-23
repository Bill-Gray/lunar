/* mpcorb.cpp: extracts orbital elements from 'mpcorb.dat' or 'astorb.dat'

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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "watdefs.h"
#include "comets.h"
#include "date.h"

/* Tools for extracting orbital elements from 'mpcorb.dat'-style files. */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

         /* Many values in MPC files are stored in an 'extended hex'     */
         /* scheme: usual '0'-'9' = 0 to 9,  'A-F' = 10 to 15.  But      */
         /* then 'G'-'Z' = 16 to 35 and 'a' to 'z' = 36 to 61.           */

static int extended_hex_to_int( const char ichar)
{
   int rval;

   if( ichar >= '0' && ichar <= '9')
      rval = ichar - '0';
   else if( ichar >= 'A' && ichar <= 'Z')
      rval = ichar - 'A' + 10;
   else    /* if( ichar >= 'a' && ichar <= 'z') */
      rval = ichar - 'a' + 10 + 26;
   return( rval);
}

static long extract_mpc_epoch( const char *epoch_buff)
{
   const long year = 100 * extended_hex_to_int( epoch_buff[0])
                    + 10 * extended_hex_to_int( epoch_buff[1])
                       +   extended_hex_to_int( epoch_buff[2]);
   const int month = extended_hex_to_int( epoch_buff[3]);
   const int day = extended_hex_to_int( epoch_buff[4]);

   return( dmy_to_day( day, month, year, CALENDAR_GREGORIAN));
}

static void do_remaining_element_setup( ELEMENTS *elem)
{
   double ma;

   elem->mean_anomaly  *= PI / 180.;
   elem->arg_per       *= PI / 180.;
   elem->asc_node      *= PI / 180.;
   elem->incl          *= PI / 180.;
   elem->q = elem->major_axis * (1. - elem->ecc);
   derive_quantities( elem, SOLAR_GM);
   elem->angular_momentum = sqrt( SOLAR_GM * elem->q * (1. + elem->ecc));
   ma = elem->mean_anomaly;
   if( ma > PI)            /* find _nearest_ perihelion time */
      ma -= PI + PI;
   elem->perih_time = elem->epoch - ma * elem->t0;
   elem->is_asteroid = 1;
   elem->central_obj = 0;
   elem->gm = SOLAR_GM;
}

/* extract_mpcorb_dat( ) extracts the data from an 'mpcorb.dat'-formatted
buffer,  putting the orbital elements into 'elem' if 'elem' is non-NULL.
The epoch of the elements is returned.  (You can call the function with
elem=NULL to verify that the buffer is in the right format,  or to find
out what the epoch is.)
   Note that there's a very similar function for 'astorb'-formatted
data below.       */

long extract_mpcorb_dat( ELEMENTS *elem, const char *buff)
{
   long epoch_jd = 0;

   if( strlen( buff) > 200 && buff[47] == ' ' && buff[82] == '.' &&
                buff[25] == ' ' && buff[29] == '.' && buff[36] == ' ')
      {
      epoch_jd = extract_mpc_epoch( buff + 20);
      if( elem)        /* just doing a format check */
         {
         elem->epoch = (double)epoch_jd - .5;
         elem->mean_anomaly = atof( buff + 26);
         elem->arg_per      = atof( buff + 37);
         elem->asc_node     = atof( buff + 48);
         elem->incl         = atof( buff + 59);
         elem->ecc          = atof( buff + 69);
         elem->major_axis   = atof( buff + 92);
         do_remaining_element_setup( elem);
         if( buff[10] == '.')
            elem->abs_mag = atof( buff + 8);
         else
            elem->abs_mag = 0.;
         if( buff[16] == '.')
            elem->slope_param = atof( buff + 14);
         else
            elem->slope_param = 0.;
         }
      }
   return( epoch_jd);
}

/* extract_astorb_dat( ) is essentially the same as the above
   extract_mpcorb_dat( ) function,  except using 'astorb.dat'-formatted
   data.       */

static long extract_long( const char *ibuff)
{
   long rval = 0;

   while( *ibuff == ' ')      /* skip leading spaces */
      ibuff++;
   while( *ibuff != ' ')
      {
      if( *ibuff != '.')
         rval = rval * 10 + (long)( *ibuff - '0');
      ibuff++;
      }
   return( rval);
}

long extract_astorb_dat( ELEMENTS *elem, const char *buff)
{
   long epoch_jd;
   ELEMENTS telem;

   if( strlen( buff) > 267
       && sscanf( buff + 41, "%lf %lf", &telem.abs_mag, &telem.slope_param) == 2)
      {
      epoch_jd = atoi( buff + 106);
      telem.mean_anomaly = (double)extract_long( buff + 115) / 1.e+6;
      telem.arg_per      = (double)extract_long( buff + 126) / 1.e+6;
      telem.asc_node     = (double)extract_long( buff + 137) / 1.e+6;
      telem.incl         = (double)extract_long( buff + 148) / 1.e+6;
      telem.ecc          = (double)extract_long( buff + 160) / 1.e+8;
      telem.major_axis = atof( buff + 169);
      epoch_jd = dmy_to_day( epoch_jd % 100L,           /* day */
                             (epoch_jd / 100L) % 100L,  /* month */
                             epoch_jd / 10000L,         /* year */
                             CALENDAR_GREGORIAN);
      telem.epoch = (double)epoch_jd - .5;
      do_remaining_element_setup( &telem);
      if( elem)
         *elem = telem;
      }
   else
      epoch_jd = 0;
   return( epoch_jd);
}
