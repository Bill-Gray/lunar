/* lunar2.cpp: functions for modest-precision lunar coords

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

/* Implements a simplified lunar ephemeris via the method described
in Meeus' _Astronomical Algorithms_.  The actual series coefficients
are stored in 'vsop.bin',  which must be read into a buffer before
you call these functions.  At the time I wrote all this -- early
1990s -- memory was a scarce resource,  so I packed bytes as much
as possible.  I suppose all this would still be a good idea on some
embedded systems.  It's probably not the way I would do things now.
The code is not easy to follow.  But it _does_ all work,  and runs
fast and has a small footprint.

   The lunar longitude/distance terms each consume twelve bytes in the
binary file:  one each for a d, m, mp, and f coefficient,  and then the
longitude (sl) and radius (sr) amplitudes are stored as 32-bit integers.

   The lunar latitude terms are similar,  but only eight bytes:  still
the d, m, mp, and f coefficients,  but then there's just a latitude
term stored as a 32-bit integer.  */

#define LON_R_TERM_SIZE          12u
#define LAT_TERM_SIZE             8u
#define N_TERMS                  60
#define LUNAR_LON_DIST_OFFSET    59354u
#define LUNAR_LAT_OFFSET         LUNAR_LON_DIST_OFFSET + LON_R_TERM_SIZE * N_TERMS
#define LUNAR_FUND_OFFSET        LUNAR_LAT_OFFSET + LAT_TERM_SIZE * N_TERMS

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "watdefs.h"
#include "lunar.h"
#include "get_bin.h"

#define Lp  fund[0]
#define D   fund[1]
#define M   fund[2]
#define Mp  fund[3]
#define F   fund[4]
#define A1  fund[5]
#define A2  fund[6]
#define A3  fund[7]
#define T   fund[8]

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int DLL_FUNC lunar_fundamentals( const void FAR *data, const double t,
                                        double DLLPTR *fund)
{
   int i, j;
   const char FAR *tptr = (const char FAR *)data + LUNAR_FUND_OFFSET;
   double tpow;

   assert( get32bits( tptr) == 0x6ed5a0b1);
   assert( get32bits( data) == 0x00260000);
   for( i = 0; i < 5; i++)
      {
      fund[i] = get_double( tptr);
      tptr += 8;
      tpow = t;
      for( j = 4; j; j--, tpow *= t, tptr += 8)
         fund[i] += tpow * get_double( tptr);
      }

   A1 = 119.75 + 131.849 * t;
   A2 =  53.09 + 479264.290 * t;
   A3 = 313.45 + 481266.484 * t;
   T = t;
   for( i = 0; i < N_FUND - 1; i++)      /* convert to radians */
      {
      fund[i] = fmod( fund[i], 360.);
      if( fund[i] < 0.) fund[i] += 360.;
      fund[i] *= PI / 180.;
      }
   return( 0);
}

int DLL_FUNC lunar_lon_and_dist( const void FAR *data, const double DLLPTR *fund,
                 double DLLPTR *lon, double DLLPTR *r, const long precision)
{
   int i, j;
   const signed char *tptr = (const signed char *)data + LUNAR_LON_DIST_OFFSET;
   double sl_sum = 0., sr_sum = 0., e;

   assert( get32bits( tptr) == 0x00010000);
   e = 1. - .002516 * T - .0000074 * T * T;
   for( i = N_TERMS; i; i--, tptr += LON_R_TERM_SIZE)
      {
      const int32_t sl = get32bits( tptr + 4);
      const int32_t sr = get32bits( tptr + 8);

      if( labs( sl) > precision || labs( sr) > precision)
         {
         double term;
         const signed char d = tptr[0], m = tptr[1], mp = tptr[2], f = tptr[3];
         const double arg = (double)d * D + (double)m * M
                          + (double)mp * Mp + (double)f * F;

         if( sl)
            {
            term = (double)sl * sin( arg);
            for( j = abs( m); j; j--)
               term *= e;
            sl_sum += term;
            }
         if( sr)
            {
            term = (double)sr * cos( arg);
            for( j = abs( m); j; j--)
               term *= e;
            sr_sum += term;
            }
         }
      }
   if( precision < 3959L)
      sl_sum += 3958. * sin( A1) + 1962. * sin( Lp - F) + 318. * sin( A2);
   *lon = (Lp * 180. / PI) + sl_sum * 1.e-6;
   while( *lon < 0.)
      *lon += 360.;
   while( *lon > 360.)
      *lon -= 360.;
   *r = 385000.56 + sr_sum / 1000.;
   return( 0);
}

double DLL_FUNC lunar_lat( const void FAR *data, const double DLLPTR *fund,
                                           const long precision)
{
   int i, j;
   const signed char *tptr = (const signed char FAR *)data + LUNAR_LAT_OFFSET;
   double rval = 0., e;

   assert( get32bits( tptr) == 0x01000000);
   e = 1. - .002516 * T - .0000074 * T * T;
   for( i = N_TERMS; i; i--, tptr += LAT_TERM_SIZE)
      {
      const int32_t sb = get32bits( tptr + 4);

      if( labs( sb) > precision)
         {
         double term;
         const signed char d = tptr[0], m = tptr[1], mp = tptr[2], f = tptr[3];
         const double arg = (double)d * D + (double)m * M
                          + (double)mp * Mp + (double)f * F;

         term = (double)sb * sin( arg);
         for( j = abs( m); j; j--)
            term *= e;
         rval += term;
         }
      }
   if( precision < 2236L)
      rval += -2235. * sin( Lp) + 382. * sin( A3) + 175. * sin( A1 - F) +
               175. * sin( A1 + F) + 127. * sin(Lp - Mp) - 115. * sin(Lp+Mp);
   return( rval * 1.e-6);
}
