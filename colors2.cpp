/* colors2.cpp: functions for color conversions

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


int tycho_to_johnson_colors( const double bt_minus_vt, double *results);
double johnson_b_minus_v_from_tycho( double b_v_t);

/* The following function takes a (B-V)T value,  i.e.,  a B-V color in
the Tycho magnitude system,  and computes from it V-VT,  the difference
between Johnson V and Tycho V;  the corresponding B-V in the Johnson
scheme;  and V-Hp,  the difference between Johnson V and Hp (Hipparcos
"magnitude") for the input color.

   I got the raw data for this from Brian Skiff,  who posted them on
the Minor Planet Mailing List (MPML) with the following comment:

   "...To get standard V and B from Tycho-2, it is probably best to use the
relation shown by Mike Bessell in the July 2000 PASP... Bessell does not give
an algebraic relation, but instead shows a cubic spline fit with a look-up
table...  I have copied out Bessell's table below as a flat ASCII list."

   In what follows,  values for (B-V)T away from the lookup table
points is computed using (again) a cubic spline.  The functions vary
with sufficient slowness that this ought to be accurate down to the
.001 mag level.

    References :

https://groups.io/g/mpml/message/771

Bessell 2000, PASP 112, 961 (July),  RELATION BETWEEN BT-VT AND
HIPPARCOS/TYCHO DATA FOR B-G MAIN-SEQUENCE STARS AND K-M GIANTS   */

#define LOOKUP_SIZE 46

int tycho_to_johnson_colors( const double bt_minus_vt, double *results)
{
   int i, table_loc = (int)( (bt_minus_vt + .25) / .05) - 1;
   double dx, coeff[4];
   static const short lookup_tbl[LOOKUP_SIZE * 3] = {
      /*  BT-VT     V-VT  del(B-V)  V-Hp */
      /* -0.250 */   38,   31,     66,
      /* -0.200 */   30,   21,     51,
      /* -0.150 */   22,   11,     36,
      /* -0.100 */   15,    5,     21,
      /* -0.050 */    8,    2,      6,
      /* -0.000 */    1,   -5,    -11,
      /*  0.050 */   -5,  -10,    -25,
      /*  0.100 */  -12,  -17,    -38,
      /*  0.150 */  -18,  -20,    -48,
      /*  0.200 */  -24,  -21,    -58,
      /*  0.250 */  -29,  -23,    -69,
      /*  0.300 */  -35,  -25,    -79,
      /*  0.350 */  -40,  -25,    -87,
      /*  0.400 */  -45,  -26,    -94,
      /*  0.450 */  -50,  -30,   -101,
      /*  0.500 */  -54,  -35,   -108,
      /*  0.550 */  -59,  -45,   -114,
      /*  0.600 */  -64,  -51,   -120,
      /*  0.650 */  -68,  -60,   -127,
      /*  0.700 */  -72,  -68,   -131,
      /*  0.750 */  -77,  -76,   -134,
      /*  0.800 */  -81,  -85,   -137,
      /*  0.850 */  -85,  -94,   -142,
      /*  0.900 */  -89, -104,   -147,
      /*  0.950 */  -93, -113,   -151,
      /*  1.000 */  -98, -122,   -155,
      /*  1.050 */ -102, -131,   -158,
      /*  1.100 */ -106, -142,   -157,
      /*  1.150 */ -110, -154,   -160,
      /*  1.200 */ -115, -166,   -162,
      /*  1.250 */ -119, -178,   -164,
      /*  1.300 */ -124, -189,   -166,
      /*  1.350 */ -128, -199,   -166,
      /*  1.400 */ -133, -210,   -165,
      /*  1.450 */ -138, -222,   -164,
      /*  1.500 */ -143, -234,   -161,
      /*  1.550 */ -148, -245,   -157,
      /*  1.600 */ -154, -256,   -153,
      /*  1.650 */ -160, -266,   -148,
      /*  1.700 */ -165, -277,   -143,
      /*  1.750 */ -172, -288,   -137,
      /*  1.800 */ -178, -299,   -131,
      /*  1.850 */ -185, -309,   -125,
      /*  1.900 */ -191, -320,   -119,
      /*  1.950 */ -199, -331,   -112,
      /*  2.000 */ -206, -342,   -106 };

   if( bt_minus_vt < -.25 || bt_minus_vt > 2.)
      return( -1);         /* out of table range */
   if( table_loc < 0)
      table_loc = 0;
   if( table_loc >= LOOKUP_SIZE - 4)
      table_loc = LOOKUP_SIZE - 4;
   dx = ((bt_minus_vt + .25) / .05) - (double)table_loc;
   coeff[0] = (dx - 1.) * (dx - 2.) * (dx - 3.) / -6.;
   coeff[1] = dx * (dx - 2.) * (dx - 3.) / 2.;
   coeff[2] = dx * (dx - 1.) * (dx - 3.) / -2.;
   coeff[3] = dx * (dx - 1.) * (dx - 2.) / 6.;
   for( i = 0; i < 3; i++)
      {
      const short *tptr = lookup_tbl + i + table_loc * 3;

      results[i] = .001 * ((double)tptr[0] * coeff[0]
                         + (double)tptr[3] * coeff[1]
                         + (double)tptr[6] * coeff[2]
                         + (double)tptr[9] * coeff[3]);
      }
   results[1] += bt_minus_vt;  /* change a 'delta' into an 'absolute' */
   return( 0);
}

#ifdef TEST_FUNC
#include <stdio.h>
#include <stdlib.h>

   /* This function derived from data on p 57,  _Intro & Guide to the Data_ */
double johnson_b_minus_v_from_tycho( double b_v_t)
{
   double delta = 0.;

   if( b_v_t < -.2 || b_v_t > 1.8)
      return( 99.);        /* no reasonable transformation possible */
   if( b_v_t < .1)
      delta = -.006 + .006 * (b_v_t + .2) / .3;
   else if( b_v_t < .5)
      delta = .046 * (b_v_t - .1) / .4;
   else if( b_v_t < 1.4)
      delta = .046 - .054 * (b_v_t - .5) / .9;
   else if( b_v_t < 1.8)
      delta = -.008 - .024 * (b_v_t - 1.4) / .4;
   return( .85 * b_v_t + delta);
}

int main( const int argc, const char **argv)
{
   double ovals[3];

   if( argc == 2 && !tycho_to_johnson_colors( atof( argv[1]), ovals))
      {
      printf( "V-VT = %.4f\n", ovals[0]);
      printf( "B-V = %.4f\n", ovals[1]);
      printf( "B-V = %.4f (from original formula)\n",
                               johnson_b_minus_v_from_tycho( atof( argv[1])));
      printf( "V-Hp = %.4f\n", ovals[2]);
      }
}
#endif
