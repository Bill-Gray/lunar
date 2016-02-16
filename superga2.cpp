/* superga2.cpp: code to compute a supergalactic-to-J2000 matrix

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

/* Running the following code produces a J2000 to supergalactic
   coordinate matrix,  which looks as follows:

 0.37501548  0.34135896  0.86188018
-0.89832046 -0.09572714  0.42878511
 0.22887497 -0.93504565  0.27075058

   Written in reply to an inquiry from Christopher Watson;  see

http://groups.yahoo.com/group/amastro/message/13964

   The code in 'supergal.cpp' gets a similar result,  using equatorial values
for the position of the supergalactic plane.  However,  the RA/dec values
given are not exact (given to one-second precision,  or 15-arcsec precision,
in RA.)  This method ought to be "exact".

   At the end,  it gives the coordinates for the SG north pole and zero
point in J2000 to full precision,  in decimal and base-60,  as

  2.82067484   02h 49m 14.4294s   (RA of zero point, SGL = SGB = 0)
 59.52834978  +59 31' 42.0592"    (dec of zero point, SGL = SGB = 0)
 18.91693936   18h 55m  0.9817s   (RA of SG north pole, SGB = 90)
 15.70893553  +15 42' 32.1679"    (dec of SG north pole)
*/

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define NORTH_POLE_LAT (6.32 * PI / 180.)
#define NORTH_POLE_LON (47.37 * PI / 180.)
#define ZERO_POINT_LON (137.37 * PI / 180.)

               /* The following matrix comes from _The Hipparcos & Tycho */
               /* Catalogues:  Introduction & Guide to the Data_, p 92:  */
static const double galactic_to_j2000[9] = {
      -.0548755604, -.8734370902, -.4838350155,
       .4941094279, -.4448296300,  .7469822445,
      -.8676661490, -.1980763734,  .4559837762 };

int main( const int unused_argc, const char **unused_argv)
{
   double matrix[9], rotation[9];
   int i, j;

            /* vector pointing toward (supergalactic) lat=lon=0: */
   rotation[0] = cos( ZERO_POINT_LON);
   rotation[1] = sin( ZERO_POINT_LON);
   rotation[2] = 0.;

            /* vector pointing toward the north supergalactic pole: */
   rotation[6] = cos( NORTH_POLE_LON) * cos( NORTH_POLE_LAT);
   rotation[7] = sin( NORTH_POLE_LON) * cos( NORTH_POLE_LAT);
   rotation[8] = sin( NORTH_POLE_LAT);

   rotation[3] = -cos( NORTH_POLE_LON) * sin( NORTH_POLE_LAT);
   rotation[4] = -sin( NORTH_POLE_LON) * sin( NORTH_POLE_LAT);
   rotation[5] = cos( NORTH_POLE_LAT);

   for( i = 0; i < 3; i++)
      for( j = 0; j < 3; j++)
         matrix[j + i * 3] =
                  rotation[i * 3] * galactic_to_j2000[j]
                + rotation[i * 3 + 1] * galactic_to_j2000[j + 3]
                + rotation[i * 3 + 2] * galactic_to_j2000[j + 6];

   printf( "%11.8lf %11.8lf %11.8lf\n", matrix[0], matrix[1], matrix[2]);
   printf( "%11.8lf %11.8lf %11.8lf\n", matrix[3], matrix[4], matrix[5]);
   printf( "%11.8lf %11.8lf %11.8lf\n", matrix[6], matrix[7], matrix[8]);
   printf( "\n");

   for( i = 0; i < 4; i++)
      {
      double oval, tval, sec;
      int deg, min;

      j = (i / 2) * 6;
      if( i % 2 == 0)
         oval = atan2( -matrix[j + 1], -matrix[j]) * 12. / PI + 12.;
      else
         oval = asin( matrix[j + 2]) * 180. / PI;
      deg = (int)oval;
      tval = (oval - (double)deg) * 60.;
      min = (int)tval;
      sec = (tval - (double)min) * 60.;
      printf( "%12.8lf   %02d %02d %7.4lf\n", oval, deg, min, sec);
      }
   return( 0);
}
