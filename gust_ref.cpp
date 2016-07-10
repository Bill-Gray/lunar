/* gust_ref.cpp: function to compute Uranian satellite matrix

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

/* The original GUST86 code transformed Uranicentric coords to B1950.
Rather than do that,  then precess to J2000,  I wrote the snippet below
to compute the full Uranicentric-to-J2000 matrix.  Also,  the original
code did quite a bit of math to convert alf and del -- the RA and dec
of the Uranian pole in B1950 -- into a rotation matrix.  If you just
store the matrix as nine precomputed doubles,  you save yourself
some code and time and math.

   Spoiler alert:  the resulting output is

    {  0.9753206898, -0.2207422915,  0.0047321138},
    {  0.0619432123,  0.2529905682, -0.9654837185},
    {  0.2119259083,  0.9419493686,  0.2604204221},

   If you look at the gust86_posn() in gust86.c,  you will see this matrix.
So there is not really any reason to build this code.
*/

int main( const int unused_argc, const char **unused_argv)
{
// const double alf = 76.60666666666667 * DEGREES_TO_RADIANS;
// const double del = 15.03222222222222 * DEGREES_TO_RADIANS;
         // The GUST86 method originally used the above 'alf' and 'del'
         // values to specify the rotation from the plane of Uranus'
         // equator to the B1950 equator.  The 3x3 matrix was computed
         // within this code.  To simplify life a little,  I've set the
         // matrix to be a batch of 'const' values.
         // alf and delta aren't used here,  but their sines and cosines are...
   const double sin_alf = .9728028367170815;
   const double cos_alf = .2316347143137209;
   const double sin_del = .2593622252488415;
   const double cos_del = .9657801178912150;
         // ...to build the following 3x3 matrix to convert from Uranicentric
         // coords to B1950:
   const double trans[3][3] = {
           { sin_alf,                    -cos_alf,       0. },
           { cos_alf * sin_del, sin_alf * sin_del, -cos_del },
           { cos_alf * cos_del, sin_alf * cos_del,  sin_del } };

   // Rotation matrix for converting from B1950.0 FK4 positions
   // to J2000.0 FK5 positions.

   double P[3][3] =
   {
      { 0.9999256782, -0.0111820611, -0.0048579477 },
      { 0.0111820610,  0.9999374784, -0.0000271765 },
      { 0.0048579479, -0.0000271474,  0.9999881997 }
   };
   int i, j, k;

            /* Multiply trans by P and show the resulting matrix: */
   for( i = 0; i < 3; i++)
      for( j = 0; j < 3; j++)
         {
         double elem = 0;

         for( k = 0; k < 3; k++)
            elem += trans[i][k] * P[j][k];
         if( !j)
            printf( "    { ");
         printf( "%13.10lf%s", elem, (j == 2) ? "},\n" : ", ");
         }
}
