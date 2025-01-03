/* jsats.cpp: functions for Galilean satellite posns

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

/* A little test program I wrote to test out my implementation of Lieske's E5
theory of the Galilean satellites,  as written in 'jsats.cpp'.  It works by
reading in an ASCII ephemeris of Jovicentric vectors for a given satellite,
computing the position for that satellite at that time using E5,  and showing
the difference,  both as an absolute XYZ in J2000 ecliptic coordinates and
in terms of "along-track" and "radial" components.  The latter allowed me to
see a long-term quadratic drift for Ganymede and Callisto.  (They were in
opposite directions,  so I'm confident that it's not a situation where I've
got the whole system rotating somehow.  My guess is that the mean longitudes
for those satellites should include a small quadratic term,  instead of
being just linear functions of time.)        */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "watdefs.h"
#include "lunar.h"
#include "afuncs.h"

#define JRADIUS_IN_KM 71418.0
// #define JRADIUS_IN_KM 71406.0

int main( const int argc, const char **argv)
{
   char buff[100];
   FILE *ifile;
   time_t t0 = time( NULL);
   const int sat_number = (argc == 2 ? atoi( argv[1]) : 0);

   if( sat_number < 1 || sat_number > 4)
      {
      printf( "'jsattest' needs a command-line argument from 1 to 4,\n");
      printf( "corresponding to the number of the Galilean satellite\n");
      printf( "that is being tested.\n");
      return( -1);
      }
   snprintf( buff, 7, "j%s.txt", argv[1]);
   ifile = fopen( buff, "rb");
   if( !ifile)
      {
      printf( "%s not opened\n", buff);
      return( -1);
      }
   printf( "All data in kilometers,  in J2000 ecliptic coords\n");
   printf( "Compiled %s %s; run %.24s\n", __DATE__, __TIME__, asctime( gmtime( &t0)));
   printf( "   JDE         dx        dy        dz             ");
   printf( "x          y          z         radial    along\n");
   while( fgets( buff, sizeof( buff), ifile))
      if( strlen( buff) > 56 && !memcmp( buff + 37, "00:00:00.0000 (CT)", 18))
         {
         const double jd = atof( buff );
         double loc[15], *tptr = loc + (sat_number - 1) * 3;
         double precess_matrix[9];
         double j2000_loc[3], x, y, z, r;

         if( !fgets( buff, sizeof( buff), ifile))
            {
            printf( "Read error");
            return( -2);
            }
         x = atof( buff);
         y = atof( buff + 24);
         z = atof( buff + 47);
         setup_ecliptic_precession( precess_matrix, 2000. + (jd - 2451545.) / 365., 2000.);
         calc_jsat_loc( jd, loc, 15, 0L);
         precess_vector( precess_matrix, tptr, j2000_loc);
         j2000_loc[0] *= JRADIUS_IN_KM;
         j2000_loc[1] *= JRADIUS_IN_KM;
         j2000_loc[2] *= JRADIUS_IN_KM;
         j2000_loc[0] -= x;
         j2000_loc[1] -= y;
         j2000_loc[2] -= z;
         r = sqrt( x * x + y * y);
         printf( "%10.2f %10.3f%10.3f%10.3f  %11.2f%11.2f%11.2f %10.2f %10.2f\n", jd,
              j2000_loc[0],
              j2000_loc[1],
              j2000_loc[2], x, y, z,
              (j2000_loc[0] * x + j2000_loc[1] * y) / r,
              (j2000_loc[1] * x - j2000_loc[0] * y) / r);
         }
   fclose( ifile);
   return( 0);
}
