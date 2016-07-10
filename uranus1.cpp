/* uranus1.cpp: test code for Uranian satellite coordinate functions,
   linked "normally" (compare to uranus2.cpp)

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

/*
** Uranus.cpp
** Functions relating to Uranus and its satellites.  Written by Chris
** Marriott for SkyMap,  with some modifications by me (Bill J. Gray).
**
** Created:   17-JUL-98
**
** $History: Uranus.cpp $
**
** *****************  Version 1  *****************
** User: Chris        Date: 2/10/98    Time: 6:55
** Created in $/SkyMap 4.0
** Moved from solar system DLL into main project.
**
** *****************  Version 2  *****************
** User: Chris        Date: 18/07/98   Time: 4:49p
** Updated in $/SkyMap 4.0/SolarSys
** New file to compute data relating to Uranus and its satellites.
*/

#include "gust86.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


//   UranusSatellite
//   Compute the positions of Uranus's satellites.

#define AU_IN_KM 1.495978707e+8

static void subtract_test_data( double *ivals, const int satellite,
                           const int is_velocity)
{
            /* Following are state vectors for all five satellites */
            /* for JDE 2451539.5 = 27 Dec 1999,  as given by Horizons. */

#ifdef COMPARE_TO_HORIZONS
   static const double test_data[60] = {
   /* X            Y            Z            vx           vy           vz */
   130609.219205, -63609.804468, 123914.883497,
                         332755.476940,   19395.593811, -339740.309375,
   169031.741269, -89894.366595, 186084.712454,
                         299664.428730,    7278.709236, -267485.829018,
  -270571.386731, -33418.471015, 341847.941085,
                        237817.716674,  -102816.931847,  177706.820567,
   528174.319532, -51132.708805,-242110.833416,
                       -100437.670099,    87763.820120, -237466.213432,
  -128420.756530,  15056.079281,  12341.354683,
                         33339.723882, -163288.106059,  552424.102310 };
#else
            /* ...and as given by this program.  If you run */
            /* 'uranus1 2451539.5',  you should get the following.  Run */
            /* 'uranus1 2451539.5 z',  and the following will be subtracted */
            /* from the output,  and you should get all zeroes.         */
   static const double test_data[60] = {
    130608.740570,  -63609.454817,  123914.157536,
                        332753.924857,   19395.739133, -339739.097154,
    169031.139042,  -89894.351535,  186085.561827,
                        299665.211058,    7278.266868, -267484.726136,
   -270565.130051,  -33419.941877,  341847.964177,
                        237818.078944, -102815.738968,  177702.440261,
    528175.409357,  -51133.015840, -242110.654650,
                       -100437.576682,   87763.875080, -237466.804419,
   -128422.243906,   15056.226282,   12341.164393,
                         33339.740820, -163282.893418,  552432.372423 };
#endif
   int i;

   for( i = 0; i < 3; i++)
      {
      double test_value = test_data[satellite * 6 + is_velocity * 3 + i];

      test_value /= AU_IN_KM;       /* cvt km to AU... */
      if( is_velocity)
         test_value /= 86400.;      /* and days to seconds */
      ivals[i] -= test_value;
      }
}

int main( int argc, char **argv)
{
   const int test_differences = (argc == 1 || argc == 3);
   const double jde = (argc == 1 ? 2451539.5 : atof( argv[1]));
   int nSat;                     // loop counter
   const char *sat_names[5] = {
         "Ariel", "Umbriel", "Titania", "Oberon", "Miranda" };

   // Process each satellite in turn.

   for (nSat=0; nSat<5; nSat++)
   {
      double dRect[6];               // satellite coordinates

      // Retrieve the coordinates of the satellite relative to the
      // centre of Saturn. These are equatorial coordinates for the
      // mean ecliptic and epoch of J2000.0.  Positions in AU,
      // velocities in AU/second.  Printed out in km and km/s.

      gust86_posn( jde, nSat, dRect );
      if( test_differences)
         subtract_test_data( dRect, nSat, 0);
      printf( "%d %-8s: %14.6f %14.6f %14.6f\n", nSat, sat_names[nSat],
                                               dRect[0] * AU_IN_KM,
                                               dRect[1] * AU_IN_KM,
                                               dRect[2] * AU_IN_KM);
      if( test_differences)
         subtract_test_data( dRect + 3, nSat, 1);
      printf( "            %14.6f %14.6f %14.6f\n",
                                              dRect[3] * AU_IN_KM,
                                              dRect[4] * AU_IN_KM,
                                              dRect[5] * AU_IN_KM);
   }
   return( 0);
}
