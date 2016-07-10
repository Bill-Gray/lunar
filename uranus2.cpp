/* uranus2.cpp: test program for Uranian satellite functions,
   loading those functions from a DLL (compare to uranus1.cpp)

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
** This is almost exactly identical to 'uranus1.cpp'.  The only difference
** is that in this code,  I experimented with dynamic loading of the DLL.
** That is to say,  one can (in Windows) compile this code without having
** it be statically linked to the Uranian satellite library.
**
** The benefit of this comes in my Guide software.  Here,  the user would
** almost never really be interested in Uranian satellites;  those objects
** could easily never need to be computed,  and loading up code for them
** would be a waste of computer cycles.  The following demonstrates how one
** can access the Uranian ephemerides on an "as needed" basis.
**
** As it turns out,  I've yet to take advantage of this particular capability.
** But figuring out how it worked has proven useful in other cases;  Guide
** does not (for example) default to having access to the rather large
** ASCOM (scope control) or SDP4 (artificial satellite ephemeris) code;
** instead,  those are loaded only when one actually wants to control a
** telescope through ASCOM or examines artificial satellites.
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "windows.h"

//   UranusSatellite
//   Compute the positions of Uranus's satellites.

#define AU_IN_KM 1.495978707e+8

static void subtract_test_data( double *ivals, const int satellite,
                           const int is_velocity)
{
            /* Following are state vectors for all five satellites */
            /* for JDE 2451539.5 = 27 Dec 1999,  as given by Horizons... */

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

static HINSTANCE load_gust86_lib( const int unload)
{
   static HINSTANCE h_gust86_lib = (HINSTANCE)0;
   static int first_time = 1;

   if( unload)
      {
      if( h_gust86_lib)
         FreeLibrary( h_gust86_lib);
      h_gust86_lib = NULL;
      first_time = 1;
      }
   else if( first_time)
      {
      h_gust86_lib = LoadLibrary( "gust86.dll");
      first_time = 0;
      }
   return( h_gust86_lib);
}

typedef void (__stdcall *gust86_fn)( const double jde, const int isat, double *r);

int gust86_posn( const double jde, const int isat, double r[6])
{
   static gust86_fn func = (gust86_fn)0;
   int rval = 0;
   HINSTANCE h_gust86_lib;

   if( !r)         /* flag to unload library */
      {
      load_gust86_lib( -1);
      return( 0);
      }
   h_gust86_lib = load_gust86_lib( 0);
   if( h_gust86_lib && !func)
      func = (gust86_fn)GetProcAddress( h_gust86_lib, (LPCSTR)1);
   if( func)
      (*func)( jde, isat, r);
   else
      rval = -1;
   if( !h_gust86_lib)
      rval = -2;
   return( rval);
}

int main( const int argc, const char **argv)
{
   const int test_differences = (argc == 1 || argc == 3);
   const double jde = (argc == 1 ? 2451539.5 : atof( argv[1]));
   int nSat;                     // loop counter

   // Process each satellite in turn.

   setvbuf( stdout, NULL, _IONBF, 0);
   for (nSat=0; nSat<5; nSat++)
   {
      double dRect[6];               // satellite coordinates

      // Retrieve the coordinates of the satellite relative to the
      // centre of Saturn. These are equatorial coordinates for the
      // mean ecliptic and epoch of J2000.0.  Positions in AU,
      // velocities in AU/second.  Printed out in km and km/s.

      memset( dRect, 0, 6 * sizeof( double));
      gust86_posn( jde, nSat, dRect );

      if( test_differences)
         subtract_test_data( dRect, nSat, 0);
      printf( "%d: %14.6lf %14.6lf %14.6lf\n", nSat, dRect[0] * AU_IN_KM,
                                               dRect[1] * AU_IN_KM,
                                               dRect[2] * AU_IN_KM);
      if( test_differences)
         subtract_test_data( dRect + 3, nSat, 1);
      printf( "   %14.6lf %14.6lf %14.6lf\n", dRect[3] * AU_IN_KM,
                                              dRect[4] * AU_IN_KM,
                                              dRect[5] * AU_IN_KM);
   }
   return( 0);
}
