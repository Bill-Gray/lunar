/* test_ref.cpp: test code for refraction

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
#include <stdlib.h>
#include "watdefs.h"
#include "afuncs.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int main( const int argc, const char **argv)
{
   int i, diff_mode = 0;
   double pressure_mb = 1013.,  temp_kelvin = 293., relative_humidity = .2;
   double wavelength_microns = .574, height_in_meters = 100.;
// extern double minimum_refractive_altitude;

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'p':
               pressure_mb = atof( argv[i] + 2);
               break;
            case 't':
               temp_kelvin = atof( argv[i] + 2) + 273.;
               break;
            case 'h':
               relative_humidity = atof( argv[i] + 2);
               break;
            case 'l':                  /* allow entry in nm */
               wavelength_microns = atof( argv[i] + 2) / 1000.;
               break;
            case 'a':
               height_in_meters = atof( argv[i] + 2);
               break;
            case 'd':
               diff_mode = 1;
               break;
            default:
               if( i > 1 || atof( argv[1]) == 0.)
                  {
                  printf( "? Didn't understand argument '%s'\n", argv[i]);
                  printf( "Command line arguments are:\n");
                  printf( "   -p(number)   Set pressure in millibars (default = 1013)\n");
                  printf( "   -t(temp)     Set temperature in degrees C (default = 20C)\n");
                  printf( "   -h(fraction) Set relative humidity fraction (default = .2)\n");
                  printf( "   -l(wavelen)  Set wavelength in nanometers (default = 574)\n");
                  printf( "   -a(ht)       Set altitude in meters (default = 100)\n");
                  return( -1);
                  }
               break;
            }
   printf( "Pressure %.2f millibars; temp %.1f C; relative humidity %.1f%%\n",
                   pressure_mb, temp_kelvin - 273., relative_humidity * 100.);
   printf( "Wavelength %.1f nm; altitude %.1f meters\n",
                   wavelength_microns * 1000., height_in_meters);
   for( i = (argc == 1 ? 0 : -1); i < 90; i++)
      if( i < 5 || i % 5 == 0)
         {
         const double rad_to_arcmin = 180. * 60. / PI;
         double observed_alt = (i == -1 ? atof( argv[1]) : (double)i) * PI / 180.;
         const double primitive_mul = (pressure_mb / 1010.) * (283. / temp_kelvin);
         double primitive = primitive_mul * refraction( observed_alt);
         double saasta = saasta_refraction( observed_alt, pressure_mb, temp_kelvin,
                             relative_humidity);
         double integrated =
               integrated_refraction( PI / 4., observed_alt, wavelength_microns,
                      height_in_meters, 100. * relative_humidity, temp_kelvin,
                      pressure_mb);

         if( i == -1 && observed_alt != 0.)
            {
            printf( "Primitive refraction:  %9.5f arcminutes\n",
                                       primitive * rad_to_arcmin) ;
            primitive = reverse_refraction( observed_alt - primitive);
            printf( "Primitive refraction reversed: %9.5f\n",
                         primitive * primitive_mul * rad_to_arcmin);
            printf( "Saasta refraction:     %9.5f\n", saasta * rad_to_arcmin);
            saasta = reverse_saasta_refraction( observed_alt - saasta,
                             pressure_mb, temp_kelvin,
                             relative_humidity);
            printf( "Saasta refraction (reversed):     %9.5f\n",
                             saasta * rad_to_arcmin);
            printf( "Integrated refraction: %9.5f' = %.3f\"\n",
                             integrated * rad_to_arcmin,
                             integrated * rad_to_arcmin * 60.);
            integrated = reverse_integrated_refraction( PI / 4.,
                              observed_alt - integrated, wavelength_microns,
                              height_in_meters, 100. * relative_humidity,
                              temp_kelvin, pressure_mb);
            printf( "Reverse integrated:    %9.5f' = %.3f\"\n",
                             integrated * rad_to_arcmin,
                             integrated * rad_to_arcmin * 60.);
//          printf( "Minimum altitude: %f meters\n",
//                           minimum_refractive_altitude);
            }
         else if( i > -1)
            {
            if( !i)
               printf( "Alt   Primit   Saasta   Integrated\n");
            if( diff_mode)
               {
               primitive -= integrated;
               saasta -= integrated;
               }
            printf( "%2d: %9.5f %9.5f %9.5f\n", i,
                          primitive * rad_to_arcmin,
                          saasta * rad_to_arcmin,
                          integrated * rad_to_arcmin);
            }
         }
   return( 0);
}
