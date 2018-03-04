/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"

      /* Unit test code for COSPAR functions.  I've used this     */
      /* after making changes to 'cospar.txt' or 'cospar.cpp'     */
      /* just to verify that only the things that were _supposed_ */
      /* to change,  actually changed.                            */

int main( const int argc, const char **argv)
{
   double matrix[9], prev_matrix[9];
   int i, j, system_number, rval;
   clock_t t0 = clock( );

   if( argc == 2)          /* you can specify the COSPAR file name from */
      load_cospar_file( argv[1]);         /* the command line */
   setvbuf( stdout, NULL, _IONBF, 0);
   for( i = 0; i < 9; i++)
      prev_matrix[i] = 0.;
   for( i = -5; i < 2000; i++)
      {
      for( system_number = 0; system_number < 4; system_number++)
         {
         rval = calc_planet_orientation( i, system_number,
                        2451000. + (double)i * 1000., matrix);
         if( rval && rval != -1)
            printf( "rval %d\n", rval);
         else if( !rval)
            if( memcmp( matrix, prev_matrix, 9 * sizeof( double)))
               {
               printf( "Planet %d, system %d\n", i, system_number);
               for( j = 0; j < 9; j += 3)
                  printf( "%11.8f %11.8f %11.8f\n",
                             matrix[j], matrix[j + 1], matrix[j + 2]);
               for( j = 0; j < 9; j++)
                  prev_matrix[j] = matrix[j];
               }
         }
      }
   printf( "Total time: %f\n",
            (double)( clock() - t0) / (double)CLOCKS_PER_SEC);
   return( 0);
}
