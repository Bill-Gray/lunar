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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "watdefs.h"
#include "lunar.h"

int main( const int argc, const char **argv)
{
   int i;
   double loc[3], jd;

   if( argc != 2)
      {
      printf( "'ssattest' takes a JD on the command line,  and outputs\n");
      printf( "coordinates for the eight main satellites of Saturn.\n");
      return( -1);
      }
   jd = atof( argv[1]);
   printf( "Date: %.5f\n", jd);
   for( i = 0; i < 8; i++)
      {
      jd = atof( argv[1]);
      calc_ssat_loc( jd, loc, i, 0L);
      printf( "%d: %9.6f %9.6f %9.6f\n", i, loc[0] * 100., loc[1] * 100.,
                  loc[2] * 100.);
      }
   return( 0);
}
