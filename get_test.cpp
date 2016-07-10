/* get_test.cpp: tests/exercises the get_time_from_string( ) function

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
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "watdefs.h"
#include "date.h"

/* Test routine for get_time_from_string().  This function is supposed
to puzzle out a huge range of bizarre inputs,  such as "March 4" (fourth
day of March) vs. "4 Mar" and so on.  'get_test.txt' contains some example
inputs;  'get_test' reads those inputs,  prints out the resulting times,
and you can then check to make sure that the approved output resulted
by comparing it to 'get_out.txt'.  Exactly one line -- the one corresponding
to 'now-4d' -- should differ. */

int main( int argc, char **argv)
{
   FILE *ifile = fopen( argc == 1 ? "get_test.txt" : argv[1], "rb");

   if( ifile)
      {
      char buff[80];
      double jd = 0;
      int i;
      unsigned time_format = 0;

      setvbuf( stdout, NULL, _IONBF, 0);
      while( fgets( buff, sizeof( buff), ifile))
         {
         for( i = 0; buff[i] >= ' '; i++)
            ;
         buff[i] = '\0';
         if( !memcmp( buff, "format", 6))
            {
            sscanf( buff + 6, "%x", &time_format);
            printf( "%s\n", buff);
            }
         else if( *buff != ';')
            {
            char obuff[80];
            int is_ut;

            jd = get_time_from_string( jd, buff, (int)time_format, &is_ut);
            full_ctime( obuff, jd, FULL_CTIME_MILLISECS | CALENDAR_JULIAN_GREGORIAN);
            printf( "%s %d %s\n", obuff, is_ut, buff);
            }
         else
            printf( "%s\n", buff);
         }
      fclose( ifile);
      }
   return( ifile ? 0 : -1);
}
