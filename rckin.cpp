#include <stdio.h>
#include <ctype.h>
/* rckin.cpp: processes/reformats JPL 'rckin' files

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

#include <string.h>
#include <stdlib.h>

/* Code to convert files such as 'rckin.ura091.txt',  etc.  from
   ftp://ssd.jpl.nasa.gov/pub/eph/satellites/rckin
   into the form used in 'rocks.cpp'.

   Of late,  JPL hasn't been updating these files.  I gather this is
because they are using integrated ephemerides for almost all "rocks".
So this may prove to be of historical interest only.  (Although when
I asked,  Bob Jacobson kindly provided updates to some of the files.)
*/

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( argv[1], "rb");
   char buff[200];
   char rock_name[80];

   if( !ifile)
      {
      printf( "%s not opened\n", argv[1]);
      return( -1);
      }

   while( fgets( buff, sizeof( buff), ifile))
      {
      char *tptr = strchr( buff + 19, '\'');

      if( tptr)
         *tptr = '\0';
      tptr = strchr( buff + 19, ',');
      if( tptr)
         *tptr = '\0';
      tptr = buff + 19;
      if( !memcmp( buff, "  RCKNAM", 8))
         {
         int i;

         strcpy( rock_name, tptr);
         for( i = 1; rock_name[i]; i++)
            rock_name[i] = tolower( rock_name[i]);
         }
      if( !memcmp( buff, "  RCKNUM", 8))
         printf( "\n   {  %d,             /* %s %s*/\n", atoi( buff + 17),
                    rock_name, (argc < 3 ? " " : argv[2]));
      if( !memcmp( buff, "  RCKEP", 7))
         {
         strcat( tptr, ",");
         printf( "      %-38s/* element epoch Julian date     */\n", tptr);
         }
      if( !memcmp( buff, "  RCKELT", 8))
         {
         const int field_no = atoi( buff + 9);
         const char *comment[10] = { NULL,
                      "a = semi-major axis (km)      ",
                      "h = e sin(periapsis longitude)",
                      "k = e cos(periapsis longitude)",
                      "l = mean longitude (deg)      ",
                      "p = tan(i/2) sin(node)        ",
                      "q = tan(i/2) cos(node)        ",
                      "apsis rate (deg/sec)          ",
                      "mean motion (deg/sec)         ",
                      "node rate (deg/sec)           " };

         if( field_no == 4 || field_no == 7 || field_no == 8 || field_no == 9)
            strcat( tptr, " * PI / 180.");
         strcat( tptr, ",");
         printf( "      %-38s/* %s*/\n", tptr,
                     comment[atoi( buff + 9)]);
         }
      if( !memcmp( buff, "  CTRRA", 7))
         {
         strcat( tptr, " * PI / 180.,");
         printf( "      %-38s/* Laplacian plane pole ra (deg) */\n", tptr);
         }
      if( !memcmp( buff, "  CTRDEC", 8))
         {
         strcat( tptr, " * PI / 180. },");
         printf( "      %-38s/* Laplacian plane pole dec (deg)*/\n", tptr);
         }
      }
   fclose( ifile);
   return( 0);
}
