#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include "mpc_func.h"

/* This code will read a file of observations containing (XXX) observations,
preceded by a line defining the lat/lon/alt as

COM Long. 239 18 45 E, Lat. 33 54 11 N, Alt. 100m, Google Earth

   or similar;  see 'mpc_func.cpp' for details.  (XXX) observations following
that line will be converted from (XXX) to roving observer data at (247).  To
be used for programs that know about (247) but don't know about the (XXX)
lat/lon/alt format.  (Which doesn't describe any of my programs.  But it
may be useful for other programs.  I may make this an on-line service;
paste in your COD XXX-style observations,  and it'll put out COD 247
roving observer observations.)  Compile with

gcc -Wall -Wextra -pedantic -Werror -o xxx2rov xxx2rov.c liblunar.a -lm  */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( argv[1], "rb");
   char buff[300];
   mpc_code_t cinfo;

   assert( argc == 2);
   assert( ifile);
   cinfo.code[0] = '\0';
   while( fgets( buff, sizeof( buff), ifile))
      {
      if( get_xxx_location_info( &cinfo, buff) == -2)
         {
         fprintf( stderr, "ERROR : malformed lat/lon/alt line\n%s\n", buff);
         exit( -1);
         }
      else if( !memcmp( buff, "COD XXX", 7))
         memcpy( buff + 4, "247", 3);
      else if( strlen( buff) > 80 && buff[80] < ' ' && !memcmp( buff + 77, "XXX", 3))
         {
         memcpy( buff + 77, "247", 3);
         buff[14] = 'V';
         printf( "%s", buff);
         buff[14] = 'v';
         sprintf( buff + 32, "1 %9.5f  %+9.5f   %4.0f",
                     cinfo.lon * 180. / PI, cinfo.lat * 180. / PI, cinfo.alt);
         buff[61] = ' ';
         }
      printf( "%s", buff);
      }
   fclose( ifile);
   return( 0);
}
