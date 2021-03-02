#include <stdio.h>
#include <stdlib.h>

#define CSI "\x1b\x5b"

/* Code to extract chunks of the light pollution map at

https://www.nasa.gov/feature/goddard/2017/new-night-lights-maps-open-up-possible-real-time-applications/

   and output a roughly 120x120-km (40x40) pixel area as an image
on an xterm,  using color control codes.

   The 8 MByte JPG is 13500 x 6750 pixels,  roughly three km to the pixel :

https://www.nasa.gov/specials/blackmarble/2016/globalmaps/BlackMarble_2016_3km.jpg

   Run 'convert BlackMarble_2016_3km.jpg /tmp/night.pgm',  and you have
something this program can use.  I doubt that anyone else is apt to use
this code,  but it's helped me to visualize the map.  */

int main( const int argc, const char **argv)
{
   const int xsize = 13500, ysize = xsize / 2;
   FILE *ifile = fopen( "/tmp/night.pgm", "rb");
   int x0, y0;
   int x, y, size = (argc > 3 ? atoi( argv[3]) : 20);
   long header_size = 18;
   unsigned char buff[200];

   if( !ifile)
      printf( "This program gets data from /tmp/night.pgm,  and cannot\n"
              "open that file.  See 'pollute.c' for details.\n");
   if( argc < 3)
      printf( "Usage : ./pollute longitude latitude\n");
   if( !ifile || argc < 3)
      return( -1);
   x0 = (int)( ( 180. + atof( argv[1])) * (double)xsize / 360.);
   y0 = (int)( ( 90. - atof( argv[2])) * (double)ysize / 180.);
   x0 %= xsize;
   for( y = y0 - size; y < y0 + size; y++)
      {
      long loc = header_size + (long)y * (long)xsize + (long)( x0 - size);

      fseek( ifile, loc, SEEK_SET);
      if( fread( (char *)buff, size * 2, 1, ifile))
         {
         for( x = 0; x < size * 2; x++)
            if( x == size && y == y0)
               printf( CSI "48;2;255;128;0m  ");
            else
               printf( CSI "48;2;%d;%d;%dm  ", buff[x], buff[x], buff[x]);
         printf( CSI "49m\n");     /* reset default bkgrnd */
         }
      }
   return( -1);
}
