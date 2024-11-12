#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "mpc_func.h"
#ifdef __has_include
   #if __has_include("cgi_func.h")
       #include "cgi_func.h"
   #else
       #error   \
         'cgi_func.h' not found.  This project depends on the 'lunar'\
         library.  See www.github.com/Bill-Gray/lunar .\
         Clone that repository,  'make'  and 'make install' it.
       #ifdef __GNUC__
         #include <stop_compiling_here>
            /* Above line suppresses cascading errors. */
       #endif
   #endif
#else
   #include "cgi_func.h"
#endif

/* Code to pack/unpack designations for 'packed.htm' (q.v) via CGI.
Compiles with

gcc -Wall -Wextra -pedantic -o desigcgi desigcgi.c liblunar.a -lcurl -lm  */

int main( void)
{
   char field[30], buff[100], tbuff[100];
   int rval;

   printf( "Content-type: text/html\n\n");
   printf( "<pre>");
   avoid_runaway_process( 300);
   rval = initialize_cgi_reading( );
   if( rval <= 0)
      {
      printf( "<p> <b> CGI data reading failed : error %d </b>", rval);
      printf( "This isn't supposed to happen.  Please <a href='contact.htm'>report this.</a></p>\n");
      return( 0);
      }
   while( !get_cgi_data( field, buff, NULL, sizeof( buff)))
      {
      if( !strcmp( field, "packed") && strlen( buff) > 3)
         {
         rval = unpack_unaligned_mpc_desig( tbuff, buff);
         printf( "'%s' unpacks to '%s'\n", buff, tbuff);
         if( rval == -1)
            printf( "(This does not appear to be a valid packed designation)\n");
         }
      if( !strcmp( field, "unpacked") && strlen( buff) > 3)
         {
         rval = create_mpc_packed_desig( tbuff, buff);
         printf( "'%s' packs to '%s'\n", buff, tbuff);
         if( rval == -1)
            printf( "(This was not successfully packed.)\n");
         }
      }
   return( 0);
}
