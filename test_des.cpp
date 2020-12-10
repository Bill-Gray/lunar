#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "mpc_func.h"

/* Unit test code for packing and unpacking MPC designations. Runs
through all designations in 'test_des.txt' both ways to make sure the
functions in 'mpc_fmt.cpp' and 'unpack.cpp' return the expected output. */

int main( void)
{
   FILE *ifile = fopen( "test_des.txt", "rb");
   char buff[200];
   size_t i;

   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      if( *buff != '#')
         {
         char tbuff[100];
         const int rval = unpack_mpc_desig( tbuff, buff);

         for( i = 0; buff[i] >= ' '; i++)
            ;
         buff[i] = '\0';
         if( strcmp( buff + 16, tbuff) || rval != atoi( buff + 13))
            printf( "UNPACKING MISMATCH\n'%s'\n'%s'\n", buff + 16, tbuff);
         create_mpc_packed_desig( tbuff, buff + 16);
         if( memcmp( tbuff, buff, 12))
            printf( "PACKING MISMATCH\n'%.12s'\n'%.12s'\n", tbuff, buff);
         }
   return( 0);
}
