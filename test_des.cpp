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
   int n_errors_found = 0;
   int testing = 3, n_unpacked = 0, n_packed = 0;

   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      if( *buff != '#')
         {
         char tbuff[100];
         const int rval = unpack_mpc_desig( tbuff, buff);

         for( i = 0; buff[i] >= ' '; i++)
            ;
         buff[i] = '\0';
         if( testing & 1)
            {
            if( strcmp( buff + 16, tbuff) || rval != atoi( buff + 13))
               {
               n_errors_found++;
               printf( "UNPACKING MISMATCH\n'%s'\n'%s'\n", buff + 16, tbuff);
               }
            else
               n_unpacked++;
            }
         if( testing & 2)
            {
            create_mpc_packed_desig( tbuff, buff + 16);
            if( memcmp( tbuff, buff, 12))
               {
               n_errors_found++;
               printf( "PACKING MISMATCH\n'%.12s'\n'%.12s'\n", tbuff, buff);
               }
            else
               n_packed++;
            }
         }
      else if( !memcmp( buff, "# Test ", 7))
         testing = buff[7] - '0';
   printf( "%d packed correctly; %d unpacked correctly\n",
                  n_packed, n_unpacked);
   if( !n_errors_found)
      printf( "No errors found\n");
   return( 0);
}
