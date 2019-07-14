/* See 'themis.cpp'.  The ephems for the Themis satellites cover about 18
days each.  My plan is to set up a cron job to run once a week or so,
downloading THEMIS ephems,  running through 'themis' to get those ephems
in the form used by 'eph2tle',  then using 'eph2tle' to generate 18 TLEs.

   Then we run this little program to enable TLEs to accumulate,  for
historical purposes.  Run as

./them_cat old_file.tle new_file.tle

   it will generate a new "temporary" TLE that contains data from 'old_file'
right up to the point where 'new_file' starts.  Then we finish off with the
'new_file' data.  If this all runs correctly,  we move the "temporary" file
over to replace the old file.       */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>

static FILE *err_fopen( const char *filename, const char *permits)
{
   FILE *rval = fopen( filename, permits);

   if( !rval)
      {
      fprintf( stderr, "'%s' not opened ", filename);
      perror( NULL);
      exit( -1);
      }
   return( rval);
}

int main( const int argc, const char **argv)
{
   FILE *old_file = err_fopen( argv[1], "rb");
   FILE *new_file = err_fopen( argv[2], "rb");
   FILE *temp_file = err_fopen( "ickywax.ugh", "wb");
   double jd0_new, end_jd_new;
   char buff[200], new_end[200];
   const time_t t0 = time( NULL);
   int i;
   bool found_end = false;

   assert( argc >= 3);
   while( fgets( buff, sizeof( buff), new_file)
                  && memcmp( buff, "# Ephem range:", 14))
      ;
   sscanf( buff + 15, "%lf %lf", &jd0_new, &end_jd_new);
   for( i = 0; i < 3; i++)       /* skip lines */
      if( !fgets( buff, sizeof( buff), new_file))
         perror( "Error reading 'new' file");
   strcpy( new_end, buff);

   while( !found_end && fgets( buff, sizeof( buff), old_file))
      {
      if( !memcmp( buff, "# Ephem range:", 14))
         sprintf( buff + 28, "%f %f\n", end_jd_new, 1.);
      if( !memcmp( buff, "# Ephemeris end:", 16))
         strcpy( buff, new_end);
      if( !memcmp( buff, "# Last updated with", 19))
         strcpy( buff + 42, ctime( &t0));

      fputs( buff,  temp_file);
      if( !memcmp( buff, "# MJD ", 6) && atof( buff + 6) == jd0_new - 1.)
         found_end = true;
      }
   assert( found_end);
   for( i = 0; i < 3 && fgets( buff, sizeof( buff), old_file); i++)
      fputs( buff, temp_file);
   fclose( old_file);
   while( fgets( buff, sizeof( buff), new_file))
      fputs( buff, temp_file);
   fclose( temp_file);
   fclose( new_file);
   unlink( argv[1]);
   rename( "ickywax.ugh", argv[1]);
   return( 0);
}
