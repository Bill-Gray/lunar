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
#include <stdlib.h>
#include <string.h>
#include "cgi_func.h"
#include "watdefs.h"
#include "stringex.h"

/* Code to invoke the 'astcheck' routine from an HTML form.
You'll see a _lot_ of overlap between this and 'sat_id2.cpp',
the code to do satellite checking from an HTML form,  and the
'fo_serve.cpp' code to run Find_Orb from an HTML form...
there's a lot of recycling going on.

   In setting up cgicheck on your server,  you'll need to make
sure that 'cgicheck.htm' is present,  of course,  along with a copy
of ObsCodes.html for parallax data.  You should have either mpcorb.dat
or the Lowell astorb.dat files,  or both,  present (note that the code
expects the lowercase filenames).  The results are similar with either
database of orbital elements,  except that 'astorb' gives you current
ephemeris uncertainties.
*/

int astcheck_main( const int argc, const char **argv);    /* astcheck.c */

int main( const int unused_argc, const char **unused_argv)
{
   const char *argv[20];
   const size_t max_buff_size = 40000;       /* room for 500 obs */
   char *buff;
   char boundary[100], field[30];
   const char *temp_obs_filename = "temp_obs.txt";
   int argc = 2;
   FILE *lock_file = fopen( "lock.txt", "w");
   size_t bytes_written = 0;
#ifndef _WIN32
   extern char **environ;
#endif
   extern int verbose;
   double search_radius = 2.;    /* default to looking two degrees */

   INTENTIONALLY_UNUSED_PARAMETER( unused_argv);
   INTENTIONALLY_UNUSED_PARAMETER( unused_argc);
#ifndef _WIN32                   /* If things take more than 60 seconds, */
   avoid_runaway_process( 60);   /* assume failure and give an error msg */
#endif         /* _WIN32            to that effect                       */
   printf( "Content-type: text/html\n\n");
   printf( "<html> <body> <pre>\n");
   if( !lock_file)
      {
      printf( "Server is busy.  (At least as currently written,  we can only run\n");
      printf( "a single batch of observations at a time.  However,  a given batch of\n");
      printf( "observations <i>usually</i> doesn't take all that long... come back in\n");
      printf( "a minute or two and try again.)\n");
      printf( "</pre> </body> </html>");
      return( 0);
      }
   fprintf( lock_file, "We're in\n");
#ifndef _WIN32
   for( size_t i = 0; environ[i]; i++)
      fprintf( lock_file, "%s\n", environ[i]);
#endif
   if( !fgets( boundary, sizeof( boundary), stdin))
      {
      printf( "<b> No info read from stdin</b>");
      printf( "This isn't supposed to happen.\n");
      return( 0);
      }
   buff = (char *)malloc( max_buff_size + 100);
   while( get_multipart_form_data( boundary, field, buff, NULL, max_buff_size) >= 0)
      {
      if( !strcmp( field, "TextArea") || !strcmp( field, "upfile"))
         {
         if( strlen( buff) > 70)
            {
            FILE *ofile = fopen( temp_obs_filename,
                               (bytes_written ? "ab" : "wb"));

            bytes_written += fwrite( buff, 1, strlen( buff), ofile);
            fclose( ofile);
            }
         }
      else if( !strcmp( field, "radius"))
         {
         const char *verbosity = strchr( buff, 'v');

         search_radius = atof( buff);
         if( verbosity)
            verbose = atoi( verbosity + 1) + 1;
         }
      else if( !strcmp( field, "catalog"))
         {
         if( buff[0] == '1')              /* use MPCORB */
            argv[argc++] = "-M";
         }
      else if( !strcmp( field, "uncertainties"))
         argv[argc++] = "-e";
      else if( !strcmp( field, "unn_only"))
         argv[argc++] = "-u";
      }
   free( buff);
   if( verbose)
      printf( "Searching to %f degrees;  %u bytes read from input\n",
                     search_radius, (unsigned)bytes_written);
   argv[0] = "cgicheck";
   argv[1] = temp_obs_filename;
   snprintf_err( field, sizeof( field), "-r%.2f", search_radius * 3600.);  /* cvt degrees to arcsec */
   argv[argc++] = field;
   argv[argc] = NULL;
   astcheck_main( argc, argv);
   printf( "</pre> </body> </html>");
   return( 0);
}
