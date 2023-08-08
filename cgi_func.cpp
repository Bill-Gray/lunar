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
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "cgi_func.h"
#include "watdefs.h"

/* Server-based versions of Find_Orb,  sat_id,  and astcheck (see
'fo_serve.cpp',  'sat_id2.cpp',  cgicheck.cpp' respectively) all use
get_multipart_form_data( ) to parse CGI data.  I use multipart for
these because those services require the option of file uploads.

   The server-based code to list GNSS (navigation) satellites,
or satellite elements for a specific date,  however,  can (and do)
use enctype="application/x-www-form-urlencoded".  This is a little
simpler and makes it easy for people to construct URLs to access
particular bits.

   All of my server-based code uses avoid_runaway_process() to ensure that
if something went wrong and the code got hung up,  it would abort after a
decent length of time and give an error message.  (Which really should be
different for each program;  at present,  you get the Find_Orb message.)
*/

#if defined( __linux) || defined( __unix__) || defined( __APPLE__)
#include <sys/time.h>         /* these allow resource limiting */
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>

static const char *email_address_mangled_for_no_spam =
   "p&#x202e;&ocirc;&#xe7;.&ouml;tulp&#x165;c&eacute;j&ocirc;&#x159;p&#x40;ot&uacute;l&#x202c;m";
         /* pluto at projectpluto address with diacritical marks and bidi  */
         /* text... should baffle all but the most devoted spammers.  Put */
         /* your mangled address here.   */

static void sighandler( int signum, siginfo_t *info, void *unused_ucontext)
                                 __attribute__ ((noreturn));

static void sighandler( int signum, siginfo_t *info, void *unused_ucontext)
{
   INTENTIONALLY_UNUSED_PARAMETER( unused_ucontext);
   if( signum == SIGXCPU)
      {
      printf( "\n\n<h1> Ran out of time </h1>\n");
      printf( "<p> Find_Orb was unable to determine an orbit for these\n");
      printf( "observations.  That may just be because the problem exceeded\n");
      printf( "the capabilities of this server;  to avoid runaway processes,\n");
      printf( "there's an intentional fifteen-second cap on getting a solution.\n");
      printf( "<p> You should probably try again with a shorter\n");
      printf( "arc,  or else use the 'native' Windows or Linux or OS/X\n");
      printf( "flavor of Find_Orb.\n");
      }
   else
      {
      printf( "\n\n<h1> Unknown signal </h1>\n");
      printf( " <p> Got signal %d from process %ld.\n", signum,
                      (unsigned long)info->si_pid);
      printf( "I'm still learning a bit about these signals.  Please\n");
      printf( "copy/paste this message and e-mail it to\n");
      printf( "%s\n", email_address_mangled_for_no_spam);
      }
   exit( 0);
}

void avoid_runaway_process( const int max_time_to_run)
{
   struct rlimit r;
   struct sigaction act;

   r.rlim_cur = max_time_to_run;
   r.rlim_max = max_time_to_run + 5;
   setrlimit( RLIMIT_CPU, &r);

   memset(&act, 0, sizeof(act));

   act.sa_sigaction = sighandler;
   act.sa_flags = SA_SIGINFO;
   sigaction(SIGXCPU, &act, NULL);
}
#else
   /* In Win32,  you can't do this sort of process limiting     */
   /* (I think),  so we just have a dummy function.             */

void avoid_runaway_process( const int max_time_to_run)
{
   INTENTIONALLY_UNUSED_PARAMETER( max_time_to_run);
}
#endif         /* _WIN32 */

static int get_urlencoded_piece( const char **idata,
              char *buff, size_t max_buff, const char end_char)
{
   int c;
   const char *tptr = *idata;

   max_buff--;
   while( (c = *tptr++) > 13 && max_buff-- && c != end_char)
      {
      if( c == '+')
         c = ' ';
      else if( c == '%')
         {
         int i, c1;

         c = 0;
         for( i = 0; i < 2; i++)
            {
            c *= 16;
            c1 = *tptr++;
            if( c1 >= '0' && c1 <= '9')
               c += c1 - '0';
            else if( c1 >= 'A' && c1 <= 'F')
               c += c1 - 'A' + 10;
            else                    /* wasn't an hex digit: */
               return( -1);         /* not supposed to happen */
            }
         }
      *buff++ = (char)c;
      }
   if( c != end_char)
      tptr--;
   *buff =  '\0';
   *idata = tptr;
   return( c);
}

int get_urlencoded_form_data( const char **idata,
                              char *field, const size_t max_field,
                              char *buff, const size_t max_buff)
{
   int c;

   if( get_urlencoded_piece( idata, field, max_field, '=') != '=')
      return( -1);
   c = get_urlencoded_piece( idata, buff, max_buff, '&');
   return( c != '&' && c > 13);
}

#define OVERRUN         -3

int get_multipart_form_data( const char *boundary, char *field,
                char *buff, char *filename, const size_t max_len)
{
   char *tptr, *endptr;
   size_t bytes_read = 0, blen = 0;

   if( filename)
      *filename = '\0';
   while( boundary[blen] >= ' ')
      blen++;
   if( fgets( buff, (int)max_len, stdin)
                  && (tptr = strstr( buff, "name=\"")) != NULL
                  && (endptr = strchr( tptr + 6, '"')) != NULL)
      {
      char *filename_ptr = strstr( tptr, "filename=\"");

      if( strlen( buff) == max_len - 1)
         return( OVERRUN);
      *endptr = '\0';
      strcpy( field, tptr + 6);
      if( filename && filename_ptr
             && (endptr = strchr( filename_ptr + 10, '"')) != NULL)
         {
         *endptr = '\0';
         strcpy( filename, filename_ptr + 10);
         }
      if( fgets( buff, (int)max_len, stdin))
         {
         if( strlen( buff) == max_len - 1)
            return( OVERRUN);
         while( bytes_read < max_len - 1 &&
                 fgets( buff + bytes_read, (int)( max_len - bytes_read), stdin)
                 && memcmp( buff + bytes_read, boundary, blen))
            bytes_read += strlen( buff + bytes_read);
         }
      }
   else
      return( -1);
   while( bytes_read && (buff[bytes_read - 1] == 10
                             || buff[bytes_read - 1] == 13))
      bytes_read--;
   buff[bytes_read] = '\0';
   return( (int)bytes_read);
}

/* Overall,  we have three situations :

   (1) If the QUERY_STRING environment variable is set,  our data
comes from the URL (i.e.,  METHOD=GET).  We can parse the
QUERY_STRING data with get_urlencoded_form_data() (q.v.).

   (2) Otherwise,  we look at the CONTENT_TYPE environment variable.
If it's

CONTENT_TYPE=application/x-www-form-urlencoded

   we can just read in a line from stdin and parse that,  just as
if it had been the QUERY_STRING variable.  Other than the source of
the text coming from stdin instead of an environment variable,  it's
exactly the same procedure.

   (3) If,  instead,  the CONTENT_TYPE variable looks like

CONTENT_TYPE=multipart/form-data; boundary=--------...more_text_here

   then we have to use the get_multipart_form_data( ) function above.

   (4) If it's none of these,  something else is going on and this
code (at present) doesn't know what to do about it. */

static char *url_string = NULL;
static const char *url_pointer = NULL;
static char *ibuff = NULL, *boundary = NULL;
static int method, ilen;

#define METHOD_GET               1
#define METHOD_PUT_URL           2
#define METHOD_PUT_MULTIPART     3

int initialize_cgi_reading( void)
{
   method = 0;
   url_string = getenv( "QUERY_STRING");
   if( url_string && *url_string)      /* METHOD=GET */
      {
      url_pointer = url_string;
      method = METHOD_GET;
      }
   else
      {
      const char *eptr = getenv( "CONTENT_TYPE");

      if( !eptr)
         return( -1);
      if( !memcmp( eptr, "application/x-www-form-urlencoded", 33))
         method = METHOD_PUT_URL;
      else if( !memcmp( eptr, "multipart/form-data;", 20))
         method = METHOD_PUT_MULTIPART;
      else
         return( -2);
      eptr = getenv( "CONTENT_LENGTH");
      assert( eptr);
      ilen = atoi( eptr);
      assert( ilen);
      ibuff = (char *)malloc( ilen + 1);
      if( !fgets( ibuff, ilen + 1, stdin))
         return( -3);
      if( method == METHOD_PUT_URL)
         url_pointer = ibuff;
      else           /* must be multipart */
         {
         boundary = (char *)malloc( strlen( ibuff) + 1);
         strcpy( boundary, ibuff);
         }
      }
   return( method);
}

int get_cgi_data( char *field, char *data, char *filename, const size_t max_data)
{
   int rval;
   const size_t max_field = 30;

   switch( method)
      {
      case METHOD_GET:
      case METHOD_PUT_URL:
         if( filename)
            *filename = '\0';
         rval = get_urlencoded_form_data( &url_pointer, field, max_field,
                                 data, max_data);

         break;
      case METHOD_PUT_MULTIPART:
         rval = get_multipart_form_data( boundary, field, data,
                                          filename, max_data);
         if( rval > 0)
            rval = 0;
         break;
      default:
         rval = -1;
         break;
      }
   return( rval);
}

void free_cgi_data( void)
{
   if( ibuff)
      free( ibuff);
}
