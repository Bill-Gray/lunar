#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifdef __WATCOMC__
   #include <stdbool.h>
#endif

#include "stringex.h"

#ifdef NO_BUILT_IN_SNPRINTF
     /* MSVC/C++ (up to VS2015) has no snprintf.  Yes,  you read that     */
     /* correctly.  MSVC has an _snprintf which doesn't add a '\0' at the */
     /* end if max_len bytes are written.  You can't pass a NULL string   */
     /* to determine the necessary buffer size.  The following, however,  */
     /* is a "good enough" replacement:  for a non-NULL string, the       */
     /* output will be "correct" (the '\0' will always be present and in  */
     /* the right place).  The only deviation from "proper" sprintf() is  */
     /* that the maximum return value is max_len;  you can never know how */
     /* many bytes "should" have been written.                            */
     /*   Somewhat ancient MSVCs don't even have vsnprintf;  you have to  */
     /* use vsprintf and work things out so you aren't overwriting the    */
     /* end of the buffer.  So be advised:  if using old MSVCs,  be sure  */
     /* your buffer is big enough for the 'full' result,  even though it  */
     /* will then get truncated correctly.                                */

int snprintf( char *string, const size_t max_len, const char *format, ...)
{
   va_list argptr;
   int rval;

   va_start( argptr, format);
#if defined(_MSC_VER) && _MSC_VER <= 1100
   rval = vsprintf( string, format, argptr);
#else
   rval = vsnprintf( string, max_len, format, argptr);
#endif
   va_end( argptr);
   return( rval);
}
#endif    /* #ifdef _MSC_VER */

/* I find myself frequently snprintf()-ing at the end of a string,  with
something like snprintf( str + strlen( str), sizeof( str) - strlen(str), ...).
snprintf_append( ) aborts if the buffer would overflow;  it is used if we
think the output should never be that big.  snprintf_append_trunc( ) simply
stops at max_len bytes,  with a '\0' terminator at max_len - 1,  and is
used if truncation may legitimately be needed. */

static int _sn_append( const bool no_truncation, char *string,
                        const size_t max_len, const char *format, va_list argptr)
{
   int rval;
   const size_t ilen = strlen( string);

   assert( ilen <= max_len);
#if defined(_MSC_VER) && _MSC_VER <= 1100
   rval = vsprintf( string + ilen, format, argptr);
#else
   rval = vsnprintf( string + ilen, max_len - ilen, format, argptr);
#endif
   va_end( argptr);
   if( rval < 0 || rval + ilen >= max_len)
      if( no_truncation)
         {
         fprintf( stderr, "snprintf_append: %ld/%ld/%ld\n%s\n%s\n",
                     (long)ilen, (long)rval, (long)max_len, string, format);
         assert( 0);
         }
   return( rval + (int)ilen);
}

int snprintf_append_trunc( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...)
{
   va_list argptr;

   va_start( argptr, format);
   return( _sn_append( false, string, max_len, format, argptr));
}

int snprintf_append( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...)
{
   va_list argptr;

   va_start( argptr, format);
   return( _sn_append( true, string, max_len, format, argptr));
}

#if !defined( __APPLE__) && !defined( __OpenBSD__)

/* strlcat() and strlcpy() appear in some BSDs,  but I don't think they
appear anyplace else (though in my opinion,  they should).  In my
humble opinion,  they should only be used when truncation may
legitimately happen.  If strlen( src) should never exceed dsize-1,
use strlcat_err() or strlcpy_err(),  shown below.  (If 'dest' is an
automatic variable,  such that dsize == sizeof( dest),  use the
strlcat_error() or strlcpy_error() macros;  see stringex.h.)  Otherwise,
the result will be silently truncated and you won't know there was a
problem (until,  probably,  a bit further down the line).

The following is from
http://www.openbsd.org/cgi-bin/cvsweb/src/lib/libc/string/ . */

/* $OpenBSD: strlcpy.c,v 1.16 2019/01/25 00:19:25 millert Exp $   */

/*
 * Copyright (c) 1998, 2015 Todd C. Miller <millert@openbsd.org>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * Copy string src to buffer dst of size dsize.  At most dsize-1
 * chars will be copied.  Always NUL terminates (unless dsize == 0).
 * Returns strlen(src); if retval >= dsize, truncation occurred.
 */
size_t strlcpy( char *dst, const char *src, const size_t dsize)
{
   const char *osrc = src;
   size_t nleft = dsize;

   /* Copy as many bytes as will fit. */
   if (nleft != 0) {
      while (--nleft != 0) {
         if ((*dst++ = *src++) == '\0')
            break;
      }
   }

   /* Not enough room in dst, add NUL and traverse rest of src. */
   if (nleft == 0) {
      if (dsize != 0)
         *dst = '\0';      /* NUL-terminate dst */
      while (*src++)
         ;
   }

   return(src - osrc - 1); /* count does not include NUL */
}

/*
 * Appends src to string dst of size dsize (unlike strncat, dsize is the
 * full size of dst, not space left).  At most dsize-1 characters
 * will be copied.  Always NUL terminates (unless dsize <= strlen(dst)).
 * Returns strlen(src) + MIN(dsize, strlen(initial dst)).
 * If retval >= dsize, truncation occurred.
 */

size_t strlcat( char *dst, const char *src, const size_t dsize)
{
   const char *odst = dst;
   const char *osrc = src;
   size_t n = dsize;
   size_t dlen;

   /* Find the end of dst and adjust bytes left but don't go past end. */
   while (n-- != 0 && *dst != '\0')
      dst++;
   dlen = dst - odst;
   n = dsize - dlen;

   if (n-- == 0)
      return(dlen + strlen(src));
   while (*src != '\0') {
      if (n != 0) {
         *dst++ = *src;
         n--;
      }
      src++;
   }
   *dst = '\0';

   return(dlen + (src - osrc));  /* count does not include NUL */
}
#endif    /* #if !defined( __APPLE__) && !defined( __OpenBSD__) */

/* Same as strlcpy() and strlcat(),  except that if truncation occurs,
we abort.  strlcpy()/strlcat() should only be used in situations where
truncation is entirely to be expected;  if truncation indicates a bug,
use the _err() versions instead.  */

size_t strlcpy_err( char *dst, const char *src, const size_t dsize)
{
   const size_t rval = strlcpy( dst, src, dsize);

   if( rval >= dsize)
      {
      fprintf( stderr, "strlcpy overflow: dsize = %ld, rval %ld, '%s', '%s'\n",
                     (long)dsize, (long)rval, dst, src);
      assert( rval >= dsize);
      exit( -1);
      }
   return( rval);
}

size_t strlcat_err( char *dst, const char *src, const size_t dsize)
{
   const size_t rval = strlcat( dst, src, dsize);

   if( rval >= dsize)
      {
      fprintf( stderr, "strlcat overflow: dsize = %ld, rval %ld, '%s', '%s'\n",
                     (long)dsize, (long)rval, dst, src);
      assert( rval >= dsize);
      exit( -1);
      }
   return( rval);
}

/* 'Traditional' snprintf() silently truncates the output if it would
go past max_len bytes.  snprintf_err() is in all respects identical to
snprintf(),  except that truncation is considered to be an error and
the program will abort.    */

int snprintf_err( char *string, const size_t max_len,      /* miscell.cpp */
                                   const char *format, ...)
{
   va_list argptr;
   int rval;

   va_start( argptr, format);
#if defined(_MSC_VER) && _MSC_VER <= 1100
   rval = vsprintf( string, format, argptr);
#else
   rval = vsnprintf( string, max_len, format, argptr);
#endif
   va_end( argptr);
   if( (size_t)rval >= max_len || rval < 0)
      {
      fprintf( stderr, "snprintf_err : %d/%ld bytes\n", rval, (long)max_len);
      fprintf( stderr, "Format string '%s'\n", format);
      assert( 0);
      exit( -1);
      }
   return( rval);
}
