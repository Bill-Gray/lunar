#ifndef STRINGEX_H_INCLUDED
#define STRINGEX_H_INCLUDED

/* Some extensions to the usual string/snprintf functions.  These
silently truncate buffer overruns,  or abort so you know you
have such overruns.

First, the OpenBSD strlcxx functions.  See comments in 'stringex.cpp'. */

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */


#if !defined( __APPLE__) && !defined( __OpenBSD__)
size_t strlcpy( char *dst, const char *src, const size_t dsize);
size_t strlcat( char *dst, const char *src, const size_t dsize);
#endif

/* The above truncate silently.  On occasion,  one wants that.  Most of
the time,  such truncation indicates a bug,  and the program should
abort (at least in debug builds and,  arguably,  always).  The
'_err()' variants will do so. */

size_t strlcpy_err( char *dst, const char *src, const size_t dsize);
size_t strlcat_err( char *dst, const char *src, const size_t dsize);

/* Frequently,  'dsize' is simply sizeof( dst).  The following macros
can result in slightly easier to read code in such cases.      */

#define strlcat_error( a, b) strlcat_err( a, b, sizeof( a))
#define strlcpy_error( a, b) strlcpy_err( a, b, sizeof( a))

/* Older MSVC/C++ lacks snprintf().  See 'stringex.cpp' for details. */

#if defined(_MSC_VER) && _MSC_VER < 1900
   #define NO_BUILT_IN_SNPRINTF
#endif

#ifdef NO_BUILT_IN_SNPRINTF
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

int snprintf_append( char *string, const size_t max_len,
                                   const char *format, ...)
#ifdef __GNUC__
         __attribute__ (( format( printf, 3, 4)))
#endif
;

int snprintf_err( char *string, const size_t max_len,
                                   const char *format, ...)
#ifdef __GNUC__
         __attribute__ (( format( printf, 3, 4)))
#endif
;

#ifdef __cplusplus
}
#endif  /* #ifdef __cplusplus */
#endif  /* #ifndef STRINGEX_H_INCLUDED */
