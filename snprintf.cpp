#include <stdio.h>
#include <stdlib.h>

#if defined(_MSC_VER) && _MSC_VER < 1900
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
#include <stdarg.h>

int snprintf( char *string, const size_t max_len, const char *format, ...)
{
   va_list argptr;
   int rval;

   va_start( argptr, format);
#if _MSC_VER <= 1100
   rval = vsprintf( string, format, argptr);
#else
   rval = vsnprintf( string, max_len, format, argptr);
#endif
   string[max_len - 1] = '\0';
   va_end( argptr);
   return( rval);
}
#endif    /* #ifdef _MSC_VER */
