#include <string.h>
#include "stringex.h"
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include "mpc_func.h"

static int64_t ten_to_the_nth( int n)
{
   int64_t rval = 1;

   while( n--)
      rval *= 10;
   return( rval);
}

int text_search_and_replace( char *str, const char *oldstr,
                                     const char *newstr)
{
   size_t ilen = strlen( str), rval = 0;
   const size_t oldlen = strlen( oldstr);
   const size_t newlen = strlen( newstr);

   while( ilen >= oldlen)
      if( !memcmp( str, oldstr, oldlen))
         {
         memmove( str + newlen, str + oldlen, ilen - oldlen + 1);
         memcpy( str, newstr, newlen);
         str += newlen;
         ilen -= oldlen;
         rval++;
         }
      else
         {
         str++;
         ilen--;
         }
   return( (int)rval);
}

/* See comments for get_ra_dec() in 'mpc_fmt.h' for info on   */
/* the meaning of 'precision' in this function.               */

void output_angle_to_buff( char *obuff, double angle, int precision)
{
   int n_digits_to_show = 0;
   int64_t power_mul, fraction;
   size_t i, full_len = 12;

   if( (precision >= 100 && precision <= 116) /* decimal quantity, dd.dd... */
        || ( precision >= 200 && precision <= 215)) /* decimal ddd.dd... */
      {
      const int two_digits = (precision <= 200);

      n_digits_to_show = precision % 100;
      power_mul = ten_to_the_nth( n_digits_to_show);
      fraction = (int64_t)( angle * (two_digits ? 1. : 15.) * (double)power_mul + .5);

      snprintf( obuff, 4, (two_digits ? "%02u" : "%03u"),
                  (int)( fraction / power_mul));
      fraction %= power_mul;
      }
   else
      switch( precision)
         {
         case -1:       /* hh mm,  integer minutes */
         case -2:       /* hh mm.m,  tenths of minutes */
         case -3:       /* hh mm.mm,  hundredths of minutes */
         case -4:       /* hh mm.mmm,  milliminutes */
         case -5:       /* hh mm.mmmm */
         case -6:       /* hh mm.mmmmm */
         case -7:       /* hh mm.mmmmmm */
            {
            n_digits_to_show = -1 - precision;
            power_mul = ten_to_the_nth( n_digits_to_show);
            fraction = (int64_t)( angle * 60. * (double)power_mul + .5);
            snprintf( obuff, 6, "%02u %02u", (unsigned)( fraction / ((int64_t)60 * power_mul)),
                                         (unsigned)( fraction / power_mul) % 60);
            fraction %= power_mul;
            }
            break;
         case 0:        /* hh mm ss,  integer seconds */
         case 1:        /* hh mm ss.s,  tenths of seconds */
         case 2:        /* hh mm ss.ss,  hundredths of seconds */
         case 3:        /* hh mm ss.sss,  thousands of seconds */
         case 4: case 5: case 6:    /* possible extra digits in ephems */
         case 7: case 8: case 9:
         case 307:      /* hhmmsss,  all packed together:  tenths */
         case 308:      /* hhmmssss,  all packed together: hundredths */
         case 309:      /* milliseconds (or milliarcseconds) */
         case 310:      /* formats 307-312 are the 'super-precise' formats */
         case 311:
         case 312:      /* microseconds (or microarcseconds) */
         case 400:      /* (RA) ddd mm ss,  integer arcseconds */
         case 401:      /* (RA) ddd mm ss.s,  tenths of arcseconds */
         case 402:      /* (RA) ddd mm ss.ss,  centiarcsec */
            {
            const char *format = "%02u %02u %02u";

            if( precision >= 400)
               {
               format = "%03u %02u %02u";
               precision -= 400;
               angle *= 15.;
               }
            n_digits_to_show = precision % 306;
            power_mul = ten_to_the_nth( n_digits_to_show);
            fraction = (int64_t)( angle * 3600. * (double)power_mul + .5);
            snprintf( obuff, 10, format,
                     (unsigned)( fraction / ((int64_t)3600 * power_mul)),
                     (unsigned)( fraction / ((int64_t)60 * power_mul)) % 60,
                     (unsigned)( fraction / power_mul) % 60);
            fraction %= power_mul;
            if( precision > 306)          /* remove spaces: */
               text_search_and_replace( obuff, " ", "");
            }
            break;
         default:                  /* try to show the angle,  but indicate */
            if( angle > -1000. && angle < 1000.)   /* the format is weird  */
               snprintf( obuff, 10, "?%.5f", angle);
            else
               strlcpy_error( obuff, "?");
            fraction = 0;   /* not really necessary;  evades nuisance GCC warning */
            break;
         }
                  /* Formats not used in astrometry -- they don't fit the */
                  /* field size for punched-card data,  but are used in ephems */
   if( precision >= 4 && precision <= 11)    /* 'overlong' dd mm ss.ssss... */
      full_len += (size_t)( precision - 3);
   if( precision >= 209 && precision <= 215) /* 'overlong' ddd.ddddd... */
      full_len += (size_t)( precision - 208);
   if( precision >= 110 && precision <= 116) /* 'overlong' dd.dddddd... */
      full_len += (size_t)( precision - 109);
   if( n_digits_to_show)
      {
      char format[8];

      if( precision < 307 || precision > 312)   /* omit decimal point for */
         strcat( obuff, ".");                    /* super-precise formats */
      assert( n_digits_to_show > 0 && n_digits_to_show < 20);
      snprintf_err( format, sizeof( format), "%%0%dld", n_digits_to_show);
      snprintf_append( obuff, full_len + 1, format, (long)fraction);
      }
   for( i = strlen( obuff); i < full_len; i++)
      obuff[i] = ' ';
   obuff[full_len] = '\0';
}

void output_signed_angle_to_buff( char *obuff, const double angle,
                               const int precision)
{
   *obuff++ = (angle < 0. ? '-' : '+');
   output_angle_to_buff( obuff, fabs( angle), precision);
}
