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

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"
#include "mpc_func.h"
#include "stringex.h"

#if defined(_MSC_VER) && _MSC_VER < 1900
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

/* MPC has,  at least thus far,  only assigned MPC codes that are an uppercase
letter followed by two digits.  Look at 'rovers.txt',  and you'll see that
user-added "MPC codes" are not so limited;  e.g.,  "2.2" and "0.6" are codes
used for the 2.2-m and 0.6-m telescopes at Mauna Kea.  */

bool is_valid_mpc_code( const char *mpc_code)
{
   int i;

   for( i = 0; i < 3; i++)
      if( mpc_code[i] <= ' ' || mpc_code[i] > 'z')
         return( false);
   return( true);
}

double quick_atof( const char *ibuff);       /* mpc_func.cpp */

static inline int quick_atoi( const char *iptr)
{
   int rval = 0;

   while( *iptr >= '0' && *iptr <= '9')
      rval = (rval * 10) + (*iptr++ - '0');
   return( rval);
}

static inline int get_two_digits( const char *iptr)
{
   return( (int)iptr[0] * 10 + (int)iptr[1] - (int)'0' * 11);
}

/* The date/time in an 80-column MPC report line is stored in columns 16-32.
MPC expects this to be in the form YYYY MM DD.dddddd,  with the date usually
given to five digits = .864 second precision;  the sixth digit should only
be given if you're really sure the time is good to that level.  In some
cases (mostly older,  poorly timed observations),  fewer digits are given.
A few ancient observations are just given to the nearest day.

Find_Orb also supports some NON-STANDARD,  FIND_ORB ONLY (i.e.,  MPC will
reject anything you do in these formats).  These formats support greater
precision and spare you the need to convert dates from HH:MM:SS.sss format
to decimal-day format,  and/or the need to convert from JD or MJD format
to YYYY MM SS format.  To demonstrate,  the following lines give the same
date/time in the various formats,  give or take a little rounding error :

            2013 02 13.141593     (MPC's expected format)
            2456336.641592653     (Julian Day format)
            M056336.141592653     (MJD = Modified Julian Day)
            K130213.141592653     (CYYMMDD.dddddd)
            K130213:032353605     (CYYMMDD HH:MM:SS.sss)

The last two use the MPC "mutant hex" convention of K for 21st century
years and J for twentieth-century years,  followed by a two-digit year.
After the two-digit month and two-digit day of month,  one can have
a decimal point and decimal day _or_ a colon and HHMMSS,  and optionally,
up to millisecond precision.  (Though you can and should supply fewer
digits if -- as will almost always be the case -- your timing is not
really that exact.)  The last one is a little confusing,  given the
placement of the ':':  it actually means "2013 02 13 03:23:53.605".

Note that these non-standard formats allow precision up to 10^-9 day
(86.4 microseconds) or,  for the last format,  one millisecond.  This
precision should be good enough for anybody.  Indeed,  that last digit
is just past the precision limit of a 64-bit float.  (We could get
around this by having dates relative to J2000 instead of using MJD
dates,  and/or by using 80-bit "long doubles".  See my get_time.cpp
code in the 'lunar' library for details.  But I don't see it as a big
issue,  at least not yet.)

For each of these formats,  a very little bit of format checking is
done (make sure digits are in certain key places,  and that the full
line is exactly 80 bytes).  Some malformed records _can_ slip past!

Note also that MPC-formatted times are on the UTC scale,  but don't
provide a way to express observations made on the sheer lunacy that is
a leap second.  The HH:MM:SS formats below _should_ let you do this,
but (as yet) do not.  (The other formats flat out can't express that
extra 86401st second in a day.)

      53    K130213:0323.870      (CYYMMDD HH:MM.mmm... not supported)
      52    K130213:0323.87       (CYYMMDD HH:MM.mm... not supported)
      51    K130213:0323.8        (CYYMMDD HH:MM.m... not supported)
      50    K130213:0323          (CYYMMDD HH:MM)

      49    M056336.641592653     (MJD, 10^-9 day)
      48    M056336.64159265      (MJD, 10^-8 day)
      47    M056336.6415926       (MJD, 10^-7 day)
      46    M056336.641592        (MJD, 10^-6 day)
      45    M056336.64159         (MJD, 10^-5 day)
      44    M056336.6415          (MJD, 10^-4 day)
      43    M056336.641           (MJD,  .001 day)
      42    M056336.64            (MJD,   .01 day)
      41    M056336.6             (MJD,    .1 day... not supported)
      40    M056336.              (MJD,     1 day... not supported)

      39    K130213.141592653     (High-prec, 10^-9 day)
      38    K130213.14159265      (High-prec, 10^-8 day)
      37    K130213.1415926       (High-prec, 10^-7 day)
      36    K130213.141592        (High-prec, 10^-6 day)
      35    K130213.14159         (High-prec, 10^-5 day)
      34    K130213.1415          (High-prec, 10^-4 day)
      33    K130213.141           (High-prec,  .001 day)
      32    K130213.14            (High-prec,   .01 day)
      31    K130213.1             (High-prec,    .1 day... not supported)
      30    K130213.              (High-prec,     1 day... not supported)

      23    K130213:032353605     (CYYMMDD HH:MM:SS.sss)
      22    K130213:03235360      (CYYMMDD HH:MM:SS.ss)
      21    K130213:0323536       (CYYMMDD HH:MM:SS.s)
      20    K130213:032353        (CYYMMDD HH:MM:SS)

      19    2456336.641592653     (Julian Day, 10^-9 day)
      18    2456336.64159265      (Julian Day, 10^-8 day)
      17    2456336.6415926       (Julian Day, 10^-7 day)
      16    2456336.641592        (Julian Day, 10^-6 day)
      15    2456336.64159         (Julian Day, 10^-5 day)
      14    2456336.6415          (Julian Day, 10^-4 day)
      13    2456336.641           (Julian Day, 10^-3 day)
      12    2456336.64            (Julian Day, 10^-2 day... not supported)
      11    2456336.6             (Julian Day, 10^-1 day... not supported)

       6    2013 02 13.141593     (MPC's expected format, 10^-6 day)
       5    2013 02 13.14159      (MPC's expected format, 10^-5 day)
       4    2013 02 13.1415       (MPC's expected format, 10^-4 day)
       3    2013 02 13.141        (MPC's expected format, 10^-3 day)
       2    2013 02 13.14         (MPC's expected format, 10^-2 day)
       1    2013 02 13.1          (MPC's expected format, 10^-1 day)
       0    2013 02 13.           (MPC's expected format, 10^-0 day) */

double extract_date_from_mpc_report( const char *buff, unsigned *format)
{
   double rval = 0.;
   int year = 0, month = 0;
   size_t start_of_decimals = 0;
   unsigned format_found = 0;
   char tbuff[18];
   const size_t len = strlen( buff);
   unsigned i, bit, digits_mask = 0;

   if( len < 80 || len > 82)       /* check for correct length */
      return( 0.);
   if( buff[12] != ' ' && buff[12] != '*' && buff[12] != '-')
      return( 0.);
   if( !is_valid_mpc_code( buff + 77))
      return( 0.);
   memcpy( tbuff, buff + 15, 17);
   for( i = 0, bit = 1; i < 17; i++, bit <<= 1)
      if( isdigit( tbuff[i]))
         digits_mask |= bit;
   tbuff[17] = '\0';
   if( (digits_mask & 0x3ff) != 0x36f    /* i.e.,  'dddd zd zd' */
                   && (digits_mask & 0x2df) == 0x24f
                   && tbuff[7] == ' ' && tbuff[10] == '.')
      {                        /* no leading zero given for month,  day, */
      if( tbuff[5] == ' ')     /* or both;  fill it in with zero(es)  */
         tbuff[5] = '0';
      if( tbuff[8] == ' ')
         tbuff[8] = '0';
      digits_mask |= 32 + 256;
      }
   if( tbuff[4] == ' ')
      {                       /* standard YYYY MM DD.dddddd format */
      if( (digits_mask & 0x3ff) == 0x36f    /* i.e.,  'dddd dd dd' */
                            && tbuff[7] == ' ' && tbuff[10] == '.')
         {
         int divisor = 1000000;

         year = quick_atoi( tbuff);
         month = get_two_digits( tbuff + 5);
//       rval = atof( tbuff + 8);
                     /* atof( ) is a little slow,  so we use a little more */
         for( i = 16; i > 11 && tbuff[i] == ' '; i--)  /* code in exchange */
            divisor /= 10;                             /* for better speed */
         rval = (double)get_two_digits( tbuff + 8) +
                           (double)quick_atoi( tbuff + 11) / (double)divisor;
         format_found = 0;
         start_of_decimals = 11;
         }
      }
   else if( *tbuff >= 'H' && *tbuff <= 'K')  /* 18th through 21st century */
      {                                          /* CYYMMDD format */
      if( (tbuff[7] == '.' || tbuff[7] == ':')
               && (digits_mask & 0x3ff) == 0x37e)  /* i.e, 'Zdddddd.dd' */
         {
         year = (*tbuff - 'J') * 100 + 1900 +
                    get_two_digits( tbuff + 1);
         month = get_two_digits( tbuff + 3);
         rval = atof( tbuff + 5);
         if( tbuff[7] == ':')
            {
            rval += (double)get_two_digits( tbuff + 8) / hours_per_day
               + (double)get_two_digits( tbuff + 10) / minutes_per_day;
            if( tbuff[12] != ' ')
               {
               rval += (double)get_two_digits( tbuff + 12) / seconds_per_day;
               tbuff[13] = '.';
               rval += atof( tbuff + 13) / seconds_per_day;
               format_found = 20;      /* formats 20-23;  see above */
               start_of_decimals = 14;
               }
            else
               format_found = 50;      /* HH:MM format;  see above */
            }
         else     /* decimal formats 32-40 */
            {
            format_found = 30;
            start_of_decimals = 8;
            }
         }
      }
   else if( tbuff[7] == '.')        /* MJD or JD format */
      {
      if( (digits_mask & 0x3fe) == 0x37e)   /* i.e., 'zdddddd.dd' */
         {
         if( *tbuff == 'M')    /* MJD */
            {
            format_found = 40;
            rval = 2400000.5 + atof( tbuff + 1);
            }
         else
            {
            format_found = 10;
            rval = atof( tbuff);      /* plain ol' JD */
            }
         start_of_decimals = 8;
         }
      }
   if( format)
      {
      if( start_of_decimals)
         while( isdigit( tbuff[start_of_decimals++]))
            format_found++;
      *format = format_found;
      }

   if( month >= 1 && month <= 12 && rval > 0. && rval < 99.)
      rval += (double)dmy_to_day( 0, month, year,
                                    CALENDAR_JULIAN_GREGORIAN) - .5;

             /* Radar obs are always given to the nearest UTC second. So  */
             /* some rounding is usually required with MPC microday data. */
   if( rval && (buff[14] == 'R' || buff[14] == 'r'))
      {
      const double time_of_day = rval - floor( rval);
      const double resolution = 1. / seconds_per_day;
      const double half = .5 / seconds_per_day;

      rval += half - fmod( time_of_day + half, resolution);
      }
   return( rval);
}

/* get_ra_dec() looks at an RA or dec from an MPC report and returns
its precision.  It interprets the formats used by MPC,  plus a lot of
"extended" formats that can be useful if your input data is in other
formats and/or has extra digits.  The return value has the following
meanings ('z' = 'hours or degrees',  i.e.,  this format can apply to
both RAs and decs.)
   hh mm ss.sss       3    (MPC uses this for 'precise' RAs)
   zz mm ss.ss        2    (MPC uses for most RAs & 'precise' decs)
   zz mm ss.s         1    (MPC uses for most decs & low-precision RA)
   zz mm ss           0    (Used by MPC, only rarely)
   zz mm             -1    (Used _very_ rarely by MPC)
   zz mm.m           -2    (Maybe used by MPC once or twice)
   zz mm.mm          -3       The following are for Find_Orb only
   zz mm.mmm         -4
   zz mm.mmmm        -5
   zz mm.mmmmm       -6
   zz mm.mmmmmm      -7
   zz.               100
   zz.z              101
   zz.zz             102
   zz.zzz            103
   zz.zzzz           104... can go up to nine places = 109 in RA,
                           or to 108 = eight places in dec
   ddd.              200  (used for RA only)
   ddd.d             201
   ddd.dd            202... can go up to eight places = 208
   HHMMSSs           307 (RA to .1 second)
   ddmmSSs           307 (dec to .1 arcsecond)
   HHMMSSss          308 (RA to .01 second)
   ddmmSSss          308 (dec to .01 arcsecond)
           ... and so forth until...
   ddmmSSsssss       311 (dec to 10 microarcseconds)
   HHMMSSssssss      312 (RA to one microsecond)

   ddd mm ss.        400   (used for RA only)
   ddd mm ss.s       401   (used for RA only)
   ddd mm ss.ss      402   (used for RA only)
   Undetermined      -99

   Please note that (at least thus far) I've only seen the first six
cases used by MPC,  and they will probably balk at anything sent in
anything but the first four formats.

   The remaining formats have been quite useful in situations where I
find data in a non-MPC format;  I can leave it in decimal degrees or
whatnot,  and can accommodate extra digits for super-high-accuracy
data.  That's why formats 307-312 were added;  they accommodate some
highly precise VLBA astrometry that really _is_ good to the tens of
microarcseconds level.  See

http://iau-comm4.jpl.nasa.gov/plan-eph-data/vlbaobs.html

   And Gaia,  at least in theory,  will be of a similar level of
accuracy;  we need to Be Prepared for that.
   The precision is stored and used to "recreate" the RA/dec in the
original form.  (It could,  and probably should,  also be used to
weight the observations.)                                             */

#define BAD_RA_DEC_FMT           -99

static double get_ra_dec( const char *ibuff, int *format, double *precision)
{
   char buff[13];
   double rval;
   unsigned i = 0, n_digits = 0;
   const bool is_dec = (*ibuff == '-' || *ibuff == '+');
   const bool is_negative = (*ibuff == '-');

   *precision = 1.;   /* in arcseconds */
   if( is_dec)
      ibuff++;
   memcpy( buff, ibuff, 12);
   buff[12] = '\0';
   rval = quick_atof( buff);
   while( isdigit( buff[i]))
      i++;
   if( i > 7)        /* "packed" highly precise RA/decs described above */
      {
      unsigned tval;
      double factor = 1. / 3600.;

      *format = 300 + i;
      n_digits = i - 6;
      buff[6] = '\0';
      tval = atoi( buff);
      rval = (double)( tval / 10000)
           + (double)( (tval / 100) % 100) / 60.
           + (double)( tval % 100) / 3600.;
      for( i = 6; i < 12 && isdigit( ibuff[i]); i++)
         {
         factor *= .1;
         rval += (double)( ibuff[i] - '0') * factor;
         }
      buff[6] = ibuff[6];
//    debug_printf( "Extended: '%s' = %.8f\n", ibuff, rval);
      }
   else if( buff[2] == '.')        /* decimal degrees or hours */
      {
      *precision = 3600.;
      for( i = 3; isdigit( buff[i]) && i < 12; i++)
         n_digits++;
      *format = 100 + n_digits;
      }
   else if( buff[3] == '.')        /* decimal degrees for RA,  ddd.ddd.... */
      {
      *precision = 3600. / 15.;
      rval /= 15.;
      for( i = 4; isdigit( buff[i]) && i < 12; i++)
         n_digits++;
      *format = 200 + n_digits;
      }
   else if( buff[2] == ' ' && i == 2)  /* zz mm(.mmm...) or zz mm ss(.sss...)  */
      {
      rval += quick_atof( buff + 3) / 60.;
      if( buff[5] == ' ' && isdigit( buff[7]))  /* i.e., seconds are given */
         {
         rval += quick_atof( buff + 6) / 3600.;
         for( i = 9; i < 12 && isdigit( buff[i]); i++)
            n_digits++;
         *format = n_digits;
         }
      else           /* minutes: */
         {
         for( i = 6; i < 12 && isdigit( buff[i]); i++)
            n_digits++;
         *format = -((int)n_digits + 1);
         *precision = 60.;
         }
      }
   else if( buff[3] == ' ' && i == 3 && !is_dec)    /* ddd mm ss(.s) RA */
      {
      *precision = 1. / 15.;
      *format = 2;
      rval += quick_atof( buff + 4) / 60. + quick_atof( buff + 7) / 3600.;
      rval /= 15;
      for( i = 10; isdigit( buff[i]) && i < 12; i++)
         n_digits++;
      *format = 400 + n_digits;
      }
   else
      *format = BAD_RA_DEC_FMT;
   while( n_digits--)
      *precision *= .1;
   if( is_negative)
      rval = -rval;
   return( rval);
}

/* Extracts both the RA and the dec from an MPC-formatted report,  subject
to all the weirdnesses described above.  Return value is zero if both
RA and dec are successfully found;  -1 if the RA is bad;  -2 if the dec
is bad;  -3 if neither was read. */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int get_ra_dec_from_mpc_report( const char *ibuff,
                       int *ra_format, double *ra, double *ra_precision,
                       int *dec_format, double *dec, double *dec_precision)
{
   int rval = 0, format;
   double prec;

   *ra  = get_ra_dec( ibuff + 32, &format, &prec) * (PI / 12.);
   if( ra_precision)
      *ra_precision = prec * 15.;     /* cvt minutes/seconds to arcmin/arcsec */
   if( format == BAD_RA_DEC_FMT)
      rval = -1;
   if( ra_format)
      *ra_format = format;

   *dec =  get_ra_dec( ibuff + 44, &format, &prec) * (PI / 180.);
   if( dec_precision)
      *dec_precision = prec;
   if( format == BAD_RA_DEC_FMT)
      rval -= 2;
   if( dec_format)
      *dec_format = format;
   return( rval);
}

static const char *net_codes[] = {
    /* http://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
         G. V. Williams, 2012, ``Minor Planet Astrophotometry'', PhD
         thesis, Open University. [2012PhDT.........7W] */
           "aUSNO-A1",
           "bUSNO-SA1",
           "cUSNO-A2",
           "dUSNO-SA2",
           "eUCAC-1",
           "fTycho-1",
           "gTycho-2",
           "hGSC-1.0",
           "iGSC-1.1",
           "jGSC-1.2",
           "kGSC-2.2",
           "lACT",
           "mGSC-ACT",
           "nSDSS-DR8",       /* was TRC */
           "oUSNO-B1",
           "pPPM",
           "qUCAC-4",     /* also UCAC2-beta for earlier obs */
           "rUCAC-2",
           "sUSNO-B2",
           "tPPMXL",
           "uUCAC-3",
           "vNOMAD",
           "wCMC-14",
           "xHIP-2",
           "yHIP-1",
           "zGSC-?",        /* no version specified */
           "AAC",
           "BSAO 1984",
           "CSAO",
           "DAGK 3",
           "EFK4",
           "FACRS",
           "GLick Gaspra Catalogue",
           "HIda93 Catalogue",
           "IPerth 70",
           "JCOSMOS/UKST Southern Sky Catalogue",
           "KYale",
           "L2MASS", /* used for WISE & PanSTARRS astrometry */
           "MGSC-2.3",
           "NSDSS-DR7",
           "OSST-RC1",
           "PMPOSC3",
           "QCMC-15",
           "RSST-RC4",
           "SURAT-1",
           "TURAT-2",
           "UGaia-DR1",
           "VGaia-DR2",
           "WGaia-DR3",
           "XGaia-EDR3",
           "XGaia-3E",
           "YUCAC-5",
           "ZATLAS-2",
           NULL };

const char *byte_code_to_net_name( const char byte_code)
{
   size_t i;

   for( i = 0; net_codes[i]; i++)
      if( byte_code == net_codes[i][0])
         return( net_codes[i] + 1);
   return( NULL);
}

/* Code to get around the fact that people (probably shouldn't,  but do)
specify NETs in a variety of nonstandard ways,  such as

NET Gaia DR1.0
NET Gaia DR1
NET Gaia-DR1
NET Gaiadr1
NET Gaia1

   If the names match after ignoring '.0', '-',  and spaces and
upper/lower case,  and dropping the 'DR' for Gaia,  we've almost
assuredly got the right catalog.   */

static void reduce_net_name( char *obuff, const char *ibuff)
{
   while( *ibuff)
      if( *ibuff == '-' || *ibuff == ' ')
         ibuff++;
      else if( ibuff[0] == '.' && ibuff[1] == '0')
         ibuff += 2;
      else if( ibuff[0] == 'D' && ibuff[1] == 'R')
         ibuff += 2;
      else
         *obuff++ = toupper( *ibuff++);
   *obuff = '\0';
}

char net_name_to_byte_code( const char *net_name)
{
   char net1[80], net2[80], rval = 0;
   size_t i;

   reduce_net_name( net1, net_name);
   for( i = 0; net_codes[i]; i++)
      {
      reduce_net_name( net2, net_codes[i] + 1);
      if( !strcmp( net1, net2))
         rval = net_codes[i][0];
      }
   if( !rval)     /* didn't find it */
      rval = '?';
   return( rval);
}

/* Currently used to determine that,  e.g., Jupiter XLVIII = Jupiter 48. */

static int extract_roman_numeral( const char *str)
{
   int rval = 0;

   while( *str)
      {
      switch( *str)
         {
         case 'I':
            rval += (str[1] == 'V' || str[1] == 'X' ? -1 : 1);
            break;
         case 'V':
            rval += 5;
            break;
         case 'X':
            rval += (str[1] == 'L' || str[1] == 'C' ? -10 : 10);
            break;
         case 'L':
            rval += 50;
            break;
         case 'C':
            rval += (str[1] == 'D' || str[1] == 'M' ? -100 : 100);
            break;
         case 'D':
            rval += 500;
            break;
         case 'M':
            rval += 1000;
            break;
         default:
            return( -1);
            break;
         }
      str++;
      }
   return( rval);
}

/* Packs desigs such as 'Uranus XXIV' or 'Earth I' or 'Saturn DCCXLVIII',
or 'Uranus 24' or 'Earth 1' or 'Saturn 748'.  Packed designations are,
including padding,  exactly twelve bytes long plus a null terminator. */

#define PACKED_DESIG_LEN   13

static int pack_permanent_natsat( char *packed, const char *fullname)
{
   size_t i;

   for( i = 0; i < 8; i++)
      {
      extern const char *planet_names_in_english[];      /* unpack.cpp */
      const size_t len = strlen( planet_names_in_english[i]);

      if( !strncmp( fullname, planet_names_in_english[i], len)
               && fullname[len] == ' ')
         {
         int sat_number = atoi( fullname + len + 1);

         if( !sat_number)
            sat_number = extract_roman_numeral( fullname + len + 1);
         if( sat_number > 0 && sat_number < 1000)
            {
            snprintf_err( packed, PACKED_DESIG_LEN, "%c%03dS       ", *fullname, sat_number);
            return( (int)i);
            }
         }
      }
   return( -1);
}

/* Packs desigs such as 'S/1999 J 1' or 'S/2013 N 42'.   */

static int pack_provisional_natsat( char *packed, const char *fullname)
{
   if( *fullname == 'S' && fullname[1] == '/' && strlen( fullname) > 9
            && fullname[6] == ' ' && fullname[8] == ' '
            && strchr( "MVEMJSUNP", fullname[7]))
      {
      const int year = atoi( fullname + 2);
      const int num = atoi( fullname + 9);

      if( year > 1900 && year < 2100 && num > 0 && num < 620)
         {
         snprintf_err( packed, PACKED_DESIG_LEN, "    S%c%02d%c%c%c0",
               'A' + (year / 100) - 10, year % 100, fullname[7],
                     int_to_mutant_hex_char( num / 10),
                     '0' + num % 10);
         return( 0);
         }
      }
   return( -1);
}

/* Returns either a negative value for an error code,  or the
location of the decimal point for a valid coordinate */

inline int get_satellite_coordinate( const char *iptr, double coord[1])
{
   char tbuff[12];
   const char sign_byte = *iptr;
   int rval = 0;

   memcpy( tbuff, iptr, 11);
   tbuff[11] = '\0';
   if( sign_byte != '+' && sign_byte != '-')
      {
      rval = SATELL_COORD_ERR_BAD_SIGN;
      *coord = atof( tbuff);
      }
   else
      {
      char *tptr;
      int n_bytes_read;

      if( sscanf( tbuff + 1, "%lf%n", coord, &n_bytes_read) != 1
                     || n_bytes_read < 7)
         rval = SATELL_COORD_ERR_BAD_NUMBER;
      else if( (tptr = strchr( tbuff, '.')) == NULL)
         rval = SATELL_COORD_ERR_NO_DECIMAL;
      else
         rval = (int)( tptr - tbuff);
      if( sign_byte == '-')
         *coord = -*coord;
      }
   return( rval);
}

/* Returns spacecraft offset in AU,  in J2000 ecliptic coords */

int DLL_FUNC get_satellite_offset( const char *iline, double xyz[3])
{
   size_t i;
   int error_code = 0;
   const int observation_units = (int)iline[32] - '0';
   double r2 = 0.;
   const double earth_radius_in_au = 6378.14 / AU_IN_KM;
   const double min_radius = 1.01 * earth_radius_in_au;
   const size_t slen = strlen( iline);

   assert( 80 <= slen && slen < 83);  /* allow for LF, CR/LF, or no line end */
   for( i = 0; i < 3; i++)    /* in case of error,  use 0 offsets */
      xyz[i] = 0.;
   iline += 34;      /* this is where the offsets start */
   for( i = 0; i < 3; i++, iline += 12)
      {
      int decimal_loc;

      decimal_loc = get_satellite_coordinate( iline, xyz + i);
      if( decimal_loc < 0 && !error_code)
         error_code = decimal_loc;
      if( observation_units == 1)         /* offset given in km */
         {
         xyz[i] /= AU_IN_KM;
         if( !error_code)
            if( decimal_loc < 6 || decimal_loc > 8)
               error_code = SATELL_COORD_ERR_DECIMAL_MISPLACED;
         if( !error_code && xyz[i] == 0.)
            error_code = SATELL_COORD_ERR_EXACTLY_ZERO;
         }
      else if( observation_units == 2)          /* offset in AU */
         {                         /* offset must be less than 100 AU */
         if( !error_code)
            if( decimal_loc != 2 && decimal_loc != 3)
               error_code = SATELL_COORD_ERR_DECIMAL_MISPLACED;
         }
      else if( !error_code)      /* don't know about this sort of offset */
         error_code = SATELL_COORD_ERR_UNKNOWN_OFFSET;
      r2 += xyz[i] * xyz[i];
      }
               /* (275) geocentric occultation obs can have offsets  */
               /* inside the earth.  All others ought to be at least */
               /* slightly outside the atmosphere.                   */
   if( !error_code && r2 < min_radius * min_radius)
      if( !strcmp( iline + 77, "275"))
         error_code = SATELL_COORD_ERR_INSIDE_EARTH;
   equatorial_to_ecliptic( xyz);
   return( error_code);
}

/* create_mpc_packed_desig( ) takes a "normal" name for a comet/asteroid,
such as P/1999 Q1a or 2005 FF351,  and turns it into the 12-byte packed
format used in MPC reports and element files.  ('packed_desig' will
_always_ be a 12-byte-long,  null-terminated string.)  Documentation of
the (somewhat cryptic) packed format is given on the MPC Web site :

https://minorplanetcenter.net//iau/info/PackedDes.html

   A 'test' main( ) at the end of this file shows the usage of this function.

   This should handle all "normal" asteroid and comet designations.
It doesn't handle natural satellites yet.  If the designation can't be
turned into a valid packed designation,  you get a dollar sign ($)
followed by the first eleven bytes of 'obj_name'.  Find_Orb,  at least,
will recognize a "packed designation" starting with $ as meaning "a
literal, unpacked name follows".

   The maximum number that can be packed under the current MPC scheme
is 620000 + 62^4.  */

int create_mpc_packed_desig( char *packed_desig, const char *obj_name)
{
   size_t i, j, len;
   int rval = -1;
   bool in_parentheses = false;
   unsigned number = 0;
   char comet_desig = 0;
   const unsigned max_number = 620000 + 62 * 62 * 62 * 62;

   while( *obj_name == ' ')
      obj_name++;
   if( *obj_name == '(')   /* possible numbered desig such as (433),  etc. */
      {
      j = 1;
      while( isdigit( obj_name[j]))
         j++;
      if( obj_name[j] == ')' && !obj_name[j + 1])
         {
         obj_name++;
         in_parentheses = true;
         }
      }
               /* Check for comet-style desigs such as 'P/1995 O1' */
               /* and such.  Leading character can be P, C, X, D, or A. */
   if( obj_name[1] == '/' && strchr( "PCXDA", *obj_name))
      {
      comet_desig = *obj_name;
      obj_name += 2;
      }

   memset( packed_desig, ' ', 12);
   packed_desig[12] = '\0';
   i = 0;
   while( isdigit( obj_name[i]))
      number = (number * 10) + (unsigned)( obj_name[i++] - '0');
   len = strlen( obj_name);
   if( in_parentheses)
      len--;
   while( len && obj_name[len - 1] == ' ')
      len--;         /* ignore trailing spaces */
   if( i == len)        /* nothing but numbers */
      if( number > 0 && number < 1000000 && i >= 5)
         in_parentheses = true;
   if( number > 0 && number < 10000 && obj_name[i]
                      && (!obj_name[i + 1] || obj_name[i + 1] == '-')
                             && strchr( "PDI", obj_name[i]))
      {        /* such as '297P',  '1I',  etc. */
      snprintf( packed_desig, 5, "%04d", number);
      packed_desig[4] = obj_name[i];
      if( obj_name[i + 1] == '-' && isupper( obj_name[i + 2]))
         {           /* fragment designation */
         if( isupper( obj_name[i + 3]))  /* two-letter fragment */
            {
            packed_desig[10] = obj_name[i + 2] + 'a' - 'A';
            packed_desig[11] = obj_name[i + 3] + 'a' - 'A';
            }
         else        /* single-letter fragment desig */
            packed_desig[11] = obj_name[i + 2] + 'a' - 'A';
         }
      return( 0);
      }

   if( i < len && obj_name[i] == ' ')
      i++;
               /* If the name starts with four digits followed by an */
               /* uppercase letter,  it's a provisional designation: */
   if( number > 999 && number < 9000 && isupper( obj_name[i]))
      {
      int sub_designator;
      bool mangled_designation = false;

      for( j = 0; j < 4; j++)
         {
         const char *surveys[4] = { "P-L", "T-1", "T-2", "T-3" };

         if( !strcmp( obj_name + i, surveys[j]))
            {
            const char *surveys_packed[4] = {
                     "PLS", "T1S", "T2S", "T3S" };

            memcpy( packed_desig + 8, obj_name, 4);
            memcpy( packed_desig + 5, surveys_packed[j], 3);
            return( 0);
            }
         }

      if( number < 6200)
         packed_desig[5] = int_to_mutant_hex_char( number / 100);
      packed_desig[6] = obj_name[2];    /* decade */
      packed_desig[7] = obj_name[3];    /* year */

      packed_desig[8] = (char)toupper( obj_name[i]);    /* prelim desigs */
      i++;                            /* are _very_ scrambled when packed */
      if( isupper( obj_name[i]))
         {
         packed_desig[11] = obj_name[i];
         i++;
         }
      else if( !comet_desig)  /* asteroid desigs _must_ have a second */
         mangled_designation = true;              /* uppercase letter */
      else
         packed_desig[11] = '0';

      sub_designator = quick_atoi( obj_name + i);
      if( sub_designator >= 0 && number < 6200)
         {
         if( sub_designator < 620)
            {
            packed_desig[10] = int_to_mutant_hex_char( sub_designator % 10);
            packed_desig[9] = int_to_mutant_hex_char( sub_designator / 10);
            }
         else if( number >= 2000 && number < 2062)
            {
            int n = (sub_designator - 620) * 25 + packed_desig[11] - 'A';

            if( packed_desig[11] > 'I')
               n--;
            if( n >= 62 * 62 * 62 * 62)
               mangled_designation = true;
            else
               {
               packed_desig[5] = '_';
               packed_desig[6] = int_to_mutant_hex_char( number - 2000);
               packed_desig[7] = packed_desig[8];    /* move half-month specifier */
               encode_value_in_mutant_hex( packed_desig + 8, 4, n);
               }
            }
         while( isdigit( obj_name[i]))
            i++;
         if( comet_desig)
            {
            packed_desig[4] = comet_desig;
            if( obj_name[i] == '-' && isupper( obj_name[i + 1]))
               {        /* Comet fragment such as C/2018 F4-A */
               i++;     /* pack as,  e.g.,  CK18F04a */
               packed_desig[11] = obj_name[i++] + 'a' - 'A';
               }
            }
         if( i == len && !mangled_designation)  /* successfully unpacked desig */
            rval = 0;
         }
      }
   else if( in_parentheses && i == len && number < max_number && number > 0
                                       && !comet_desig)
      {           /* permanently numbered asteroid */
      rval = 0;
      if( number < 620000)
         snprintf_err( packed_desig, PACKED_DESIG_LEN, "%c%04d       ",
               int_to_mutant_hex_char( number / 10000),
               number % 10000);
      else
         {
         packed_desig[0] = '~';
         encode_value_in_mutant_hex( packed_desig + 1, 4, number - 620000);
         }
      }
   else if( number < 10000 && number > 0 && comet_desig)
      {           /* permanently numbered asteroid */
      int n_fragment_letters = 0;

      if( obj_name[i] == '-' && isupper( obj_name[i + 1]))
         {
         if( i + 2 == len)
            n_fragment_letters = 1;
         else if( i + 3 == len && isupper( obj_name[i + 2]))
            n_fragment_letters = 2;
         }
      if( i == len || n_fragment_letters)
         {
         rval = 0;
         snprintf_err( packed_desig, PACKED_DESIG_LEN, "%04d%c       ", number, comet_desig);
         if( n_fragment_letters == 1)
            packed_desig[11] = obj_name[i + 1] + 'a' - 'A';
         if( n_fragment_letters == 2)
            {
            packed_desig[10] = obj_name[i + 1] + 'a' - 'A';
            packed_desig[11] = obj_name[i + 2] + 'a' - 'A';
            }
         }
      }
   else if( number >= 1957 && number < 2100 && obj_name[4] == '-' && len >= 9
                  && len <= 11)
      {
      if( isdigit( obj_name[5]) && isdigit( obj_name[6]) && isdigit( obj_name[7])
                  && isupper( obj_name[8]))
         {
         rval = 0;     /* artsat desig such as '1992-044A' or '2000-357KHA' */
         memcpy( packed_desig, obj_name, len);
         }
      }
   else if( pack_permanent_natsat( packed_desig, obj_name) >= 0)
      rval = 0;
   else if( pack_provisional_natsat( packed_desig, obj_name) >= 0)
      rval = 0;

   if( rval == -1)       /* Undeciphered.  If 7 chars or less,  assume it's */
      {                  /* an observer-supplied desig.  8 or more chars,   */
      if( comet_desig)   /* start with '$' and put in up to eleven bytes of */
         {               /* the input name. */
         obj_name -= 2;
         len += 2;
         }
      if( len <= 7)
         memcpy( packed_desig + 5, obj_name, len);
      else
         {
         *packed_desig++ = '$';
         for( j = 0; j < 11 && obj_name[j]; j++)
            packed_desig[j] = obj_name[j];
         }
      }
   return( rval);
}

#ifdef TEST_MAIN
/* Compile with :

g++ -Wall -Wextra -pedantic -DTEST_MAIN -o mpc_fmt mpc_fmt.cpp date.cpp unpack.cpp */

int main( const int argc, const char **argv)
{
   if( argc > 1)
      {
      char buff[100];
      int rval;

      if( argc == 2)
         rval = create_mpc_packed_desig( buff, argv[1]);
      else
         rval = unpack_mpc_desig( buff, argv[1]);

      printf( "%2d: '%s'\n", rval, buff);
      }
   return( 0);
}
#endif
