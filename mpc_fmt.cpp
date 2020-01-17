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
         int divisor = 1;

         year = atoi( tbuff);
         month = get_two_digits( tbuff + 5);
//       rval = atof( tbuff + 8);
                     /* atof( ) is a little slow,  so we use a little more */
         for( i = 11; i < 17 && tbuff[i] != ' '; i++)  /* code in exchange */
            divisor *= 10;                             /* for better speed */
         rval = (double)get_two_digits( tbuff + 8) +
                           (double)atoi( tbuff + 11) / (double)divisor;
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
               + (double)get_two_digits( tbuff + 10) / minutes_per_day
               + (double)get_two_digits( tbuff + 12) / seconds_per_day;
            tbuff[13] = '.';
            rval += atof( tbuff + 13) / seconds_per_day;
            format_found = 20;      /* formats 20-23;  see above */
            start_of_decimals = 14;
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
   rval = atof( buff);
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
      rval += atof( buff + 3) / 60.;
      if( buff[5] == ' ' && isdigit( buff[7]))  /* i.e., seconds are given */
         {
         rval += atof( buff + 6) / 3600.;
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
      rval += atof( buff + 4) / 60. + atof( buff + 7) / 3600.;
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
           "zGSC",        /* no version specified */
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
           "UGAIA-DR1",
           "VGAIA-DR2",
           "WUCAC-5",
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

/* "Mutant hex" uses the usual hex digits 0123456789ABCDEF for numbers
0 to 15,  followed by G...Z for 16...35 and a...z for 36...61.  MPC stores
epochs and certain other numbers using this scheme to save space.  */

static char mutant_hex( const int ival)
{
   int rval = -1;

   if( ival >= 0)
      {
      if( ival < 10)
         rval = '0';
      else if( ival < 36)
         rval = 'A' - 10;
      else if( ival < 62)
         rval = 'a' - 36;
      }
   assert( rval >= 0);
   return( (char)( rval + ival));
}

/* create_mpc_packed_desig( ) takes a "normal" name for a comet/asteroid,
such as P/1999 Q1a or 2005 FF351,  and turns it into the 12-byte packed
format used in MPC reports and element files.  Documentation of this format
is given on the MPC Web site.  A 'test main' at the end of this file
shows the usage of this function.
   This should handle all "normal" asteroid and comet designations.
It doesn't handle natural satellites.  */

int create_mpc_packed_desig( char *packed_desig, const char *obj_name)
{
   int i, j, rval = 0;
   unsigned number;
   char buff[20], comet_desig = 0;

   while( *obj_name == ' ')
      obj_name++;

               /* Check for comet-style desigs such as 'P/1995 O1' */
               /* and such.  Leading character can be P, C, X, D, or A. */
               /* Or 'S' for natural satellites.  */
   if( strchr( "PCXDAS", *obj_name) && obj_name[1] == '/')
      {
      comet_desig = *obj_name;
      obj_name += 2;
      }

               /* Create a version of the name with all spaces removed: */
   for( i = j = 0; obj_name[i] && j < 19; i++)
      if( obj_name[i] != ' ')
         buff[j++] = obj_name[i];
   buff[j] = '\0';

   memset( packed_desig, ' ', 12);
   packed_desig[12] = '\0';
   number = atoi( buff);
   i = 0;
   while( isdigit( buff[i]))
      i++;
   if( buff[i] == 'P' && buff[i + 1] == '\0' && number < 10000)
      snprintf( packed_desig, 13, "%04uP       ", number);
               /* If the name starts with four digits followed by an */
               /* uppercase letter,  it's a provisional designation: */
   else if( number > 999 && number < 9000 && isupper( buff[4]))
      {
      int sub_designator;

      for( i = 0; i < 4; i++)
         {
         const char *surveys[4] = { "P-L", "T-1", "T-2", "T-3" };

         if( !strcmp( buff + 4, surveys[i]))
            {
            const char *surveys_packed[4] = {
                     "PLS", "T1S", "T2S", "T3S" };

            memcpy( packed_desig + 8, buff, 4);
            memcpy( packed_desig + 5, surveys_packed[i], 3);
            return( rval);
            }
         }

      snprintf( packed_desig + 5, 4, "%c%02d",
                  mutant_hex( number / 100), number % 100);
      packed_desig[6] = buff[2];    /* decade */
      packed_desig[7] = buff[3];    /* year */

      packed_desig[8] = (char)toupper( buff[4]);    /* prelim desigs are */
      i = 5;                                        /* _very_ scrambled  */
      if( isupper( buff[i]))                        /* when packed:      */
         {
         packed_desig[11] = buff[i];
         i++;
         }
      else
         packed_desig[11] = '0';

      sub_designator = atoi( buff + i);
      assert( sub_designator >= 0 && sub_designator < 620);
      packed_desig[10] = mutant_hex( sub_designator % 10);
      packed_desig[9] = mutant_hex( sub_designator / 10);
      if( comet_desig)
         {
         packed_desig[4] = comet_desig;
         while( isdigit( buff[i]))
            i++;
         if( buff[i] >= 'a' && buff[i] <= 'z')
            packed_desig[11] = buff[i];
         }
      }
   else if( !buff[i] && number < 620000
               && (!comet_desig || number < 10000))
      {                         /* simple numbered asteroid or comet */
      const int number = atoi( buff);

      if( comet_desig)
         sprintf( packed_desig, "%04d%c       ", number, comet_desig);
      else
         sprintf( packed_desig, "%c%04d       ", mutant_hex( number / 10000),
               number % 10000);
      }
   else                 /* strange ID that isn't decipherable.  For this, */
      {                 /* we just copy the first eleven non-space bytes, */
      while( j < 11)                              /* padding with spaces. */
         buff[j++] = ' ';
      buff[11] = '\0';
      *packed_desig = '~';
      strcpy( packed_desig + 1, buff);
      rval = -1;
      }
   return( rval);
}
