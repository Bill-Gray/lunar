/* com_file.cpp: functions for merging comet orbital elems

Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "watdefs.h"
#include "afuncs.h"

/* 10 Apr 2003:  added code to handle cases such as
P/Linear-NEAT (2001 HT50): no "P/" or "C/" to go by.  Also,  revised the
code to understand that P/2001 HT50 is the same as C/2001 HT50,  or
D/2001 HT50,  or plain ol' 2001 HT50.  This was needed because
(for example) some comets change from "C/2003 A1" to "P/2003 A1",
once their periodic nature is determined.  So the '(letter)/' part
is now skipped over when determining if new items are replacements
for old ones. */

/* No-?:  parentheses don't match in this file */

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

#define BUFF_SIZE 212
#define MARGIN .5

/* Revised 30 Mar 2009:  this looks first for a permanent periodic comet
designation,  such as (145P).  If it finds it,  all is well;  otherwise,
it looks for designations such as (C/1999 X1) or (P/2003 RZ47).  If it
still hasn't found anything,  it looks for designations such as (2001 HT50).
I don't think such designations exist nowadays,  but there are probably
some old 'comets.dat' files that contain them.  NOTE that until 30 Mar 2009,
the code looked for the P/YYYY and C/YYYY designations _first_. */

int DLL_FUNC extract_periodic_name( const char *istr, char *ostr)
{
   int rval = 0, i;
   const char *loc;

   loc = strstr( istr, "P)");
   if( loc)          /* got a permanent numbered designation: */
      {
      const int end = (int)( loc - istr);

      for( i = end; i && istr[i] != '('; i--)
         ;
      if( istr[i] == '(')
         {
         memcpy( ostr, istr + i + 1, end - i);
         ostr[end - i] = '\0';
         rval = 1;
         }
      }

   for( i = 0; !rval && istr[i]; i++)
      if( istr[i] == 'P' || istr[i] == 'C')
         if( istr[i + 1] == '/')
            if( istr[i + 2] >= '1' && istr[i + 2] <= '9')
               {
               loc = istr + i;
               *ostr++ = *loc++;    /* copy the C or P */
               *ostr++ = *loc++;    /* copy the '/' */
               while( *loc == ' ')     /* skip possible irrelevant spaces */
                  loc++;
               while( *loc >= '0' && *loc <= '9')     /* copy in the year */
                  *ostr++ = *loc++;
               *ostr++ = ' ';
               while( *loc == ' ')  /* skip more possible irrelevant spaces */
                  loc++;
               while( *loc >= 'A' && *loc <= 'Z')  /* possible uppercase ID */
                  *ostr++ = *loc++;
               while( *loc >= '0' && *loc <= '9') /* possible numeric sub-ID */
                  *ostr++ = *loc++;
               if( *loc == '-')     /* ...followed,  possibly,  by '-(ID)'  */
                  {                 /* (happens in IMCCE 'eltnom.txt' file) */
                  *ostr++ = *loc++;
                  while( *loc >= 'A' && *loc <= 'Z')
                     *ostr++ = *loc++;
                  }
               *ostr = '\0';
               rval = 1;
               }
                     /* 10 Apr 2003:  added code to handle cases such as    */
                     /* P/Linear-NEAT (2001 HT50): no "P/" or "C/" to go by */
   for( loc = istr; !rval && *loc; loc++)
      if( *loc == '(')
         if( atoi( loc + 1) > 1000 && atoi( loc + 1) < 2300)
            {
            loc++;
            for( i = 0; loc[i] && loc[i] != ')'; i++)
               ;
            memcpy( ostr, loc, i);
            ostr[i] = '\0';
            rval = 1;
            }
   return( rval);
}

#ifdef OBSOLETE
      /* The following 'get_comet_file' function was used in Guide until */
      /* October 2012.  It used various methods to figure out what       */
      /* comets would reach 'mag_limit' over a given year,  outputting   */
      /* them to 'now_comt.dat'.  The whole thing is obsolete now,  but  */
      /* I'm  not quite ready to delete it and call it completely dead.  */

static void chop_cr_lf( char *buff)
{
   while( *buff && *buff != 10 && *buff != 13)
      buff++;
   *buff = '\0';
}

static void put_long_hundredths( char *buff, long val)
{
   int i;

   buff += 9;
   for( i = 0; val; i++)
      {
      *buff-- = (char)( val % 10L + '0');
      val /= 10L;
      if( i == 1)
         *buff-- = '.';
      }
}

int DLL_FUNC get_comet_file( const char *cd_path,
                 const double year, const double mag_limit)
{
   uint16_t FAR *limit_data;
   uint16_t max_time, min_time;
   char *buff;
   char FAR * FAR *periodics = NULL;
   static const char * const local_file = "now_comt.dat";
   static char * const big_comet_file = "asteroid\\cometg.dat";
   FILE *ifile = NULL, *ofile = NULL;
   int i, n_periodics = 0, rval = 0;
   int n_comets_in_file = 0, cometg_line_size;
   long lval1, lval2, lval3;

   buff = (char *)calloc( BUFF_SIZE, sizeof( char));
   if( !buff)
      return( -1);

   if( mag_limit >= 40.)      /* gotta use the whole thing */
      limit_data = NULL;
   else
      {
      long offset;

      ifile = fopen( "cometlim.bin", "rb");
      if( !ifile)
         {
         free( buff);
         return( -2);
         }
      fseek( ifile, 0L, SEEK_END);
      n_comets_in_file = (int)( ftell( ifile) / 16L);

      limit_data = (uint16_t FAR *)FCALLOC( n_comets_in_file,
                                    2 * sizeof( uint16_t));
      if( !limit_data)
         {
         fclose( ifile);
         free( buff);
         return( -3);
         }
      if( mag_limit > 25.)
         offset = (long)n_comets_in_file * 4L * 3L;
      else if( mag_limit > 20.)
         offset = (long)n_comets_in_file * 4L * 2L;
      else if( mag_limit > 15.)
         offset = (long)n_comets_in_file * 4L;
      else
         offset = 0L;
      fseek( ifile, offset, SEEK_SET);
#ifdef BITS_32
      fread( limit_data, n_comets_in_file, 2 * sizeof( uint16_t), ifile);
#else
      i = 0;
      while( i < n_comets_in_file)
         {
         int n_read = n_comets_in_file - i;

         if( n_read > 13)
            n_read = 13;
         fread( buff, 2 * sizeof( uint16_t), n_read, ifile);
         FMEMCPY( limit_data + i * 2, buff, n_read * 4);
         i += n_read;
         }
#endif
      fclose( ifile);
      ifile = NULL;
      }

   lval1 = (long)(( year - MARGIN) * 100. + .5);
   lval2 = (long)(( year + MARGIN) * 100. + .5);
   lval3 = (long)( mag_limit * 100. + .5);

   memset( buff, ' ', 30);
   buff[30] = '\n';
   put_long_hundredths( buff, lval1);
   put_long_hundredths( buff + 10, lval2);
   put_long_hundredths( buff + 20, lval3);

   ofile = fopen( local_file, "wb");
   if( !ofile)
      {
      rval = -4;
      goto The_End;
      }

   fwrite( buff, 31, 1, ofile);

   ifile = fopen( "comets.dat", "rb");
   if( ifile)
      {
      int n_comets;

      fgets( buff, BUFF_SIZE, ifile);
      n_comets = atoi( buff);
      while( n_comets--)
         {
         double ecc;
         int use_it = 1, is_periodic;
         char period_name[50];

         buff[152] = '\0';
         fgets( buff, BUFF_SIZE, ifile);
         chop_cr_lf( buff);
         ecc = atof( buff + 86);
         if( ecc < 1. && buff[152] != 'A')
            {
            double semi = atof( buff + 73) / (1. - ecc);
            double half_year = semi * sqrt( semi) / 2.;
            double t0 = atof( buff + 55) + atof( buff + 52) / 12.;

            if( t0 + half_year < year-MARGIN || t0 - half_year > year+MARGIN)
               use_it = 0;
            }
         is_periodic = extract_periodic_name( buff, period_name);
         if( is_periodic && period_name[1] == '/')  /* skip the 'P/' or 'C/' */
            memmove( period_name, period_name + 2, strlen( period_name));
         if( is_periodic)
            for( i = 0; i < n_periodics; i++)
               if( !FSTRCMP( period_name, periodics[i]))
                  use_it = 0;
         for( i = 151; i > 65; i--)
            if( buff[i] == '0' && buff[i + 1] == ' ')
               buff[i] = ' ';

         if( use_it)
            fprintf( ofile, "%s\n", buff);

         if( use_it && is_periodic)
            {
            if( !(n_periodics % 10))      /* gotta grow the array */
               {
               char FAR * FAR *new_arr;

               if( !n_periodics)          /* no,  gotta _make_ the array */
                  new_arr = (char FAR * FAR *)FCALLOC( 10, sizeof( char FAR *));
               else
                  new_arr = (char FAR * FAR *)FREALLOC( periodics,
                                (n_periodics + 10) * sizeof( char FAR *));
               if( !new_arr)
                  {
                  rval = -5;
                  goto The_End;
                  }
               periodics = new_arr;
               }
            periodics[n_periodics] = (char FAR *)FMALLOC( strlen( period_name) + 2);
            if( !periodics[n_periodics])
               {
               rval = -6;
               goto The_End;
               }
            FSTRCPY( periodics[n_periodics++], period_name);
            }
         }
      fclose( ifile);
      }

   min_time = (uint16_t)((year - MARGIN + 550.) * 20.);
   max_time = (uint16_t)((year + MARGIN + 550.) * 20.);

   ifile = fopen( big_comet_file + 9, "rb");      /* think locally... */
   if( !ifile)
      {
      strcpy( buff, cd_path);
      strcat( buff, big_comet_file);
      ifile = fopen( buff, "rb");        /* ...then CD */
      }
   if( !ifile)                /* real panic time! */
      {
      rval = -7;
      goto The_End;
      }

   fgets( buff, BUFF_SIZE, ifile);
   n_comets_in_file = atoi( buff);
   fgets( buff, BUFF_SIZE, ifile);
   cometg_line_size = strlen( buff);
   for( i = 0; i < n_comets_in_file; i++)
      if( !limit_data || (max_time > limit_data[i+i] &&
                          min_time < limit_data[i+i+1]))
         {
         int use_it = 1, j;
         char period_name[50];

         fseek( ifile, 6L + (long)i * (long)cometg_line_size, SEEK_SET);
         fread( buff, cometg_line_size, 1, ifile);
         buff[cometg_line_size] = '\0';
         if( extract_periodic_name( buff, period_name))
            {
            if( period_name[1] == '/')    /* skip the 'P/' or 'C/' */
               memmove( period_name, period_name + 2, strlen( period_name));
            for( j = 0; j < n_periodics; j++)
               if( !FSTRCMP( periodics[j], period_name))
                  use_it = 0;
            }

         if( use_it)
            {
            chop_cr_lf( buff);
               /* Added 21 Jan 97:  Not essential,  but it removes some */
               /* redundant zeroes in the COMETG data */
            if( cometg_line_size == 159)
               {
               memmove( buff + 137, buff + 143, strlen( buff + 142));
               for( j = 151; j > 65; j--)
                  if( buff[j] == '0' && buff[j + 1] == ' ')
                     buff[j] = ' ';
               }
            fprintf( ofile, "%s\n", buff);
            }
         }

The_End:

   if( ofile)
      fclose( ofile);
   if( ifile)
      fclose( ifile);
   if( limit_data)
      FFREE( limit_data);
   if( periodics)
      {
      for( i = 0; i < n_periodics; i++)
         if( periodics[i])
            FFREE( periodics[i]);
      FFREE( periodics);
      }
   free( buff);
   return( rval);
}
#endif           // #ifdef OBSOLETE
