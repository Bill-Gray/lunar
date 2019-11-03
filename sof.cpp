/* sof.cpp: extracts orbital elements from SOF (Standard Orbit Format)
(See https://www.projectpluto.com/orb_form.htm for an explanation)

Copyright (C) 2018, Project Pluto

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

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "watdefs.h"
#include "comets.h"
#include "date.h"

int extract_sof_data( ELEMENTS *elem, const char *buff, const char *header);

#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

#define SOF_Q_FOUND                 0x0001
#define SOF_ECC_FOUND               0x0002
#define SOF_NAME_FOUND              0x0004
#define SOF_TPERIH_FOUND            0x0008
#define SOF_TEPOCH_FOUND            0x0010
#define SOF_INCL_FOUND              0x0020
#define SOF_ASC_NODE_FOUND          0x0040
#define SOF_ARG_PERIH_FOUND         0x0080
#define SOF_ABS_MAG_FOUND           0x0100
#define SOF_SLOPE_PARAM_FOUND       0x0200

   /* If you don't have these fields,  you don't have a defined orbit. */
#define MIN_FIELDS_NEEDED (SOF_Q_FOUND | SOF_ECC_FOUND | SOF_TPERIH_FOUND \
              | SOF_INCL_FOUND | SOF_ASC_NODE_FOUND | SOF_ARG_PERIH_FOUND)

/* All dates are currently stored as YYYYMMDD[.dddd],  Gregorian. */

static double extract_jd( const char *buff)
{
   long t;
   int bytes_read;
   double rval = 0.;

   if( sscanf( buff, "%ld%n", &t, &bytes_read) == 1)
      {
      if( bytes_read == 8)    /* YYYYMMDD */
         rval = (double)dmy_to_day( t % 100, (t / 100) % 100, t / 10000,
                                    CALENDAR_GREGORIAN) - .5;
      if( buff[bytes_read] == '.')
         rval += atof( buff + bytes_read);
      }
   return( rval);
}

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int extract_sof_data( ELEMENTS *elem, const char *buff, const char *header)
{
   int fields_found = 0, rval;

// if( strlen( buff) != strlen( header))
//    return( -1);
   memset( elem, 0, sizeof( ELEMENTS));
   elem->slope_param = 0.15;
   elem->gm = SOLAR_GM;
   while( *header >= ' ')
      {
      size_t i = 0;

      while( header[i] >= ' ' && header[i] != '|')
         i++;
      if( i < 2)
         return( -1);
      if( header[1] == ' ')
         {
         switch( header[0])
            {
            case 'q':
               elem->q = atof( buff);
               fields_found |= SOF_Q_FOUND;
               break;
            case 'e':
               elem->ecc = atof( buff);
               fields_found |= SOF_ECC_FOUND;
               break;
            case 'i':
               elem->incl = atof( buff) * PI / 180.;
               fields_found |= SOF_INCL_FOUND;
               break;
            case 'O':
               elem->asc_node = atof( buff) * PI / 180.;
               fields_found |= SOF_ASC_NODE_FOUND;
               break;
            case 'o':
               elem->arg_per = atof( buff) * PI / 180.;
               fields_found |= SOF_ARG_PERIH_FOUND;
               break;
            case 'H':
               elem->abs_mag = atof( buff);
               fields_found |= SOF_ABS_MAG_FOUND;
               break;
            case 'G':
               elem->slope_param = atof( buff);
               fields_found |= SOF_SLOPE_PARAM_FOUND;
               break;
            case 'C':
               elem->central_obj = atoi( buff);
               break;
            }
         }
      else if( header[2] == ' ')
         {
         switch( header[0])
            {
            case 'T':
               {
               const double jd = extract_jd( buff);

               if( header[1] == 'p')
                  {
                  elem->perih_time = jd;
                  fields_found |= SOF_TPERIH_FOUND;
                  }
               else if( header[1] == 'e')
                  {
                  elem->epoch = jd;
                  fields_found |= SOF_TEPOCH_FOUND;
                  }
               fields_found |= SOF_Q_FOUND;
               }
               break;
            case 'O':
               elem->asc_node = atof( buff) * PI / 180.;
               fields_found |= SOF_ASC_NODE_FOUND;
               break;
            case 'o':
               elem->arg_per = atof( buff) * PI / 180.;
               fields_found |= SOF_ARG_PERIH_FOUND;
               break;
            }

         }
      if( header[i] == '|')
         i++;
      header += i;
      buff += i;
      }
   if( (fields_found & MIN_FIELDS_NEEDED) == MIN_FIELDS_NEEDED)
      {
      rval = 0;            /* success */
      derive_quantities( elem, elem->gm);
      }
   else
      {
      printf( "Got '%x'\n", (unsigned)fields_found);
      rval = -1;
      }
   return( rval);
}

#ifdef TEST_CODE

#define MAX_LEN 300

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( argv[1], "rb");
   char header_line[MAX_LEN], buff[MAX_LEN];
   char *tptr;
   size_t name_len;

   assert( ifile);
   assert( argc == 3);
   if( !fgets( header_line, sizeof( header_line), ifile))
      {
      printf( "Couldn't read header\n");
      return( -1);
      }
   tptr = strchr( header_line, '|');
   assert( tptr);
   name_len = tptr - header_line;
   while( fgets( buff, sizeof( buff), ifile))
      if( (tptr = strstr( buff, argv[2])) && tptr < buff + name_len
                  && tptr[strlen( argv[2])] == ' ')
         {
         ELEMENTS elem;
         int rval = extract_sof_data( &elem, buff, header_line);

         printf( "rval %d\n", rval);
         printf( "   Tp= %f\n", elem.perih_time);
         printf( "   q = %f\n", elem.q);
         printf( "   a = %f\n", elem.major_axis);
         printf( "   e = %f\n", elem.ecc);
         printf( "   H = %f\n", elem.abs_mag);
         printf( "   i = %f\n", elem.incl * 180. / PI);
         printf( "   Om= %f\n", elem.asc_node * 180. / PI);
         printf( "   om= %f\n", elem.arg_per * 180. / PI);
         printf( "   Epoch = %f\n", elem.epoch);
         }
   fclose( ifile);
   return( 0);
}
#endif
