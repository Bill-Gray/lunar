/* Code to read in 'mpcorb.dat' and 'ELEMENTS.COMET' and output them in a
unified SOF (Standard Orbit Format) file.  See

http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET
https://www.projectpluto.com/orb_form.htm
https://minorplanetcenter.net/iau/MPCORB.html

The plan is that astcheck, Find_Orb,  astephem,  etc. will all be able to
use SOF,  thereby evading the numberless limitations baked into the MPCORB
and ASTORB formats, and allowing comets and asteroids to be expressed in a
unified format. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "watdefs.h"
#include "date.h"
#include "comets.h"

long extract_mpcorb_dat( ELEMENTS *elem, const char *buff);

const double PI =
   3.1415926535897932384626433832795028841971693993751058209749445923;
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

static int parse_elements_dot_comet( ELEMENTS *elem, const char *buff)
{
   int rval = -1;

   if( strlen( buff) > 116 && buff[113] == '.')
      {
      char t_perih[15];

      memset( elem, 0, sizeof( ELEMENTS));
      elem->q = atof( buff + 51);
      elem->epoch = atof( buff + 44) + 2400000.5;
      elem->ecc = atof( buff + 64);
      elem->incl = atof( buff + 75) * PI / 180.;
      elem->arg_per = atof( buff + 85) * PI / 180.;
      elem->asc_node = atof( buff + 95) * PI / 180.;
      memcpy( t_perih, buff + 105, 14);
      t_perih[14] = '\0';
      elem->perih_time = get_time_from_string( 0., t_perih, 0, NULL);
      derive_quantities( elem, SOLAR_GM);
      rval = 0;
      }
   return( rval);
}

static int mutant_hex( const char ival)
{
   int rval = -1;

   if( ival >= '0' && ival <= '9')
      rval = ival - '0';
   else if( ival >= 'A' && ival <= 'Z')
      rval = ival + 10 - 'A';
   else if( ival >= 'a' && ival <= 'z')
      rval = ival + 36 - 'a';
   assert( rval >= 0 && rval < 62);
   return( rval);
}

static void decrypt_packed_desig( char *name, const char *packed)
{
   if( packed[5] == ' ')         /* numbered object */
      snprintf( name, 8, "%7d",
                      mutant_hex( *packed) * 10000 + atoi( packed + 1));
   else if( !memcmp( packed, "PLS", 3) || (packed[0] == 'T' && packed[2] == 'S'))
      {        /* Palomar-Leiden Survey or Trojan 1, 2, 3 surveys */
      snprintf( name, 9, "%.4s %c-%c", packed + 3, packed[0], packed[1]);
      }
   else                          /* provisional desig */
      snprintf( name, 13, "%2d%.2s %c%c%d%c", mutant_hex( *packed),
            packed + 1, packed[3], packed[6], mutant_hex( packed[4]),
            packed[5]);
}


const char *sof_header =
       "Name        |Tp      .       |Te      |q          |"
       "i  .      |Om .      |om .      |e         |"
       "rms |n_o  |Tlast   |H .  |G . ^\n";

static void output_sof( const ELEMENTS *elem, char *obuff)
{
   char perih_time[20], epoch_time[20];
   const int base_time_format = FULL_CTIME_YMD | FULL_CTIME_NO_SPACES
                | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_FORMAT_DAY
                | FULL_CTIME_LEADING_ZEROES;

   full_ctime( perih_time, elem->perih_time, base_time_format |
                                 FULL_CTIME_7_PLACES);
   full_ctime( epoch_time, elem->epoch, base_time_format);
   sprintf( obuff, "%s %s %11.8f ", perih_time, epoch_time, elem->q);
   snprintf( obuff + strlen( obuff), 45, "%10.6f %10.6f %10.6f %10.8f ",
                  elem->incl * 180. / PI,  elem->asc_node * 180. / PI,
                  elem->arg_per * 180. / PI, elem->ecc);
}

static FILE *err_fopen( const char *filename, const char *permits)
{
   FILE *rval = fopen( filename, permits);

   if( !rval)
      {
      fprintf( stderr, "Failed to open '%s'", filename);
      perror( " ");
      exit( -1);
      }
   return( rval);
}

int qsort_compare( const void *a, const void *b)
{
   return( strcmp( (const char *)a, (const char *)b));
}

#define MAX_ORBITS 1000000
#define MAX_OUT 200

int main( const int argc, const char **argv)
{
   const size_t reclen = strlen( sof_header);
   char buff[400], *obuff = (char *)calloc( MAX_ORBITS, reclen);
   char tbuff[MAX_OUT];
   FILE *ifile = err_fopen( (argc > 1 ? argv[1] : "mpcorb.dat"), "rb");
   ELEMENTS elem;
   long epoch;
   int i;
   size_t n_out = 0;

   while( fgets( buff, sizeof( buff), ifile)
            && memcmp( buff, "------", 6))
      ;           /* skip over header */
   printf( "%s", sof_header);
   while( fgets( buff, sizeof( buff), ifile))
      if( strlen( buff) == 203 &&
                       (epoch = extract_mpcorb_dat( &elem, buff)) > 0L)
         {
         char name[30];

         decrypt_packed_desig( name, buff);
         snprintf( tbuff, 14, "%-13s", name);
         output_sof( &elem, tbuff + 13);
         sprintf( tbuff + strlen( tbuff), "%.4s %.5s ",
                      buff + 137, buff + 117);         /* rms, number obs */
         sprintf( tbuff + strlen( tbuff), "%.8s %.5s %.5s\n",
                  buff + 194, buff + 8, buff + 14);        /* Tlast, H, G */
         assert( strlen( tbuff) == reclen);
         memcpy( obuff + n_out * reclen, tbuff, reclen);
         n_out++;
         assert( n_out < MAX_ORBITS - 1);
         }
   fclose( ifile);

   ifile = err_fopen( (argc > 2 ? argv[2] : "ELEMENTS.COMET"), "rb");
   for( i = 0; i < 2; i++)       /* ELEMENTS.COMET has two header lines */
      if( !fgets( buff, sizeof( buff), ifile))
         {
         fprintf( stderr, "Failed to read line from ELEMENTS.COMET\n");
         return( -2);
         }
   while( fgets( buff, sizeof( buff), ifile))
      if( !parse_elements_dot_comet( &elem, buff) && elem.epoch > 1721000.)
         {
         char *tptr, name[13];

         memset( name, 0, sizeof( name));
         if( buff[2] == ' ')    /* provisional desig */
            {
            memcpy( name, buff + 3, 12);
            tptr = strstr( name, " (");
            if( tptr)
               *tptr = '\0';
            }
         else           /* permanent desig */
            {
            memcpy( name, buff, 4);
            tptr = buff + 43;
            while( *tptr == ' ')       /* search for end of name */
               tptr--;
            if( tptr[-1] == '-')       /* one-character suffix */
               memcpy( name + 4, tptr - 1, 2);
            if( tptr[-2] == '-')       /* two-character suffix */
               memcpy( name + 4, tptr - 2, 3);
            }
         snprintf( tbuff, 14, "%-12s ", name);
         output_sof( &elem, tbuff + 13);
         strcat( tbuff, "           ");    /* rms, number obs */
         strcat( tbuff, "                    \n");     /* Tlast, H, G */
         assert( strlen( tbuff) == reclen);
         memcpy( obuff + n_out * reclen, tbuff, reclen);
         n_out++;
         assert( n_out < MAX_ORBITS - 1);
         }
   fclose( ifile);
   qsort( obuff, n_out, reclen, qsort_compare);
   printf( "%s", obuff);
   free( obuff);
   return( 0);
}
