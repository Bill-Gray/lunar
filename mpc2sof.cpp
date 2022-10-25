/* Code to read in 'mpcorb.dat' and 'ELEMENTS.COMET' and output them in a
unified SOF (Standard Orbit Format) file.  See

http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET
https://www.projectpluto.com/orb_form.htm
https://minorplanetcenter.net/iau/MPCORB.html

The plan is that astcheck, Find_Orb,  astephem,  etc. will all be able to
use SOF,  thereby evading the numberless limitations baked into the MPCORB
and ASTORB formats, and allowing comets and asteroids to be expressed in a
unified format.

The asteroid elements can be in either 'mpcorb.dat' or 'MPCORB.DAT'.  For
Find_Orb,  the file should be placed in ~./find_orb.   I'll make sure that
other programs using this file (astcheck,  for example) look in that
directory as well.  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "watdefs.h"
#include "date.h"
#include "comets.h"
#include "stringex.h"

long extract_mpcorb_dat( ELEMENTS *elem, const char *buff);

const double PI =
   3.1415926535897932384626433832795028841971693993751058209749445923;
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

#define IS_POWER_OF_TWO( n)    (((n) & ((n)-1)) == 0)

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

static void extract_name( char *name, const char *iline)
{
   memset( name, ' ', 12);
   name[12] = '\0';
   if( iline[173] == ')')         /* numbered object */
      {
      size_t i;

      memcpy( name + 5, iline + 166, 7);
      for( i = 0; i < 12; i++)
         if( name[i] == '(')
            name[i] = ' ';
      }
   else           /* provisional desig */
      memcpy( name, iline + 175, 12);
   if( !memcmp( name, "            ", 12))   /* temp desig */
      memcpy( name, iline + 166, 12);
}


const char *sof_header =
       "Name        |Tp      .       |Te      |q          |"
       "i  .      |Om .      |om .      |e         |"
       "rms |n_o  |Tfirst  |Tlast   |Perts  |H .  |G . ^\n";

static void output_sof( const ELEMENTS *elem, char *obuff)
{
   const size_t obuff_size = 90;
   char perih_time[20], epoch_time[20];
   const int base_time_format = FULL_CTIME_YMD | FULL_CTIME_NO_SPACES
                | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_FORMAT_DAY
                | FULL_CTIME_LEADING_ZEROES;

   full_ctime( perih_time, elem->perih_time, base_time_format |
                                 FULL_CTIME_7_PLACES);
   full_ctime( epoch_time, elem->epoch, base_time_format);
   snprintf_err( obuff, obuff_size, "%s %s %11.8f ", perih_time, epoch_time, elem->q);
   snprintf_append( obuff, obuff_size, "%10.6f %10.6f %10.6f %10.8f ",
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
      fprintf( stderr, "Comments at the top of 'mpc2sof.cpp' should tell you\n"
                       "where to find this file.\n");
      exit( -1);
      }
   return( rval);
}

#define MAX_OUT 200

int main( const int argc, const char **argv)
{
   const size_t reclen = strlen( sof_header);
   char buff[400], *obuff = NULL;
   char tbuff[MAX_OUT];
   FILE *ifile = fopen( (argc > 1 ? argv[1] : "mpcorb.dat"), "rb");
   FILE *ofile = err_fopen( "mpcorb.sof", "wb");
   ELEMENTS elem;
   long epoch;
   int i;
   size_t n_out = 0, n_written;

   if( !ifile)
      ifile = err_fopen( "MPCORB.DAT", "rb");
// while( fgets( buff, sizeof( buff), ifile) && memcmp( buff, "------", 6))
//    ;           /* skip over header */
   fprintf( ofile, "%s", sof_header);
   while( fgets( buff, sizeof( buff), ifile))
      if( strlen( buff) == 203 &&
                       (epoch = extract_mpcorb_dat( &elem, buff)) > 0L)
         {
         char name[13], tfirst_buff[20];
         double jd;

         extract_name( name, buff);
         snprintf_err( tbuff, 14, "%-13s", name);
         output_sof( &elem, tbuff + 13);
         snprintf_append( tbuff, sizeof( tbuff), "%.4s %.5s ",
                      buff + 137, buff + 117);         /* rms, number obs */
         jd = get_time_from_string( 0., buff + 194, FULL_CTIME_YMD, NULL);
         if( !memcmp( buff + 132, "days", 4))
            jd -= atof( buff + 128);
         else
            {
            int year1, year2, n_scanned;

            n_scanned = sscanf( buff + 127, "%d-%d", &year1, &year2);
            assert( n_scanned == 2);
            assert( year1 > 1700);
            assert( year2 >= year1);
            assert( year2 < 2100);
            jd -= (double)( 365 * (year2 - year1 + 1));
            }
         full_ctime( tfirst_buff, jd, FULL_CTIME_YMD | FULL_CTIME_NO_SPACES
                           | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_DATE_ONLY
                           | FULL_CTIME_LEADING_ZEROES);
         snprintf_append( tbuff, sizeof( tbuff), "%.8s %.8s %.7s %.5s %.5s\n",
                  tfirst_buff, buff + 194,             /* Tfirst, Tlast */
                  buff + 142, buff + 8, buff + 14);    /* perts, H, G */
         assert( strlen( tbuff) == reclen);
         n_out++;
         if( IS_POWER_OF_TWO( n_out))
            obuff = (char *)realloc( obuff, n_out * 2 * reclen);
         memcpy( obuff + (n_out - 1) * reclen, tbuff, reclen);
         }
   fclose( ifile);

   if( argc <= 2 || strcmp( argv[2], "i"))
      {
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
            strcat( tbuff, "                                     \n");     /* Tlast, perts, H, G */
            assert( strlen( tbuff) == reclen);
            memcpy( obuff + n_out * reclen, tbuff, reclen);
            n_out++;
            }
      fclose( ifile);
      }
   n_written = fwrite( obuff, reclen, n_out, ofile);
   assert( n_written == n_out);
   free( obuff);
   fclose( ofile);
   return( 0);
}
