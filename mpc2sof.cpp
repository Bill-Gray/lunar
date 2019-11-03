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

#define DESIG_NUMBERED          0
#define DESIG_PROVISIONAL       1
#define DESIG_PLS_OR_TS         2
#define DESIG_NUMBERED_COMET    3
#define DESIG_PROVISIONAL_COMET 4
#define DESIG_UNRECOGNIZED     -1

static int desig_type( const char *buff)
{
   int rval = DESIG_UNRECOGNIZED;

   if( buff[0] == ' ' && buff[1] == ' '
               && buff[2] == ' ')
      rval = DESIG_NUMBERED;
   else if( buff[6] == '-' && (buff[5] == 'P' || buff[5] == 'T'))
      rval = DESIG_PLS_OR_TS;
   else if( buff[1] == '/')
      rval = DESIG_PROVISIONAL_COMET;
   else if( buff[3] == 'P' || buff[3] == 'D')
      rval = DESIG_NUMBERED_COMET;
   else if( buff[4] == ' ')
      if( buff[0] == '1' || buff[0] == '2')
         rval = DESIG_PROVISIONAL;
   return( rval);
}

int qsort_compare( const void *a, const void *b)
{
   const char *astr = (const char *)a;
   const char *bstr = (const char *)b;
   const int type1 = desig_type( astr);
   const int type2 = desig_type( bstr);
   int rval = type1 - type2;

   if( !rval)
      switch( type1)
         {
         case DESIG_NUMBERED:
         case DESIG_NUMBERED_COMET:
            break;
         case DESIG_PROVISIONAL:
            rval = memcmp( astr, bstr, 6);
            if( !rval)
               {
               const int na = (astr[7] == ' ' ? 0 : atoi( astr + 7));
               const int nb = (bstr[7] == ' ' ? 0 : atoi( bstr + 7));

               rval = na - nb;
               }
            if( !rval)
               rval = astr[6] - bstr[6];
            break;
         case DESIG_PLS_OR_TS:
            rval = memcmp( astr + 6, bstr + 6, 3);
            break;
         case DESIG_PROVISIONAL_COMET:
            rval = atoi( astr + 2) - atoi( bstr + 2);
            if( !rval)
               rval = memcmp( astr + 2, bstr + 2, 6);
            if( !rval)
               if( astr[11] < 'A' && bstr[11] < 'A')
                  rval = atoi( astr + 11) - atoi( bstr + 11);
            break;
         case DESIG_UNRECOGNIZED:
            rval = memcmp( astr, bstr, 12);
            break;
         }
   if( !rval)
      rval = strcmp( astr, bstr);
   return( rval);
}

#define MAX_ORBITS 1000000
#define MAX_OUT 200

int main( const int argc, const char **argv)
{
   const size_t reclen = strlen( sof_header);
   char buff[400], *obuff = (char *)calloc( MAX_ORBITS, reclen);
   char tbuff[MAX_OUT];
   FILE *ifile = fopen( (argc > 1 ? argv[1] : "mpcorb.dat"), "rb");
   FILE *ofile = err_fopen( "mpcorb.sof", "wb");
   ELEMENTS elem;
   long epoch;
   int i;
   size_t n_out = 0;

   if( !ifile)
      ifile = err_fopen( "MPCORB.DAT", "rb");
// while( fgets( buff, sizeof( buff), ifile) && memcmp( buff, "------", 6))
//    ;           /* skip over header */
   fprintf( ofile, "%s", sof_header);
   while( fgets( buff, sizeof( buff), ifile))
      if( strlen( buff) == 203 &&
                       (epoch = extract_mpcorb_dat( &elem, buff)) > 0L)
         {
         char name[13];

         extract_name( name, buff);
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
            strcat( tbuff, "                    \n");     /* Tlast, H, G */
            assert( strlen( tbuff) == reclen);
            memcpy( obuff + n_out * reclen, tbuff, reclen);
            n_out++;
            assert( n_out < MAX_ORBITS - 1);
            }
      }
   fclose( ifile);
   qsort( obuff, n_out, reclen, qsort_compare);
   fprintf( ofile, "%s", obuff);
   free( obuff);
   fclose( ofile);
   return( 0);
}
