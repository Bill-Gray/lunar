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

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "mpc_func.h"

#if defined(_MSC_VER) && _MSC_VER < 1900
                      /* For older MSVCs,  we have to supply our own  */
                      /* snprintf().  See snprintf.cpp for details.  */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

#define MAX_DEPTH 20
#define PIECE_SIZE   25

typedef struct
{
   int depth, tags[MAX_DEPTH];
   char line[83];    /* allow possible CR, LF,  & null */
   char line2[83];
   char rms_ra[PIECE_SIZE], rms_dec[PIECE_SIZE], corr[PIECE_SIZE];
   char rms_mag[PIECE_SIZE], rms_time[PIECE_SIZE], center[PIECE_SIZE];
   int id_set, getting_lines;
   int prev_line_passed_through;
   int prev_rval, n_psv_fields;
   int *psv_tags;
} ades2mpc_t;

void *init_ades2mpc( void)
{
   ades2mpc_t *rval = (ades2mpc_t *)calloc( 1, sizeof( ades2mpc_t));

   return( rval);
}

/* The logic of processing ADES data to 80-column punched-card format is
somewhat convoluted.  We have to take a possible mix of PSV and XML ADES
with 80-column data,  plus possible extraneous text,  and convert it to
a modified version of the MPC's format.  The 'modification' is that
the times and RA/decs can be given to greater precision,  and uncertainties
expressed with comment lines such as

COM Sigmas 0.5x0.4,0.93 m:0.02 t:0.5

to indicate RA sigma = 0.5 arcsec,  dec simga = 0.4,  that the correlation
between them is 0.93,  magnitude uncertainty is 0.02,  and the timing
uncertainty is 0.5 seconds.  Note that the above comment can pass
through MPC's processing,  but the actual information is left present
for Find_Orb and other software.

   But it _does_ mean that a single PSV line may be output as two or
more lines (the above sigmas,  and finally an 80-column astrometric
record).  With XML,  we have to read in many lines,  gradually assembling
the sigmas and the 80-column output;  only when we reach the closing
</optical> tag can we cough up two or more lines.  If it's an MPC
line,  we can pass it through unaltered... _unless_ it's a Dave Tholen
style 92-column line,  with 80 "traditional" columns and sigmas in
the last twelve columns.

   The above 'ades2mpc_t' context structure stores the various bits
of the output line(s).  When we read a PSV line or hit the </optical>
XML closing tag,  we set 'getting_lines' to 1,  i.e.,  "we've got
the data we need;  start outputting corresponding lines such as
the above comment and 80-column MPC data."

   At that point,  each call to xlate_ades2mpc( ) gets passed
along to get_a_line( ),  which checks to see if we have some sigmas
or (for satellite observations) an object center or one or
(for satellite,  radar,  and roving observations) two 80-column
lines to be output.  After we've done all that,  'getting_lines'
is reset to zero,  and we're ready to accept the next observation(s). */

/* In theory,  when we reach the end of the file,  any 'state' flags that
got set (because of encountering ADES tags) should be cleared (because of
encountering closing tags).  */

int free_ades2mpc_context( void *context)
{
   ades2mpc_t *tptr = (ades2mpc_t *)context;
   const int rval = tptr->depth;

   if( tptr->psv_tags)
      free( tptr->psv_tags);
   free( context);
   return( rval);
}

   /* see 'adestags.c' for code that generates these #defines.  Tags are
sorted so that a binary search can speed up tag matching slightly. */

#define ADES_Location                      1
#define ADES_MPCID                         2
#define ADES_OffsetVal                     3
#define ADES_OpticalID                     4
#define ADES_OpticalRes                    5
#define ADES_OpticalResMag                 6
#define ADES_OpticalResiduals              7
#define ADES_Photometry                    8
#define ADES_Precision                     9
#define ADES_RadarID                      10
#define ADES_RadarResiduals               11
#define ADES_RadarValue                   12
#define ADES_ades                         13
#define ADES_aperture                     14
#define ADES_arrayScale                   15
#define ADES_artSat                       16
#define ADES_astCat                       17
#define ADES_astrometry                   18
#define ADES_band                         19
#define ADES_biasDec                      20
#define ADES_biasMag                      21
#define ADES_biasRA                       22
#define ADES_biasTime                     23
#define ADES_coinvestigators              24
#define ADES_collaborators                25
#define ADES_com                          26
#define ADES_comment                      27
#define ADES_ctr                          28
#define ADES_dec                          29
#define ADES_decStar                      30
#define ADES_delay                        31
#define ADES_deltaDec                     32
#define ADES_deltaRA                      33
#define ADES_deprecated                   34
#define ADES_design                       35
#define ADES_detector                     36
#define ADES_disc                         37
#define ADES_dist                         38
#define ADES_doppler                      39
#define ADES_exp                          40
#define ADES_fRatio                       41
#define ADES_filter                       42
#define ADES_fitOrder                     43
#define ADES_frq                          44
#define ADES_fundingSource                45
#define ADES_institution                  46
#define ADES_line                         47
#define ADES_localUse                     48
#define ADES_logSNR                       49
#define ADES_mag                          50
#define ADES_measurers                    51
#define ADES_mode                         52
#define ADES_mpcCode                      53
#define ADES_nStars                       54
#define ADES_name                         55
#define ADES_notes                        56
#define ADES_nucMag                       57
#define ADES_objectDetection              58
#define ADES_obsBlock                     59
#define ADES_obsCenter                    60
#define ADES_obsContext                   61
#define ADES_obsData                      62
#define ADES_obsID                        63
#define ADES_obsTime                      64
#define ADES_observatory                  65
#define ADES_observers                    66
#define ADES_occultation                  67
#define ADES_offset                       68
#define ADES_optical                      69
#define ADES_opticalResidual              70
#define ADES_orbID                        71
#define ADES_orbProd                      72
#define ADES_pa                           73
#define ADES_permID                       74
#define ADES_photAp                       75
#define ADES_photCat                      76
#define ADES_photMod                      77
#define ADES_photProd                     78
#define ADES_photometry                   79
#define ADES_pixelScale                   80
#define ADES_pos1                         81
#define ADES_pos2                         82
#define ADES_pos3                         83
#define ADES_posCov11                     84
#define ADES_posCov12                     85
#define ADES_posCov13                     86
#define ADES_posCov22                     87
#define ADES_posCov23                     88
#define ADES_posCov33                     89
#define ADES_precDec                      90
#define ADES_precRA                       91
#define ADES_precTime                     92
#define ADES_prog                         93
#define ADES_provID                       94
#define ADES_ra                           95
#define ADES_raStar                       96
#define ADES_radar                        97
#define ADES_radarResidual                98
#define ADES_rcv                          99
#define ADES_ref                         100
#define ADES_remarks                     101
#define ADES_resDec                      102
#define ADES_resDelay                    103
#define ADES_resDoppler                  104
#define ADES_resMag                      105
#define ADES_resRA                       106
#define ADES_rmsCorr                     107
#define ADES_rmsDec                      108
#define ADES_rmsDelay                    109
#define ADES_rmsDist                     110
#define ADES_rmsDoppler                  111
#define ADES_rmsFit                      112
#define ADES_rmsMag                      113
#define ADES_rmsPA                       114
#define ADES_rmsRA                       115
#define ADES_rmsTime                     116
#define ADES_seeing                      117
#define ADES_selAst                      118
#define ADES_selDelay                    119
#define ADES_selDoppler                  120
#define ADES_selPhot                     121
#define ADES_sigCorr                     122
#define ADES_sigDec                      123
#define ADES_sigDelay                    124
#define ADES_sigDoppler                  125
#define ADES_sigMag                      126
#define ADES_sigRA                       127
#define ADES_sigTime                     128
#define ADES_software                    129
#define ADES_stn                         130
#define ADES_subFmt                      131
#define ADES_subFrm                      132
#define ADES_submitter                   133
#define ADES_sys                         134
#define ADES_telescope                   135
#define ADES_trkID                       136
#define ADES_trkSub                      137
#define ADES_trx                         138
#define ADES_uncTime                     139

/* See 'adestags.c' for code that created following array
and the above #defines. */

static int find_tag( const char *buff, size_t len)
{
   static const char *tags[] = {
       "Location", "MPCID", "OffsetVal", "OpticalID",
       "OpticalRes", "OpticalResMag", "OpticalResiduals",
       "Photometry", "Precision", "RadarID", "RadarResiduals",
       "RadarValue", "ades", "aperture", "arrayScale", "artSat",
       "astCat", "astrometry", "band", "biasDec", "biasMag",
       "biasRA", "biasTime", "coinvestigators", "collaborators",
       "com", "comment", "ctr", "dec", "decStar", "delay",
       "deltaDec", "deltaRA", "deprecated", "design", "detector",
       "disc", "dist", "doppler", "exp", "fRatio", "filter",
       "fitOrder", "frq", "fundingSource", "institution", "line",
       "localUse", "logSNR", "mag", "measurers", "mode",
       "mpcCode", "nStars", "name", "notes", "nucMag",
       "objectDetection", "obsBlock", "obsCenter", "obsContext",
       "obsData", "obsID", "obsTime", "observatory", "observers",
       "occultation", "offset", "optical", "opticalResidual",
       "orbID", "orbProd", "pa", "permID", "photAp", "photCat",
       "photMod", "photProd", "photometry", "pixelScale", "pos1",
       "pos2", "pos3", "posCov11", "posCov12", "posCov13",
       "posCov22", "posCov23", "posCov33", "precDec", "precRA",
       "precTime", "prog", "provID", "ra", "raStar", "radar",
       "radarResidual", "rcv", "ref", "remarks", "resDec",
       "resDelay", "resDoppler", "resMag", "resRA", "rmsCorr",
       "rmsDec", "rmsDelay", "rmsDist", "rmsDoppler", "rmsFit",
       "rmsMag", "rmsPA", "rmsRA", "rmsTime", "seeing", "selAst",
       "selDelay", "selDoppler", "selPhot", "sigCorr", "sigDec",
       "sigDelay", "sigDoppler", "sigMag", "sigRA", "sigTime",
       "software", "stn", "subFmt", "subFrm", "submitter", "sys",
       "telescope", "trkID", "trkSub", "trx", "uncTime",
   NULL };

   int i, rval = -1;

   if( *buff == '/')       /* closing tag */
      {
      buff++;
      len--;
      }
   if( !memcmp( buff, "ades", 4))
      rval = 0;
   for( i = 0; tags[i] && rval == -1; i++)
      if( !memcmp( tags[i], buff, len) && !tags[i][len])
         rval = i + 1;
   return( rval);
}

static inline void move_fits_time( char *optr, const char *iptr)
{
   if( iptr[4] == '-' && iptr[7] == '-' && iptr[10] == 'T'
               && iptr[13] == ':')
      {
      *optr++ = (char)((iptr[0] * 10 + iptr[1] - '0' * 11) + 'A' - 10);
      iptr += 2;
      while( *iptr && *iptr != 'Z')
         {
         if( *iptr >= '0' && *iptr <= '9')
            *optr++ = *iptr;
         else if( *iptr == 'T')
            *optr++ = ':';
         iptr++;
         }
      }
}

static inline void place_value( char *optr, const char *iptr, const size_t ilen,
             size_t leading_places)
{
   size_t i, point_loc = ilen;

   for( i = 0; i < ilen; i++)
      if( iptr[i] == '.')
         point_loc = i;
   assert( leading_places >= point_loc);
   while( leading_places > point_loc)
      {
      *optr++ = '0';
      leading_places--;
      }
   memcpy( optr, iptr, (ilen < 9 + leading_places) ? ilen : 9 + leading_places);
}

static const char *skip_whitespace( const char *tptr)
{
   while( *tptr && isspace( *tptr))
      tptr++;
   return( tptr);
}

static int get_a_line( char *obuff, ades2mpc_t *cptr)
{
   if (cptr->rms_ra[0])
      {
      sprintf( obuff, "COM Sigmas %s", cptr->rms_ra);
      if( strcmp( cptr->rms_ra, cptr->rms_dec))
         {
         sprintf( obuff + strlen( obuff), "x%s", cptr->rms_dec);
         if( atof( cptr->corr))
            sprintf( obuff + strlen( obuff), ",%s", cptr->corr);
         }
      if( cptr->rms_mag[0])
         {
         strcat( obuff, " m:");
         strcat( obuff, cptr->rms_mag);
         cptr->rms_mag[0] = '\0';
         }
      if( cptr->rms_time[0])
         {
         strcat( obuff, " t:");
         strcat( obuff, cptr->rms_time);
         cptr->rms_mag[0] = '\0';
         }
      cptr->rms_ra[0] = '\0';
      strcat( obuff, "\n");
      }
   else if( cptr->center[0])
      {
      sprintf( obuff, "COM Offset center %s", cptr->center);
      cptr->center[0] = '\0';
      }
   else if( cptr->line[0])
      {
      strcpy( obuff, cptr->line);
      cptr->line[0] = '\0';
      }
   else if( cptr->line2[0])
      {
      strcpy( obuff, cptr->line2);
      cptr->line2[0] = '\0';
      }
   else                         /* got all the lines we'll get */
      cptr->getting_lines = 0;
   return( cptr->getting_lines);
}

static int process_ades_tag( char *obuff, ades2mpc_t *cptr, const int itag,
                 const char *tptr, size_t len)
{
   int rval = 0;
   char name[40];

   if( len < sizeof( name))
      {
      memcpy( name, tptr, len);
      name[len] = '\0';
      }
   switch( itag)
      {
      case ADES_mpcCode:
         sprintf( obuff, "COD %s\n", name);
         rval = 1;
         break;
      case ADES_name:
      case ADES_line:
      case ADES_institution:
         {
         const char *format = NULL;

         if( cptr->depth > 1)
            switch( cptr->tags[cptr->depth - 2])
               {
               case ADES_mpcCode:
               case ADES_comment:
                  format = "COM %.*s\n";
                  break;
               case ADES_observatory:
                  format = "COM Observatory name %.*s\n";
                  break;
               case ADES_observers:
                  format = "OBS %.*s\n";
                  break;
               case ADES_measurers:
                  format = "MEA %.*s\n";
                  break;
               case ADES_submitter:
                  format = "CON %.*s\n";
                  break;
               default:
                  format = "COM Mangled '%.*s'\n";
                  break;
               }
         if( format)
            sprintf( obuff, format, (int)len, tptr);
         else
            strcpy( obuff, "COM Mangled name data\n");
         rval = 1;
         }
         break;
      case ADES_stn:
         memcpy( cptr->line + 77, tptr, 3);
         break;
      case ADES_obsTime:
         move_fits_time( cptr->line + 15, tptr);
         break;
      case ADES_band:
         cptr->line[70] = *tptr;
         break;
      case ADES_mode:
         if( len == 3)
            {
            const char *modes = "CCCD VVID PPHO eENC pPMT MMIC TMER ";
            int i;

            for( i = 0; modes[i]; i += 5)
               if( !memcmp( modes + i + 1, tptr, 3))
                  cptr->line[14] = modes[i];
            }
         break;
      case ADES_prog:
         {
         const char *programs = "0123456789!\"#$%&'()*+,-./[\\]^_`{|}~";
         const int idx = atoi( tptr);

         if( idx >=0 && idx <= 34)
            cptr->line[13] = programs[idx];
         }
         break;
      case ADES_sys:
         cptr->line2[0] = ' ';
         if( !strcmp( name, "ICRF_KM"))
            cptr->line2[32] = '1';
         else if( !strcmp( name, "ICRF_AU"))
            cptr->line2[32] = '2';
         else
            {
            strcpy( obuff, "Bad <sys> tag\n");
            assert( 1);
            }
         cptr->line2[14] = 's';
         cptr->line[14] = 'S';
         break;
      case ADES_ctr:
         assert( len < sizeof( cptr->center));
         strcpy( cptr->center, name);
         break;
      case ADES_pos1:
      case ADES_pos2:
      case ADES_pos3:
         {
         int dec_loc = 40 + (itag - ADES_pos1) * 12;
         int nlen = (int)len;
         char *tptr2;

                  /* cvt scientific notation,  if any : */
         if( strchr( name, 'e') || strchr( name, 'E'))
            {
            snprintf( name, sizeof( name), "%.13f", atof( name));
            nlen = 12;
            }
         if( *name != '+' && *name != '-')   /* no sign provided; */
            {                             /* insert one */
            nlen++;
            memmove( name + 1, name, nlen);
            *name = '+';
            }
         cptr->line2[dec_loc - 6] = *name;
         tptr2 = strchr( name, '.');
         assert( tptr2);
         dec_loc -= (int)(tptr2 - name);
         if( cptr->line2[32] == '2')
            dec_loc -= 4;
         else if( cptr->line2[32] != '1')
            {
            strcpy( obuff, "Bad posn data\n");
            rval = 1;
            }
         memcpy( &cptr->line2[dec_loc + 1], name + 1,
                                        (nlen > 12 ? 11 : nlen - 1));
         }
         break;
      case ADES_ra:
         if( *tptr == '+')
            {
            tptr++;
            len--;
            }
         place_value( cptr->line + 32, tptr, len, 3);
         break;
      case ADES_dec:
         if( *tptr == '-' || *tptr == '+')
            {
            cptr->line[44] = *tptr++;
            len--;
            }
         else
            cptr->line[44] = '+';
         place_value( cptr->line + 45, tptr, len, 2);
         break;
      case ADES_astCat:
         assert( len < sizeof( name));
         cptr->line[71] = net_name_to_byte_code( name);
         break;
      case ADES_rmsRA:
         assert( len < sizeof( cptr->rms_ra));
         strcpy( cptr->rms_ra, name);
         break;
      case ADES_rmsDec:
         assert( len < sizeof( cptr->rms_dec));
         strcpy( cptr->rms_dec, name);
         break;
      case ADES_rmsCorr:
         assert( len < sizeof( cptr->corr));
         strcpy( cptr->corr, name);
         break;
      case ADES_uncTime:
         assert( len < sizeof( cptr->rms_time));
         strcpy( cptr->rms_time, name);
         break;
      case ADES_rmsMag:
         assert( len < sizeof( cptr->rms_mag));
         strcpy( cptr->rms_mag, name);
         break;
      case ADES_provID:
         if( cptr->id_set == ADES_permID)
            break;                  /* FALLTHRU */
      case ADES_permID:
         {
         char tbuff[20];

         assert( len < sizeof( name));
         create_mpc_packed_desig( tbuff, name);
         memcpy( cptr->line, tbuff, 12);
         cptr->id_set = itag;
         }
         break;
      case ADES_artSat:
         assert( len < 13);
         if( !cptr->id_set)
            {
            cptr->id_set = ADES_artSat;
            memcpy( cptr->line, tptr, len);
            }
         break;
      case ADES_trkSub:
         assert( len < 13);
         if( !cptr->id_set)
            {
            cptr->id_set = ADES_trkSub;
            if( len < 8)
               memcpy( cptr->line + 5, tptr, len);
            else
               memcpy( cptr->line + 12 - len, tptr, len);
            }
         break;
      case ADES_mag:
         memcpy( cptr->line + 65, tptr, (len < 5) ? len : 5);
         break;
      }
   return( rval);
}

static int process_psv_tag( ades2mpc_t *cptr, char *obuff, const int itag,
                                       const char *ibuff)
{
   int rval = 0;
   size_t len = 0;

   ibuff = skip_whitespace( ibuff);
   while( ibuff[len] != '|' && ibuff[len] >= ' ')
      len++;
   while( len && ibuff[len - 1] == ' ')
      len--;            /* drop trailing spaces */
   if( len && itag >= 0)
      rval = process_ades_tag( obuff, cptr, itag, ibuff, len);
   return( rval);
}

static void setup_observation( ades2mpc_t *cptr)
{
   memset( cptr->line, ' ', 80);
   strcpy( cptr->line + 80, "\n");
   memset( cptr->line2, ' ', 80);
   strcpy( cptr->line2 + 80, "\n");
   cptr->line2[0] = '\0';
   cptr->id_set = 0;
}

/* The following checks to see if the input is astrometry with "Dave
Tholen sigmas",  such as

2019 SH3      C2019 09 26.56792500 00 13.041-05 25 38.05         18.6 G      T12 0.050 0.049

   If it is,  the line is converted to two lines with the designation in the
packed form desired by MPC,  and the sigmas on a separate comment line :

COM Sigmas 0.050x0.049
     K19S03H  C2019 09 26.56792500 00 13.041-05 25 38.05         18.6 G      T12
*/

static int check_for_tholen_sigmas( ades2mpc_t *cptr, char *obuff, const char *ibuff)
{
   int rval = 0;

   if( strlen( ibuff) >= 92 && (unsigned char)ibuff[92] < ' '
          && ibuff[82] == '.' && ibuff[88] == '.'
          && ibuff[80] == ' ' && ibuff[86] == ' '
          && is_valid_mpc_code( ibuff + 77))
      {
      memcpy( obuff, ibuff, 80);
      obuff[80] = '\0';
      if( extract_date_from_mpc_report( obuff, NULL))
         {
         char packed_desig[13];

         setup_observation( cptr);
         memcpy( cptr->rms_ra, ibuff + 81, 5);
         memcpy( cptr->rms_dec, ibuff + 87, 5);
         cptr->rms_ra[5] = cptr->rms_dec[5] = '\0';
         strcpy( cptr->line, obuff);
         obuff[12] = '\0';
         create_mpc_packed_desig( packed_desig, obuff);
         memcpy( cptr->line, packed_desig, 12);
         cptr->getting_lines = rval = 1;
         }
      }
   return( rval);
}

static int process_psv_line( ades2mpc_t *cptr, char *obuff, const char *ibuff)
{
   const char *tptr = ibuff;
   int i = 0;

   while( tptr && i <= cptr->n_psv_fields)
      {
      tptr = strchr( tptr, '|');
      if( tptr)
         tptr++;
      i++;
      }
   if( i != cptr->n_psv_fields)
      return( 0);       /* it's not a PSV line */
   setup_observation( cptr);
   for( i = 0; i < cptr->n_psv_fields; i++)
      {
      process_psv_tag( cptr, obuff, cptr->psv_tags[i], ibuff);
      ibuff = strchr( ibuff, '|');
      if( i != cptr->n_psv_fields - 1)
         {
         assert( ibuff);
         ibuff++;
         }
      }
   assert( !ibuff);
   return( 1);
}

#define MIN_PSV_TAGS    4

/* Returns either zero if the line is not a PSV header,  or the number
of PSV fields (which equals the number of pipe separators plus one.) */

static int check_for_psv_header( ades2mpc_t *cptr, const char *buff)
{
   int n_psv_tags = 0;
   const char *tptr = buff;

   while( tptr)
      {
      tptr = skip_whitespace( tptr);
      if( *tptr >= 'a' && *tptr <= 'z')
         {                /* could be a PSV header field */
         n_psv_tags++;
         tptr = strchr( tptr, '|');
         if( tptr)
            tptr++;
         }
      else     /* PSV header fields must start with a lowercase letter */
         return( 0);
      }
   if( n_psv_tags < MIN_PSV_TAGS)
      return( 0);
   cptr->n_psv_fields = n_psv_tags;
   if( cptr->psv_tags)
      free( cptr->psv_tags);
   cptr->psv_tags = (int *)calloc( n_psv_tags, sizeof( int));
   tptr = buff;
   n_psv_tags = 0;
   while( tptr)
      {
      int i = 0;

      tptr = skip_whitespace( tptr);
      while( tptr[i] != '|' && tptr[i] > ' ')
         i++;
      cptr->psv_tags[n_psv_tags] = find_tag( tptr, i);
      assert( cptr->psv_tags[n_psv_tags] > 0);
      n_psv_tags++;
      assert( tptr);
      tptr = strchr( tptr, '|');
      if( tptr)
         tptr++;
      }
   assert( !tptr);
   return( n_psv_tags);
}

static int process_psv_header( ades2mpc_t *cptr, char *obuff, const char *ibuff)
{
   size_t i = 2;
   int itag, rval = 0;

   while( ibuff[i] > ' ')
      i++;
   itag = find_tag( ibuff + 2, i - 2);
   if( itag >= 0)
      {
      if( *ibuff == '#')
         {
         cptr->tags[1] = itag;
         rval = 2;         /* signal non-visible 'container' tag */
         }
      else
         {
         cptr->tags[2] = itag;
         cptr->depth = 3;
         ibuff = skip_whitespace( ibuff + i);
         i = 0;
         while( ibuff[i] >= ' ')
            i++;
         process_ades_tag( obuff, cptr, itag, ibuff, i);
         cptr->depth = 0;
         rval = 1;
         }
      }
   return( rval);
}

#define ADES_CLOSING_UNOPENED_TAG         (-2)
#define ADES_DEPTH_MAX                    (-3)
#define ADES_MALFORMED_TAG                (-4)

int xlate_ades2mpc( void *context, char *obuff, const char *buff)
{
   int rval = 0;
   const char *tptr;
   ades2mpc_t *cptr = (ades2mpc_t *)context;
   char temp_obuff[300], *orig_obuff = NULL;

   if( cptr->prev_line_passed_through)
      {
      cptr->prev_line_passed_through = 0;
      return( 0);
      }
   if( cptr->getting_lines)
      return( get_a_line( obuff, cptr));
   if( !cptr->depth && strstr( buff, "<optical>"))
      cptr->depth = 1;
   if( check_for_psv_header( cptr, buff))
      {
      cptr->depth = 1;
      return( 0);
      }
   if( cptr->psv_tags)
      {
      rval = process_psv_line( cptr, (obuff == buff ? temp_obuff : obuff), buff);
      if( !rval)        /* we've reached the end of a PSV data section */
         {
         free( cptr->psv_tags);
         cptr->psv_tags = NULL;
         cptr->n_psv_fields = 0;
         cptr->depth = 0;
         }
      }
   if( !rval)
      rval = check_for_tholen_sigmas( cptr, (obuff == buff ? temp_obuff : obuff), buff);
   if( !rval && buff[1] == ' ' && (buff[0] == '#' || buff[0] == '!'))
      {
      rval = process_psv_header( cptr, (obuff == buff ? temp_obuff : obuff), buff);
      if( rval)
         {
         if( rval == 1)
            {
            if( obuff == buff)
               strcpy( obuff, temp_obuff);
            }
         else        /* 'container' tag;  affects internal state,  but */
            rval = 0;    /* nothing is output */
         cptr->prev_line_passed_through = rval;
         return( rval);
         }
      }
   if( rval)
      cptr->getting_lines = 1;
   if( !rval && !cptr->depth && !strstr( buff, "<ades version"))
      {
      if( obuff != buff)
         strcpy( obuff, buff);
      cptr->prev_line_passed_through = 1;
      return( 1);
      }
   if( cptr->getting_lines)
      return( get_a_line( obuff, cptr));
   if( obuff == buff)            /* translating in place */
      {
      orig_obuff = obuff;
      obuff = temp_obuff;
      }
   tptr = skip_whitespace( buff);
   while( rval >= 0 && *tptr)
      {
      size_t len = 0;

      if( *tptr == '<')        /* start of tag,  we assume */
         {
         while( tptr[len] && tptr[len] != '>')
            len++;
         if( tptr[len] != '>')
            return( ADES_MALFORMED_TAG);
         else
            {
            const int tag_idx = find_tag( tptr + 1, len - 1);

            if( tag_idx >= 0)
               {
               if( tptr[1] == '/')
                  {
                  cptr->depth--;
                  if( cptr->depth < 0 || cptr->tags[cptr->depth] != tag_idx)
                     rval = ADES_CLOSING_UNOPENED_TAG;
                  }
               else
                  {
                  cptr->tags[cptr->depth++] = tag_idx;
                  if( cptr->depth == MAX_DEPTH)
                     rval = ADES_DEPTH_MAX;
                  }
               if( tag_idx == ADES_optical)
                  {
                  if( tptr[1] == '/')
                     cptr->getting_lines = rval = 1;
                  else
                     setup_observation( cptr);
                  }
               }
            }
         tptr += len + 1;
         }
      else if( cptr->depth)
         {
         const int itag = cptr->tags[cptr->depth - 1];

         while( tptr[len] && tptr[len] != '<')
            len++;
         while( len && (unsigned char)tptr[len - 1] <= ' ')
            len--;
         rval = process_ades_tag( obuff, cptr, itag, tptr, len);
         tptr += len;
         }
      tptr = skip_whitespace( tptr);
      }
   if( rval)
      {
      if( cptr->line2[0])     /* yes,  we've a valid line 2 */
         {
         memcpy( cptr->line2, cptr->line, 12);
         memcpy( cptr->line2 + 15, cptr->line + 15, 17);
         memcpy( cptr->line2 + 77, cptr->line + 77, 3);
         }
      get_a_line( obuff, cptr);
      }
   if( orig_obuff)
      strcpy( orig_obuff, temp_obuff);
   return( rval);
}

int xlate_ades2mpc_in_place( void *context, char *buff)
{
   const int rval = xlate_ades2mpc( context, buff, buff);

   return( rval);
}

int fgets_with_ades_xlation( char *buff, const size_t len,
                        void *ades_context, FILE *ifile)
{
   ades2mpc_t *cptr = (ades2mpc_t *)ades_context;
   int prev_rval = cptr->prev_rval;

   if( prev_rval)
      prev_rval = xlate_ades2mpc_in_place( ades_context, buff);
   while( !prev_rval && fgets( buff, (int)len, ifile))
      prev_rval = xlate_ades2mpc_in_place( ades_context, buff);
   while( *buff && *buff != 10 && *buff != 13)
      buff++;
   *buff = '\0';
   cptr->prev_rval = prev_rval;
   return( prev_rval);
}
