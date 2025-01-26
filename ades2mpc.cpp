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
02110-1301, USA.

   The following functions support reading of PSV and XML ADES
data.  These can be mixed with one another,  and with "traditional"
80-column astrometry.  It is the underlying code for all astrometry
parsing in my software.

   See 'adestest.cpp' for a simple example of its use.   */

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "stringex.h"
#include <stdio.h>
#include <ctype.h>
#include "watdefs.h"
#include "date.h"
#include "mpc_func.h"
#include "afuncs.h"

#define NOT_A_VALID_TIME -3.141e+17
#define MAX_DEPTH 20
#define PIECE_SIZE   26

typedef struct
{
   int depth, tags[MAX_DEPTH];
   char line[83];    /* allow possible CR, LF,  & null */
   char line2[83];
   char rms_ra[PIECE_SIZE], rms_dec[PIECE_SIZE], corr[PIECE_SIZE];
   char rms_mag[PIECE_SIZE], rms_time[PIECE_SIZE], unc_time[PIECE_SIZE];
   char full_ra[PIECE_SIZE], full_dec[PIECE_SIZE], full_mag[9];
   char notes[7], program_code[3];
   char trk_sub[14], obs_id[PIECE_SIZE], trk_id[12], passband[4];
   long double full_t2k;
   int id_set, getting_lines, spacecraft_center;
   int prev_line_passed_through;
   int prev_rval, n_psv_fields;
   int *psv_tags;
   double spacecraft_vel[3];
   bool ignore_artsat_desigs;
} ades2mpc_t;

void *init_ades2mpc( void)
{
   ades2mpc_t *rval = (ades2mpc_t *)calloc( 1, sizeof( ades2mpc_t));

   rval->full_t2k = NOT_A_VALID_TIME;
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
#define ADES_fltr                         44
#define ADES_frq                          45
#define ADES_fundingSource                46
#define ADES_institution                  47
#define ADES_line                         48
#define ADES_localUse                     49
#define ADES_logSNR                       50
#define ADES_mag                          51
#define ADES_measurers                    52
#define ADES_mode                         53
#define ADES_mpcCode                      54
#define ADES_nStars                       55
#define ADES_name                         56
#define ADES_notes                        57
#define ADES_nucMag                       58
#define ADES_objectDetection              59
#define ADES_obsBlock                     60
#define ADES_obsCenter                    61
#define ADES_obsContext                   62
#define ADES_obsData                      63
#define ADES_obsID                        64
#define ADES_obsSubID                     65
#define ADES_obsTime                      66
#define ADES_observatory                  67
#define ADES_observers                    68
#define ADES_occultation                  69
#define ADES_offset                       70
#define ADES_optical                      71
#define ADES_opticalResidual              72
#define ADES_orbID                        73
#define ADES_orbProd                      74
#define ADES_pa                           75
#define ADES_permID                       76
#define ADES_photAp                       77
#define ADES_photCat                      78
#define ADES_photMod                      79
#define ADES_photProd                     80
#define ADES_photometry                   81
#define ADES_pixelScale                   82
#define ADES_pos1                         83
#define ADES_pos2                         84
#define ADES_pos3                         85
#define ADES_posCov11                     86
#define ADES_posCov12                     87
#define ADES_posCov13                     88
#define ADES_posCov22                     89
#define ADES_posCov23                     90
#define ADES_posCov33                     91
#define ADES_precDec                      92
#define ADES_precRA                       93
#define ADES_precTime                     94
#define ADES_prog                         95
#define ADES_provID                       96
#define ADES_ra                           97
#define ADES_raStar                       98
#define ADES_radar                        99
#define ADES_radarResidual               100
#define ADES_rcv                         101
#define ADES_ref                         102
#define ADES_remarks                     103
#define ADES_resDec                      104
#define ADES_resDelay                    105
#define ADES_resDoppler                  106
#define ADES_resMag                      107
#define ADES_resRA                       108
#define ADES_rmsCorr                     109
#define ADES_rmsDec                      110
#define ADES_rmsDelay                    111
#define ADES_rmsDist                     112
#define ADES_rmsDoppler                  113
#define ADES_rmsFit                      114
#define ADES_rmsMag                      115
#define ADES_rmsPA                       116
#define ADES_rmsRA                       117
#define ADES_rmsTime                     118
#define ADES_seeing                      119
#define ADES_selAst                      120
#define ADES_selDelay                    121
#define ADES_selDoppler                  122
#define ADES_selPhot                     123
#define ADES_shapeOcc                    124
#define ADES_sigCorr                     125
#define ADES_sigDec                      126
#define ADES_sigDelay                    127
#define ADES_sigDoppler                  128
#define ADES_sigMag                      129
#define ADES_sigRA                       130
#define ADES_sigTime                     131
#define ADES_software                    132
#define ADES_stn                         133
#define ADES_subFmt                      134
#define ADES_subFrm                      135
#define ADES_submitter                   136
#define ADES_sys                         137
#define ADES_telescope                   138
#define ADES_trkID                       139
#define ADES_trkMPC                      140
#define ADES_trkSub                      141
#define ADES_trx                         142
#define ADES_uncTime                     143
#define ADES_vel1                        144
#define ADES_vel2                        145
#define ADES_vel3                        146

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
       "fitOrder", "fltr", "frq", "fundingSource", "institution",
       "line", "localUse", "logSNR", "mag", "measurers", "mode",
       "mpcCode", "nStars", "name", "notes", "nucMag",
       "objectDetection", "obsBlock", "obsCenter", "obsContext",
       "obsData", "obsID", "obsSubID", "obsTime", "observatory",
       "observers", "occultation", "offset", "optical",
       "opticalResidual", "orbID", "orbProd", "pa", "permID",
       "photAp", "photCat", "photMod", "photProd", "photometry",
       "pixelScale", "pos1", "pos2", "pos3", "posCov11",
       "posCov12", "posCov13", "posCov22", "posCov23", "posCov33",
       "precDec", "precRA", "precTime", "prog", "provID", "ra",
       "raStar", "radar", "radarResidual", "rcv", "ref",
       "remarks", "resDec", "resDelay", "resDoppler", "resMag",
       "resRA", "rmsCorr", "rmsDec", "rmsDelay", "rmsDist",
       "rmsDoppler", "rmsFit", "rmsMag", "rmsPA", "rmsRA",
       "rmsTime", "seeing", "selAst", "selDelay", "selDoppler",
       "selPhot", "shapeOcc", "sigCorr", "sigDec", "sigDelay",
       "sigDoppler", "sigMag", "sigRA", "sigTime", "software",
       "stn", "subFmt", "subFrm", "submitter", "sys", "telescope",
       "trkID", "trkMPC", "trkSub", "trx", "uncTime", "vel1",
       "vel2", "vel3",
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

/* Packs references such as 'MPEC 2021-D85' into the five-byte form
'ED085'.  */

static inline void pack_mpc_reference( char *packed, const char *ref)
{
   const size_t len = strlen( ref);
   char temp_packed[6];

   if( len >= 12 && len <= 14 && !memcmp( ref, "MPEC ", 5) && ref[9] == '-')
      {
      *packed++ = 'E';
      *packed++ = ref[10];
      ref += 11;
      packed[0] = packed[1] = '0';
      memcpy( packed + 14 - len, ref, len - 11);
      }
   else if( !memcmp( ref, "MPS ", 4))
      {
      unsigned mps_number = atoi( ref + 4);

      if( mps_number < 260000)
         {
         snprintf_err( temp_packed, sizeof( temp_packed),
                                "%c%04u", 'a' + mps_number / 10000,
                                mps_number % 10000);
         memcpy( packed, temp_packed, 5);
         }
      else
         {
         *packed = '~';
         encode_value_in_mutant_hex( packed + 1, 4, mps_number - 260000);
         }
      }
   else if( !memcmp( ref, "MPC ", 4))
      {
      unsigned mpc_number = atoi( ref + 4);

      if( mpc_number < 110000)
         {
         snprintf_err( temp_packed, sizeof( temp_packed),
               (mpc_number < 100000 ? "%05u" : "@%04u"), mpc_number % 100000);
         memcpy( packed, temp_packed, 5);
         }
      else
         {
         *packed = '#';
         encode_value_in_mutant_hex( packed + 1, 4, mpc_number - 110000);
         }
      }
   else if( *ref == '!')
      {
      memcpy( packed, "     ", 5);
      memcpy( packed, ref, len > 5 ? 5 : len);
      }
}

static inline size_t move_fits_time( char *optr, const char *iptr)
{
   char *optr0 = optr;

   if( iptr[4] == '-' && iptr[7] == '-' && iptr[10] == 'T'
               && iptr[13] == ':')
      {
      const int century = iptr[0] * 10 + iptr[1] - '0' * 11;

      *optr++ = int_to_mutant_hex_char( century);
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
   return( optr - optr0);
}

static inline int place_value( char *optr, const char *iptr, const size_t ilen,
             size_t leading_places)
{
   size_t i, point_loc = ilen;

   for( i = 0; i < ilen; i++)
      if( iptr[i] == '.')
         point_loc = i;
   assert( leading_places >= point_loc);
   if( ilen - point_loc > 8)     /* overly long;  read value and round it */
      {
      char tbuff[13];

      snprintf_err( tbuff, sizeof( tbuff), (leading_places == 2 ? "%11.8f " : "%12.8f"),
                        atof( iptr));
      memcpy( optr, tbuff, 12);
      return (1);
      }
   while( leading_places > point_loc)
      {
      *optr++ = '0';
      leading_places--;
      }
   memcpy( optr, iptr, (ilen < 9 + leading_places) ? ilen : 9 + leading_places);
   return( 0);
}

static const char *skip_whitespace( const char *tptr)
{
   while( *tptr && isspace( *tptr))
      tptr++;
   return( tptr);
}

static int get_a_line( char *obuff, const size_t obuff_size, ades2mpc_t *cptr)
{
   if (cptr->rms_ra[0])
      {
      snprintf_err( obuff, obuff_size, "COM Sigmas %s", cptr->rms_ra);
      if( strcmp( cptr->rms_ra, cptr->rms_dec))
         {
         snprintf_append( obuff, obuff_size, "x%s", cptr->rms_dec);
         if( atof( cptr->corr))
            snprintf_append( obuff, obuff_size, ",%s", cptr->corr);
         }
      if( cptr->rms_mag[0])
         {
         strlcat_err( obuff, " m:", obuff_size);
         strlcat_err( obuff, cptr->rms_mag, obuff_size);
         cptr->rms_mag[0] = '\0';
         }
      if( cptr->rms_time[0])
         {
         strlcat_err( obuff, " t:", obuff_size);
         strlcat_err( obuff, cptr->rms_time, obuff_size);
         cptr->rms_time[0] = '\0';
         }
      if( cptr->unc_time[0])
         {
         strlcat_err( obuff, " u:", obuff_size);
         strlcat_err( obuff, cptr->unc_time, obuff_size);
         cptr->unc_time[0] = '\0';
         }
      cptr->rms_ra[0] = '\0';
      strlcat_err( obuff, "\n", obuff_size);
      }
   else if( cptr->trk_sub[0] || cptr->obs_id[0] || cptr->trk_id[0])
      {
      strlcpy_err( obuff, "COM IDs", obuff_size);
      if( cptr->trk_sub[0])
         snprintf_append( obuff, obuff_size, " trkSub:%s", cptr->trk_sub);
      if( cptr->obs_id[0])
         snprintf_append( obuff, obuff_size, " obsID:%s", cptr->obs_id);
      if( cptr->trk_id[0])
         snprintf_append( obuff, obuff_size, " trkID:%s", cptr->trk_id);
      cptr->trk_sub[0] = cptr->obs_id[0] = cptr->trk_id[0] = '\0';
      strlcat_err( obuff, "\n", obuff_size);
      }
   else if( cptr->full_ra[0] || cptr->full_dec[0] || cptr->full_t2k != NOT_A_VALID_TIME)
      {
      snprintf_err( obuff, obuff_size, "COM RA/dec %s %s",
               (cptr->full_ra[0] ? cptr->full_ra : "-"),
               (cptr->full_dec[0] ? cptr->full_dec : "-"));
      if( cptr->full_t2k != NOT_A_VALID_TIME)
         snprintf_append( obuff, obuff_size, " %.15f", (double)cptr->full_t2k);
      strlcat_err( obuff, "\n", obuff_size);
      cptr->full_ra[0] = cptr->full_dec[0] = '\0';
      cptr->full_t2k = NOT_A_VALID_TIME;
      }
   else if( cptr->full_dec[0])
      {
      snprintf_err( obuff, obuff_size, "COM RA/dec - %s\n", cptr->full_dec);
      cptr->full_dec[0] = '\0';
      }
   else if( cptr->full_mag[0])
      {
      snprintf_err( obuff, obuff_size, "COM full mag %s\n", cptr->full_mag);
      cptr->full_mag[0] = '\0';
      }
   else if( cptr->passband[1] || cptr->notes[0] > ' ' || cptr->program_code[0] >= ' ')
      {
      snprintf_err( obuff, obuff_size, "COM ADES tags");
      if( cptr->passband[1])
         snprintf_append( obuff, obuff_size, " band:%s", cptr->passband);
      if( cptr->notes[0] >= ' ')
         snprintf_append( obuff, obuff_size, " notes:%s", cptr->notes);
      if( cptr->program_code[0] >= ' ')
         snprintf_append( obuff, obuff_size, " progcode:%s", cptr->program_code);
      cptr->passband[1] = cptr->notes[0] = cptr->program_code[0] = '\0';
      }
   else if( cptr->spacecraft_vel[0] || cptr->spacecraft_vel[1] || cptr->spacecraft_vel[2])
      {
      double multiplier = 1.;

      if( cptr->line2[32] == '2')   /* units are AU and days;  */
         multiplier = AU_IN_KM / seconds_per_day;
      snprintf_err( obuff, obuff_size, "COM vel (km/s) %.14s  %+13.7f%13.7f%13.7f %s",
                     cptr->line + 15,
                     cptr->spacecraft_vel[0] * multiplier,
                     cptr->spacecraft_vel[1] * multiplier,
                     cptr->spacecraft_vel[2] * multiplier,
                     cptr->line + 77);
      memset( cptr->spacecraft_vel, 0, sizeof( cptr->spacecraft_vel));
      }
   else if( cptr->line[0])
      {
      strlcpy( obuff, cptr->line, obuff_size);
      if( cptr->line2[0])     /* yes,  we've a valid line 2 */
         {
         memcpy( cptr->line2, cptr->line, 12);
         memcpy( cptr->line2 + 15, cptr->line + 15, 17);
         if( cptr->spacecraft_center != 399)
            snprintf_err( cptr->line2 + 69, 9, "%8d", cptr->spacecraft_center);
         memcpy( cptr->line2 + 77, cptr->line + 77, 3);
         }
      cptr->line[0] = '\0';
      }
   else if( cptr->line2[0])
      {
      strlcpy_err( obuff, cptr->line2, obuff_size);
      cptr->line2[0] = '\0';
      }
   else                         /* got all the lines we'll get */
      cptr->getting_lines = 0;
   return( cptr->getting_lines);
}

/* Returns 1 if it's a properly handled header tag,  0 if it's some
other tag or among the remaining unhandled header tags (I'm not
dealing with the telescope details yet,  for example.) */

static int process_ades_tag( char *obuff, ades2mpc_t *cptr, const int itag,
                 const char *tptr, size_t len)
{
   int rval = 0;
   char name[40];
   const size_t obuff_size = 221;

   assert( obuff);
   *obuff = '\0';
   if( len < sizeof( name))
      {
      memcpy( name, tptr, len);
      name[len] = '\0';
      }
   switch( itag)
      {
      case ADES_mpcCode:
         snprintf_err( obuff, obuff_size, "COD %s\n", name);
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
            snprintf_err( obuff, obuff_size, format, (int)len, tptr);
         else
            strlcpy_err( obuff, "COM Mangled name data\n", obuff_size);
         rval = 1;
         }
         break;
      case ADES_stn:
         memcpy( cptr->line + 77, tptr, 3);
         break;
      case ADES_obsTime:
         {
         const bool too_far_in_future = (atoi( tptr) > 2099);

         if( move_fits_time( cptr->line + 15, tptr) > 17 || too_far_in_future)
            {
            char *zptr = strchr( name, 'Z');

            if( zptr)
               *zptr = '\0';
            cptr->full_t2k = get_time_from_stringl( 0., name, 0, NULL);
            if( too_far_in_future)
               cptr->line[15] = 'K';
            }
         }
         break;
      case ADES_band:
         if( (*tptr == 'P' || *tptr == 'S')
                           && strchr( "grizwy", tptr[1]) && tptr[1])
            cptr->line[70] = tptr[1];     /* PanSTARRS or Sloan band */
         else if( *tptr == 'A' && (tptr[1] == 'o' || tptr[1] == 'c'))
            cptr->line[70] = tptr[1];   /* ATLAS Ao & Ac bands */
         else if( *tptr == 'G' && (tptr[1] == 'b' || tptr[1] == 'r'))
            cptr->line[70] = tptr[1];   /* Gaia b or r band */
         else
            cptr->line[70] = *tptr;
         assert( len > 0 && len < 4);    /* three-byte passcodes are allowed */
         strcpy( cptr->passband, name);
         break;
      case ADES_mode:
         if( len == 3)
            {     /* https://www.minorplanetcenter.net/iau/info/ADESFieldValues.html */
            const char *modes =  "CCCD BCMO nVID PPHO eENC pPMT"
                                " MMIC TMER CTDI EOCC ?UNK ";
            int i;

            for( i = 0; modes[i]; i += 5)
               if( !memcmp( modes + i + 1, tptr, 3))
                  cptr->line[14] = modes[i];
            }
         break;
      case ADES_deprecated:
         cptr->line[14] = 'X';
         break;
      case ADES_disc:
         if( *tptr == '*')
            cptr->line[12] = '*';
         break;
      case ADES_ref:
         if( len < sizeof( name))
            pack_mpc_reference( cptr->line + 72, name);
         break;
      case ADES_prog:
         {
         const char *programs = "0123456789!\"#$%&'()*+,-./[\\]^_`{|}~";
         const int idx = atoi( tptr);

         if( idx >=0 && idx <= 34)
            cptr->line[13] = programs[idx];
         assert( len < sizeof( cptr->program_code));
         strlcpy_err( cptr->program_code, name, sizeof( cptr->program_code));
         }
         break;
      case ADES_sys:
         cptr->line2[0] = ' ';
         cptr->line2[14] = 's';    /* assume spacecraft-based */
         if( !strcmp( name, "ICRF_KM"))
            cptr->line2[32] = '1';
         else if( !strcmp( name, "ICRF_AU"))
            cptr->line2[32] = '2';
         else if( !strcmp( name, "WGS84"))
            {
            cptr->line2[32] = '1';
            cptr->line2[14] = 'v';     /* nope,  it's a roving observer */
            }
         else if( !strcmp( name, "ITRF") || !strcmp( name, "IAU"))
            {
            strlcpy_err( obuff, "Can't handle <sys> ITRF or IAU yet\n", obuff_size);
            assert( 0);
            }
         else
            {
            strlcpy_err( obuff, "Bad <sys> tag\n", obuff_size);
            assert( 0);
            }
         cptr->line[14]  = toupper( cptr->line2[14]);
         break;
      case ADES_ctr:
         cptr->spacecraft_center = atoi( name);
         break;
      case ADES_vel1:
      case ADES_vel2:
      case ADES_vel3:
         assert( cptr->line[14] == 'S');
         cptr->spacecraft_vel[itag - ADES_vel1] = atof( name);
         break;
      case ADES_pos1:
      case ADES_pos2:
      case ADES_pos3:
         assert( cptr->line[14] == 'S' || cptr->line[14] == 'V');
         if( cptr->line[14] == 'S')    /* satellite offset */
            {
            const int sign_loc = 34 + (itag - ADES_pos1) * 12;
            int nlen = (int)len, decimal_loc = 0;
            char *tptr2;

            if( cptr->line2[32] == '2')
               decimal_loc = sign_loc + 6;
            else if( cptr->line2[32] == '1')
               decimal_loc = sign_loc + 2;
            else
               {
               strlcpy_err( obuff, "Bad posn data\n", obuff_size);
               rval = 1;
               }
                     /* cvt scientific notation,  if any : */
            if( strchr( name, 'e') || strchr( name, 'E'))
               {
               snprintf_err( name, sizeof( name), "%.13f", atof( name));
               nlen = 12;
               }
            if( *name != '+' && *name != '-')   /* no sign provided; */
               {                             /* insert one */
               nlen++;
               memmove( name + 1, name, nlen);
               *name = '+';
               }
            cptr->line2[sign_loc] = *name;
            tptr2 = strchr( name, '.');
            assert( tptr2);
            if( cptr->line2[32] == '1')
               {
               decimal_loc = sign_loc + 6;
               if( tptr2 - name >= 7)     /* beyond 100000 km */
                  decimal_loc++;
               if( tptr2 - name >= 8)     /* one to ten million km */
                  decimal_loc++;
               assert( tptr2 - name < 9);
               }
            else if( cptr->line2[32] == '2')
               {
               decimal_loc = sign_loc + 2;
               if( tptr2 - name == 3)     /* beyond 10 AU */
                  decimal_loc++;
               assert( tptr2 - name < 4);
               }
            else
               {
               strlcpy_err( obuff, "Bad posn data\n", obuff_size);
               rval = 1;
               }
            if( decimal_loc)
               {
               decimal_loc -= (int)(tptr2 - name);
               memcpy( &cptr->line2[decimal_loc + 1], name + 1,
                                           (nlen > 11 ? 10 : nlen - 1));
               }
            }
         else               /* roving observer */
            {
            const double ival = atof( name);
            size_t loc;

            if( itag == ADES_pos1)     /* East longitude */
               snprintf_err( cptr->line2 + 34, 10, "%9.5f",
                              fmod( ival + 360., 360.));
            else if( itag == ADES_pos2)      /* latitude must be signed */
               snprintf_err( cptr->line2 + 45, 10, "%+9.5f", ival);
            else           /* altitude */
               snprintf_err( cptr->line2 + 56, 6, "%5d", (int)ival);
            loc = strlen( cptr->line2);
            cptr->line2[loc] = ' ';
            }
         break;
      case ADES_ra:
         if( *tptr == '+')
            {
            tptr++;
            len--;
            }
         if( place_value( cptr->line + 32, tptr, len, 3))
            {
            memcpy( cptr->full_ra, tptr, len);     /* 'overlong' RA */
            cptr->full_ra[len] = '\0';
            }
         break;
      case ADES_dec:
         if( *tptr == '-' || *tptr == '+')
            {
            cptr->line[44] = *tptr++;
            len--;
            }
         else
            cptr->line[44] = '+';
         if( place_value( cptr->line + 45, tptr, len, 2))
            {              /* 'overlong' dec */
            cptr->full_dec[0] = cptr->line[44];
            memcpy( cptr->full_dec + 1, tptr, len);
            cptr->full_dec[len + 1] = '\0';
            }
         break;
      case ADES_astCat:
         assert( len < sizeof( name));
         cptr->line[71] = net_name_to_byte_code( name);
         assert( cptr->line[71]);
         break;
      case ADES_rmsRA:
         assert( len < sizeof( cptr->rms_ra));
         strlcpy_err( cptr->rms_ra, name, sizeof( cptr->rms_ra));
         break;
      case ADES_notes:
         assert( len < sizeof( cptr->notes));
         strlcpy_err( cptr->notes, name, sizeof( cptr->notes));
         cptr->line[13] = name[0];
         break;
      case ADES_rmsDec:
         assert( len < sizeof( cptr->rms_dec));
         strlcpy_err( cptr->rms_dec, name, sizeof( cptr->rms_dec));
         break;
      case ADES_rmsCorr:
         assert( len < sizeof( cptr->corr));
         strlcpy_err( cptr->corr, name, sizeof( cptr->corr));
         break;
      case ADES_rmsTime:
         assert( len < sizeof( cptr->rms_time));
         strlcpy_err( cptr->rms_time, name, sizeof( cptr->rms_time));
         break;
      case ADES_uncTime:
         assert( len < sizeof( cptr->unc_time));
         strlcpy_err( cptr->unc_time, name, sizeof( cptr->unc_time));
         break;
      case ADES_rmsMag:
         assert( len < sizeof( cptr->rms_mag));
         strlcpy_err( cptr->rms_mag, name, sizeof( cptr->rms_mag));
         break;
      case ADES_provID:
         if( cptr->id_set == ADES_permID)
            break;                  /* FALLTHRU */
      case ADES_permID:
         {
         char tbuff[20];
         int i = 0;

         assert( len < sizeof( name));
         while( isdigit( name[i]))
            i++;
         if( !name[i])        /* simple numbered object,  desig is nothing but */
            {                 /* digits : put parentheses around it */
            memmove( name + 1, name, i);
            name[0] = '(';
            name[i + 1] = ')';
            name[i + 2] = '\0';
            }
         create_mpc_packed_desig( tbuff, name);
         memcpy( cptr->line, tbuff, 12);
         cptr->id_set = itag;
         }
         break;
      case ADES_artSat:
         assert( len < 13);
         if( !cptr->id_set && !cptr->ignore_artsat_desigs)
            {
            cptr->id_set = ADES_artSat;
            memcpy( cptr->line, tptr, len);
            }
         break;
      case ADES_trkSub:
         assert( len < 13);
         strlcpy_err( cptr->trk_sub, name, sizeof( cptr->trk_sub));
         if( !cptr->id_set)
            {
            cptr->id_set = ADES_trkSub;
            if( len < 8)
               memcpy( cptr->line + 5, tptr, len);
            else
               memcpy( cptr->line + 12 - len, tptr, len);
            }
         break;
         break;
      case ADES_trkID:
         assert( len < sizeof( cptr->trk_id));
         strlcpy_err( cptr->trk_id, name, sizeof( cptr->trk_id));
         break;
      case ADES_obsID:
         assert( len < sizeof( cptr->obs_id));
         strlcpy_err( cptr->obs_id, name, sizeof( cptr->obs_id));
         break;
      case ADES_mag:
         memcpy( cptr->line + 65, tptr, (len < 5) ? len : 5);
         if( len > 5)
            strlcpy_err( cptr->full_mag, name, sizeof( cptr->full_mag));
         break;
      case ADES_trkMPC:    /* trkSub if it's deprecated;  */
         break;            /* currently ignored           */
      default:
         snprintf_err( obuff, obuff_size, "COM Unhandled ADES tag %d", itag);
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
   strlcpy_err( cptr->line + 80, "\n", 2);
   memset( cptr->line2, ' ', 80);
   strlcpy_err( cptr->line2 + 80, "\n", 2);
   cptr->line2[0] = '\0';
   cptr->id_set = 0;
   cptr->spacecraft_center = 399;      /* default to geocentric */
   cptr->full_t2k = NOT_A_VALID_TIME;
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
         strlcpy_err( cptr->line, obuff, sizeof( cptr->line));
         obuff[12] = '\0';
         if( !create_mpc_packed_desig( packed_desig, obuff))
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
/*    if( cptr->psv_tags[n_psv_tags] <= 0)
         fprintf( stderr, "Tag %d '%s' failed in '%s'\n", n_psv_tags, tptr, buff); */
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
         rval = process_ades_tag( obuff, cptr, itag, ibuff, i);
         cptr->depth = 0;
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
   const size_t obuff_size = 221;

   if( cptr->prev_line_passed_through)
      {
      cptr->prev_line_passed_through = 0;
      return( 0);
      }
   if( cptr->getting_lines)
      return( get_a_line( obuff, obuff_size, cptr));
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
               strlcpy_err( obuff, temp_obuff, obuff_size);
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
         strlcpy_err( obuff, buff, obuff_size);
      cptr->prev_line_passed_through = 1;
      return( 1);
      }
   if( cptr->getting_lines)
      return( get_a_line( obuff, obuff_size, cptr));
   if( obuff == buff)            /* translating in place */
      {
      orig_obuff = obuff;
      obuff = temp_obuff;
      *temp_obuff = '\0';
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
      else
         return rval;
      tptr = skip_whitespace( tptr);
      }
   if( rval)
      get_a_line( obuff, obuff_size, cptr);
   if( orig_obuff)
      strlcpy_err( orig_obuff, temp_obuff, obuff_size);
   return( rval);
}

int xlate_ades2mpc_in_place( void *context, char *buff)
{
   const int rval = xlate_ades2mpc( context, buff, buff);

   return( rval);
}

/* 'fgets' with trailing LF, CR/LF,  or CR trimmed.  For a CR terminator,
we may actually read in some extra bytes and therefore have to fseek()
backward to the point right after the CR in question... fortunately,
such files are quite rare as of 2020 (they were more common a couple of
decades earlier.) */

static char *fgets_trimmed( char *buff, const int len, FILE *ifile)
{
   char *rval = fgets( buff, len, ifile);

   if( rval)
      {
      int i = 0;

      while( buff[i] != 10 && buff[i] != 13 && buff[i])
         i++;
      if( buff[i] == 13 && buff[i + 1] != 10)   /* CR terminated */
         fseek( ifile, 1L + (long)i - (long)strlen( buff), SEEK_CUR);
      while( i && buff[i - 1] == ' ')
         i--;           /* drop trailing spaces */
      buff[i] = '\0';
      }
   return( rval);
}

int fgets_with_ades_xlation( char *buff, const size_t len,
                        void *ades_context, FILE *ifile)
{
   ades2mpc_t *cptr = (ades2mpc_t *)ades_context;
   int prev_rval = cptr->prev_rval;

   if( prev_rval)
      prev_rval = xlate_ades2mpc_in_place( ades_context, buff);
   while( !prev_rval && fgets_trimmed( buff, (int)len, ifile))
      prev_rval = xlate_ades2mpc_in_place( ades_context, buff);
   while( *buff && *buff != 10 && *buff != 13)
      buff++;
   *buff = '\0';
   cptr->prev_rval = prev_rval;
   return( prev_rval);
}

void ades_artsat_desigs( void *ades_context, const bool ignore_artsat_desigs)
{
   ades2mpc_t *cptr = (ades2mpc_t *)ades_context;

   cptr->ignore_artsat_desigs = ignore_artsat_desigs;
}
