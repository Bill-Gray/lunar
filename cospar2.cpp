/* cospar2.cpp: computes planet/satellite orientations
(STRONLY ADVISE USING 'COSPAR.CPP" INSTEAD.  This version hasn't
been updated for a while and has various weaknesses compared to the
approach taken in 'cospar.cpp'.)

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

#include <math.h>
#include <string.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"

/* Older version of the COSPAR code to compute orientations of planets
and natural satellites.  The newer version is in 'cospar.cpp'.  This
version is being kept around for reference,  and in case somebody has
been using it and needs a roadmap linking the new 'cospar.cpp' (which
gets most of its data from the file 'cospar.txt') to the older
implementation (which had a slew of coefficients built into it). */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define J2000 2451545.0
#define HALF_PI (PI / 2.)
#define TRIG_TERM struct trig_term
#define N_TRIG_TERMS 14

TRIG_TERM
   {
   char planet_no;
   unsigned short arg;
   long omega;
   short ra_term, dec_term;
   long w_term;
   };

/* 15 Nov 98:  Added some data for Saturn systems I and II,  based on  */
/* software from Jeff Beish.  This is described in BEISH.DOC.          */

static const TRIG_TERM tterms[N_TRIG_TERMS] = {
   {  7, 35785,     52316L,    700,  -510,  -480} ,           /* Neptune */
   { 11, 28390,   4850700L,     94,    40,    85} ,           /* Io */
   { 15, 17740, -36505500L,  13560, -1530,-13480} ,           /* Mimas */
   { 15, 31645,    506200L,      0,     0,-44850} ,           /* Mimas */
   { 17, 31645,    506200L,      0,     0,  2230} ,           /* Tethys */
   { 17, 30000,  -7225900L,   9660, -1090, -9600} ,           /* Tethys */
   { 19, 34520,  -1016300L,   3100,  -350, -3080} ,           /* Rhea */
   { 20,  2980,    -52100L,   2660,  -300, -2640} ,           /* Titan */
   { 28, 16951, -15916280L,   1790, -1080, -1420} ,           /* Phobos */
   { 29,  5347,   -662965L,   2980, -1780, -2580} ,           /* Deimos */
   { 30, 17785,     52316L, -32350, 22550, 22250} ,           /* Triton */
   { 30, 35570,    104632L,  -6280,  2100,  6730} ,           /* Triton */
   { 30, 17355,    156948L,  -2080,   550,  2050} ,           /* Triton */
   { 30, 35140,    209264L,   -740,   160,   740} };          /* Triton */


#define LUNAR_INCL 1.54242
#define LUNAR_POLE_DEC ((90 - LUNAR_INCL) * PI / 180.)
#define N_POLES 33

static int calc_planet_pole_posn( int planet_no, const double t_centuries,
                                                  double *radec)
{
                     /*     Sun      Mercury  Venus  Earth  Mars     Jupiter
                                    Saturn  Uranus  Neptune   Pluto */
   static const long ra0[N_POLES]   = {286126L, 281010L,  92760L, 0L, 317681L,
        268050L, 40589L, 257311L, 299360L, 313020L,
        269995L,                                   /* Lunar */
        268050L, 268080L, 268200L, 268720L,        /* Io Eu Ga Ca */
        40660L, 40660L, 40660L, 40660L, 40380L, 36410L, 0L, 318160L,
        257430L, 257430L, 257430L, 257430L, 257430L,   /* Ar Um Ti Ob Mi */
        317680L, 316650L,                              /* Phobos Deimos */
        299360L, 313020L,                              /* Triton Charon */
        355000L };                                     /* Phoebe */

   static const short radot[N_POLES]  = { 0, -33, 0, 0, -108, -9, -36, 0, 0, 0,
        3,                                  /* Lunar */
        -9, -9, -9, -9,                     /* Lu Io Eu Ga Ca */
        -36, -36, -36, -36, -36, -36, 0, -3949,
        0, 0, 0, 0, 0,                                 /* Ar Um Ti Ob Mi */
        -108, -108,                                    /* Phobos Deimos */
        0, 0, 0 };                                  /* Triton Charon Phoebe */

   static const long dec0[N_POLES]  = { 63870L,  61450L,  -67160L, 0L,  52886L,
        64490L, 83537L, -15175L,  43460L,   9090L,
        66539L,                                 /* Lunar */
        64500L, 64510L, 64570L, 64830L,         /* Io Eu Ga Ca */
        83520L, 83520L, 83520L, 83520L, 83550L, 83940L, 0L, 75030L,
        -15080L, -15080L, -15080L, -15080L, -15080L,  /* Ar Um Ti Ob Mi */
        52900L, 53520L,                               /* Phobos Deimos */
        41170L, 9090L, 68700L };              /* Triton Charon Phoebe */

   static const short decdot[N_POLES] = { 0, -5, 0, 0, -61, 3, -4, 0, 0, 0,
        13,                                    /* Lunar */
        3, 3, 3, 3,                            /* Io Eu Ga Ca */
        -4, -4, -4, -4, -4, -4, 0, -1143,
        0, 0, 0, 0, 0,                         /* Ar Um Ti Ob Mi */
        -61, -61,                              /* Phobos Deimos */
        0, 0, 0 };                              /* Triton Charon Phoebe */

   radec[0] = ((double)ra0[planet_no] +
                       (double)radot[planet_no] * t_centuries) * PI / 180000.;
   radec[1] = ((double)dec0[planet_no] +
                       (double)decdot[planet_no] * t_centuries) * PI / 180000.;
   return( 0);
}

/* We have data for objects 0 (Sun) through 14 (Callisto).  Two more  */
/* systems are provided for Jupiter System II and III,  plus two extra */
/* systems for Saturn. */

#define TOTAL_SYSTEMS (N_POLES + 4)

int DLL_FUNC calc_planet_orientation( int planet_no, int system_no, double jd,
                                                         double *matrix)
{
   static const long omega0[TOTAL_SYSTEMS] = { 84200L, 329680L,  19800L,
           0L, 176901L,   /* Earth, Mars */
           67100L,       /* Jupiter system I */
           227900L,      /* Saturn system I */
           203810L, 253180L, 236770L,           /* Uran Nept Plut */
           38321L,                              /* Luna */
           200390L, 35670L, 44040L, 259730L,    /* Io Euro Gany Call */
           337460L, 2820L, 10450L, 357000L,     /* Mima Ence Teth Dion */
           235160L, 189640L, 0L, 350200L,        /* Rhea Tita Hype Jape */
           156220L, 108050L, 77740L, 6770L, 30700L, /* Ar Um Ti Ob Mi */
           35060L, 79410L,                          /* Phobos Deimos */
           296530L, 56770L,                         /* Triton Charon */
           304700L,                                 /* Phoebe */
           43300L,       /* Jupiter system II */
           284950L,      /* Jupiter system III */
           104900L,     /* Saturn system II */
           38900L};    /* Saturn system III */
   static const double rates[TOTAL_SYSTEMS] = { 360. / 25.38, 6.1385025, 1.4813688,
           0, 350.891983,   /* Earth,  Mars */
           877.900,      /* Jupiter system I */
           844.3,            /* Saturn system I */
           -501.1600928, 536.3128492, -56.3623195, /* Sa Ur Ne Pl */
           13.17635815,                        /* earth's moon */
           203.488955432, 101.3747235, 50.3176081, 21.5710715, /* galilean */
           381.9945550, 262.7318996, 190.6979085, 131.5349316, /* Mi En Te Di */
           79.6900478, 22.5769768, 0., 4.5379572,    /* Rhea Tita Hype Japet */
           -142.8356681, -86.8688923, -41.3514316,   /* Ariel Umbri Titania */
           -26.7394932, -254.6906892,                /* Oberon Miranda */
           1128.844585, 285.1618970,                /* Phobos Deimos */
           -61.2572637, -56.3623195,                /* Triton Charon */
           930.8338720,                             /* Phoebe */
           870.270,      /* Jupiter system II */
           870.536,      /* Jupiter system III */
           812.0,            /* Saturn system II */
           810.7939024 };    /* Saturn system III */

   double ang, ra, dec;
   double t_cen = (jd - J2000) / 36525.;
   int i;
   double ra_dec[2];
   int idx = planet_no;
   static int prev_planet_no = -1, prev_system_no = -1;
   static double prev_jd = -1., prev_matrix[9];

   if( planet_no == 40) /* Phoebe... have to remap here to skip 8 jovians: */
      planet_no = idx = 32;

   set_identity_matrix( matrix);
   if( planet_no >= N_POLES || planet_no == 21)
      return( -1);  /* for Hyperion,  others,  return an identity matrix */

   if( planet_no == prev_planet_no && system_no == prev_system_no
                           && jd == prev_jd)
      {
      memcpy( matrix, prev_matrix, 9 * sizeof( double));
      return( 0);
      }

   prev_planet_no = planet_no;
   prev_system_no = system_no;
   prev_jd = jd;

   if( planet_no == 3)        /* handle earth with "normal" precession: */
      {
      setup_precession( matrix, 2000., 2000. + t_cen * 100.);
      for( i = 3; i < 6; i++)
         matrix[i] = -matrix[i];
      spin_matrix( matrix, matrix + 3, green_sidereal_time( jd));
      memcpy( prev_matrix, matrix, 9 * sizeof( double));
      return( 0);
      }

         /* For everybody else,  we use TD.  Only the earth uses UT. */
         /* (This correction added 5 Nov 98,  after G Seronik pointed */
         /* out an error in the Saturn central meridian.)             */

   jd += td_minus_ut( jd) / seconds_per_day;   /* correct from UT to TD */

   if( planet_no == 5 && system_no > 1)   /* Jupiter, system II & III */
      idx = TOTAL_SYSTEMS - 5 + (system_no - 1);
   if( planet_no == 6 && system_no > 1)   /* Saturn, system II & III */
      idx = TOTAL_SYSTEMS - 3 + (system_no - 1);
   ang = (double)omega0[idx] / 1000.;
   jd -= J2000;
   ang += jd * rates[idx] + 180.;
   ang *= PI / 180.;
   calc_planet_pole_posn( planet_no, t_cen, ra_dec);
   ra = ra_dec[0];
   dec = ra_dec[1];
   for( i = 0; i < N_TRIG_TERMS; i++)
      if( tterms[i].planet_no == planet_no)
         {
         double iang = ((double)tterms[i].arg / 100. +
                        (double)tterms[i].omega * t_cen / 1000.) * PI / 180.;
         double sin_iang = sin( iang);

         ra += (PI / 180.) * sin_iang * (double)tterms[i].ra_term * .001;
         dec += (PI / 180.) * cos( iang) * (double)tterms[i].dec_term * .001;
         ang += (PI / 180.) * sin_iang * (double)tterms[i].w_term * .001;
         }
   if( idx == 10)          /* Lunar */
      {
      const double edata[15] =  {
                              125.045,  -0.0529921,
                              250.089,  -0.1059842,
                              260.008,  13.0120009,
                              176.625,  13.3407154,
                              357.529,    .9856003,
                              311.589,  26.4057084,
                              134.963,  13.0649930, 0 };

      ang += -1.4e-12 * (PI / 180.) * jd * jd;
      for( i = 0; edata[i + i]; i++)
         {
         const double e_ang =
                       (PI / 180.) * (edata[i + i] + edata[i + i + 1] * jd);
         const double sin_e_ang = sin( e_ang);
         const double ra_factor[7] = { -3.8787, -.1204, .0700, -.0172,
                        0., .0072, 0. };
         const double dec_factor[7] = { 1.5419, .0239, -.0278, .0068,
                        0., -.0029, .0009 };
         const double ang_factor[7] = { 3.5610, .1208, -.0642, .0158,
                        .0252, -.0066, -.0047 };

         ra += (PI / 180.) * ra_factor[i] * sin_e_ang;
         ang += (PI / 180.) * ang_factor[i] * sin_e_ang;
         dec += (PI / 180.) * dec_factor[i] * cos( e_ang);
         }
      }

   polar3_to_cartesian( matrix, ra - PI / 2., 0.);
   polar3_to_cartesian( matrix + 3, ra - PI, HALF_PI - dec);
   polar3_to_cartesian( matrix + 6, ra, dec);

   spin_matrix( matrix, matrix + 3, ang);
   if( rates[idx] < 0.)       /* retrograde case */
      for( i = 3; i < 6; i++)
         matrix[i] *= -1.;

   memcpy( prev_matrix, matrix, 9 * sizeof( double));
   return( 0);
}
