/* elp82dat.cpp: computes lunar coords using ELP-82

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
#include "lunar.h"

/* Code to compute lunar positions using the ELP 2000-82B analytical
theory.  Available from VizieR :

https://cdsarc.cds.unistra.fr/viz-bin/cat/VI/79       */



/* 12 Dec 1998:  Following an e-mail from Luc Desamore,  I 'offset' the */
/* time used in the ELP function by: */

/*  -0.000091 (n + 26)(year-1955)^2 seconds, n being here -23.8946.  */

/* This shift (in seconds) corresponds to a difference in Delta-T between */
/* that used for VSOP,  etc.,  and that used for ELP-82.                  */

static double elp_time_offset( const double t_centuries)
{
   const double n = -23.8946;
   const double x = t_centuries * 100. + 45.;
   const double seconds_per_century = 86400. * 36525.;

   return( (-.000091 / seconds_per_century) * (n + 26.) * x * x);
}


#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

#define DMS_TO_RAD( D, M, S)   ((D + M / 60. + S / 3600.) * (PI / 180.))
#define SEC_TO_RAD( X) (X * (PI / 180.) / 3600.)

#define W1_0       DMS_TO_RAD( 218., 18., 59.95571)
#define W1_1       SEC_TO_RAD( 1732559343.73604)
#define W2_0       DMS_TO_RAD(  83.,  21., 11.67475)
#define W2_1       SEC_TO_RAD( 14643420.2632)
#define W3_0       DMS_TO_RAD( 125.,   2., 40.39816)
#define W3_1       SEC_TO_RAD( -6967919.3622)
#define T_0        DMS_TO_RAD( 100., 27., 59.22059)
#define T_1        SEC_TO_RAD( 129597742.2758)
#define OHP_0      DMS_TO_RAD( 102., 56., 14.42753)
#define OHP_1      SEC_TO_RAD( 1161.2283)

#define D_0      (W1_0 - T_0 + PI)
#define D_1      (W1_1 - T_1)
#define LPRIME_0 (T_0 - OHP_0)
#define LPRIME_1 (T_1 - OHP_1)
#define L_0      (W1_0 - W2_0)
#define L_1      (W1_1 - W2_1)
#define F_0      (W1_0 - W3_0)
#define F_1      (W1_1 - W3_1)
#define P          SEC_TO_RAD( 5029.0966)

#define W1         fund[0]
#define W2         fund[1]
#define W3         fund[2]
#define T          fund[3]
#define OHP        fund[4]
#define D          fund[5]
#define L_PRIME    fund[6]
#define L          fund[7]
#define F          fund[8]
#define ZETA       fund[17]
#define D_LINEAR   fund[18]
#define LP_LINEAR  fund[19]
#define L_LINEAR   fund[20]
#define F_LINEAR   fund[21]
#define T_LINEAR   fund[22]
#define N_FUND_COEFFS     23

#define A0         384747980.674

static void compute_lunar_polynomials( const double t_cen, double *fund,
                    const double *icoeffs)
{
   int i;
   const double *tptr = icoeffs;

   for( i = 0; i < 5; i++, tptr += 5)
      fund[i] = tptr[0] + t_cen * (tptr[1] + t_cen * (tptr[2] + t_cen *
                                 (tptr[3] + t_cen *  tptr[4])));
   for( i = 0; i < 7; i++, tptr++)
      fund[i + 9] = tptr[0] + t_cen * tptr[7];
                     /* from page 8:  compute Delaunay arguments */
   D = W1 - T + PI;
   L_PRIME = T - OHP;
   L = W1 - W2;
   F = W1 - W3;
                     /* Zeta is given by (2), p. 7: */
   ZETA = W1_0 + (W1_1 + P) * t_cen;
                 /* For series ELP4 to ELP36,  we need linear versions: */
   D_LINEAR = D_0 + t_cen * D_1;
   LP_LINEAR = LPRIME_0 + t_cen * LPRIME_1;
   L_LINEAR = L_0 + t_cen * L_1;
   F_LINEAR = F_0 + t_cen * F_1;
   T_LINEAR = T_0 + t_cen * T_1;
   fund[16] = 0.;                      /* wound up unused */
   for( i = 0; i < N_FUND_COEFFS; i++)
      {
      fund[i] = fmod( fund[i], 2. * PI);
      if( fund[i] < 0.)
         fund[i] += PI + PI;
      }
}

#define CHUNK_SIZE 10000

static double add_in_series( FILE *ifile, const int series_no,
                 const double *fund, const double prec, long n_terms)
{
   double rval = 0.;
   char *tptr, *ibuff;
   int i, coeffs[N_FUND_COEFFS];
   const int series_type = series_no / 3;
   size_t term_size, chunk_size;
   long lprec = (long)( prec * 100000.);   /* work in .00001-arcsec units */

   for( i = 0; i < N_FUND_COEFFS; i++)
      coeffs[i] = 0;
   switch( series_type)
      {
      case 0:        /* main problem */
         term_size = 8;
         break;
      case 1:  case 2:        /* Earth figure perturbations */
         term_size = 13;
         break;
      case 7: case 8: case 9: case 10: case 11:
         term_size = 12;
         break;
      case 3:  case 4:        /* Planetary perturbations */
      case 5:  case 6:        /* Planetary perturbations */
         term_size = 19;
         break;
      default:
#ifdef TEST_CODE
         printf( "??? series type %d\n", series_type);
#endif
         return( 0.);
      }
   chunk_size = CHUNK_SIZE / term_size;
   if( chunk_size > (size_t)n_terms)
      chunk_size = (size_t)n_terms;
   ibuff = tptr = (char *)calloc( (size_t)chunk_size, term_size);
   if( !ibuff)
#ifdef TEST_CODE
      {
      printf( "Failed to alloc %u * %u bytes\n", chunk_size, term_size);
      printf( "n_terms: %ld\n", n_terms);
      }
#else
      return( 0.);
#endif
   while( n_terms--)
      {
      long amplitude;
      double angle;

      if( ibuff == tptr)
         if( !fread( ibuff, chunk_size, term_size, ifile))
            {
            free( ibuff);
            return( 0.);
            }
      amplitude = *(int32_t *)tptr;
      if( amplitude < lprec && amplitude > -lprec)
         return( rval * 1.e-5);
      switch( series_type)
         {
         case 0:        /* main problem */
            for( i = 4; i < 8; i++)
               coeffs[i + 1] = (int)tptr[i];
            break;
         case 1:  case 2:        /* Earth figure perturbations */
         case 7: case 8: case 9: case 10: case 11: /* a hodgepodge of things */
            for( i = 8; i < 12; i++)
               coeffs[i + 10] = (int)tptr[i];
            if( series_type < 3)          /* yes,  there is a zeta term */
               coeffs[17] = tptr[12];
            break;
         case 3:  case 4:        /* Planetary perturbations */
            coeffs[9] = tptr[8];
            coeffs[10] = tptr[9];
            coeffs[22] = tptr[10];
            coeffs[11] = tptr[11];
            coeffs[12] = tptr[12];
            coeffs[13] = tptr[13];
            coeffs[14] = tptr[14];
            coeffs[15] = tptr[15];
            coeffs[18] = tptr[16];
            coeffs[20] = tptr[17];
            coeffs[21] = tptr[18];
            break;
         case 5:  case 6:        /* Planetary perturbations */
            coeffs[9] = tptr[8];
            coeffs[10] = tptr[9];
            coeffs[22] = tptr[10];
            coeffs[11] = tptr[11];
            coeffs[12] = tptr[12];
            coeffs[13] = tptr[13];
            coeffs[14] = tptr[14];
            coeffs[18] = tptr[15];
            coeffs[19] = tptr[16];
            coeffs[20] = tptr[17];
            coeffs[21] = tptr[18];
            break;
         }

      if( series_type)
         angle = (double)*(int32_t *)( tptr + 4) * (PI / 180.) / 100000.;
      else
         angle = 0.;
      for( i = 0; i < N_FUND_COEFFS; i++)
         if( coeffs[i])
            {
            if( coeffs[i] == 1)
               angle += fund[i];
            else if( coeffs[i] == -1)
               angle -= fund[i];
            else
               angle += (double)coeffs[i] * fund[i];
            coeffs[i] = 0;
            }
      if( series_no == 2)     /* main distance theory is oddball */
         rval += (double)amplitude * cos( angle);
      else
         rval += (double)amplitude * sin( angle);
      tptr += term_size;
      if( (size_t)( tptr - ibuff) == chunk_size * term_size)     /* end of line */
         tptr = ibuff;
      }
   free( ibuff);
   return( rval * 1.e-5);
}

#define ELP_DATA_HEADER struct elp_data_header

ELP_DATA_HEADER
   {
   int32_t offsets[37 * 2];
   double poly_coeffs[5 * 5 + 7 * 2];
   };

static int get_elp_values( FILE *ifile, const double t_cen,
                     const double prec0, double *ovals)
{
   int i;
   double addition;
   double fund[N_FUND_COEFFS];
   ELP_DATA_HEADER *hdr = (ELP_DATA_HEADER *)malloc( sizeof( ELP_DATA_HEADER));

   if( !hdr)
      return( -1);
   if( !fread( hdr, sizeof( ELP_DATA_HEADER), 1, ifile))
      {
      free( hdr);
      return( -2);
      }
   compute_lunar_polynomials( t_cen, fund, hdr->poly_coeffs);

               /* First longitude term has to be 'adjusted': */
   ovals[0] = fund[0] + (22639.58578 * PI / 180.) * sin( fund[7]) / 3600.;
   ovals[1] = 0.;
   ovals[2] = 385000.52719;
   for( i = 0; i < 36; i++)
      {
      double prec = prec0;
      int series_type = i / 3;

      if( prec != 0.)
         {
         if( i % 3 == 2)         /* distance term: cvt to kilometers */
            prec *= A0 / 1000.;
         else                       /* angular term: cvt to arcseconds */
            prec *= (180. * 3600. / PI);
         }
      fseek( ifile, hdr->offsets[i + i], SEEK_SET);
      if( hdr->offsets[i + i + 1])
         addition = add_in_series( ifile, i, fund, prec,
                                               hdr->offsets[i + i + 1]);
      else
         addition = 0.;
      if( series_type == 2 || series_type == 4 ||
          series_type == 6 || series_type == 8)
         addition *= t_cen;
      if( series_type == 11)
         addition *= t_cen * t_cen;
      if( (i % 3) == 2)
         ovals[2] += addition;
      else
         ovals[i % 3] += addition * (PI / 180.) / 3600.;
      }
   free( hdr);
   return( 0);
}

               /* Laskar's coeffs for precession,  p. 12: */
#define P_0          0.10180391e-4
#define P_1          0.47020439e-6
#define P_2         -0.5417367e-9
#define P_3         -0.2507948e-11
#define P_4          0.463486e-14
#define Q_0         -0.113469002e-3
#define Q_1          0.12372674e-6
#define Q_2          0.1265417e-8
#define Q_3         -0.1371808e-11
#define Q_4         -0.320334e-14

int DLL_FUNC compute_elp_xyz( FILE *ifile, const double t_cen,
                   const double prec, double *ecliptic_xyz_2000)
{
   double uvr[3];
   int i, close_file = 0, rval;
   const double adjusted_t_cen = t_cen + elp_time_offset( t_cen);

   if( !ifile)
      {
      ifile = fopen( "elp82.dat", "rb");
      if( !ifile)
         return( -1);
      close_file = 1;
      }
   fseek( ifile, 0L, SEEK_SET);
   rval = get_elp_values( ifile, adjusted_t_cen, prec, uvr);
   if( rval)
      for( i = 0; i < 4; i++)
         *ecliptic_xyz_2000++ = 0.;
      else
         {
         const double x = uvr[2] * cos( uvr[0]) * cos( uvr[1]);
         const double y = uvr[2] * sin( uvr[0]) * cos( uvr[1]);
         const double z = uvr[2] *           sin( uvr[1]);
         double p = 0., q = 0., twice_root_pq_term;
         static const double p_coeff[5] = { P_4, P_3, P_2, P_1, P_0 };
         static const double q_coeff[5] = { Q_4, Q_3, Q_2, Q_1, Q_0 };
         double matrix[9];

         for( i = 0; i < 5; i++)
            {
            p = p * adjusted_t_cen + p_coeff[i];
            q = q * adjusted_t_cen + q_coeff[i];
            }
         p *= adjusted_t_cen;
         q *= adjusted_t_cen;
         twice_root_pq_term = 2. * sqrt( 1. - p * p - q * q);
         matrix[0] = 1. - 2. * p * p;
         matrix[1] = matrix[3] = 2. * p * q;
         matrix[2] = p * twice_root_pq_term;
         matrix[6] = -matrix[2];
         matrix[7] = q * twice_root_pq_term;
         matrix[5] = -matrix[7];
         matrix[4] = 1. - 2. * q * q;
         matrix[8] = matrix[0] - 2. * q * q;
         for( i = 0; i < 9; i += 3)
            *ecliptic_xyz_2000++ =
                           matrix[i] * x + matrix[i + 1] * y + matrix[i + 2] * z;
         *ecliptic_xyz_2000++ = uvr[2];         /* give the radius,  too */
         }
   if( close_file)
      fclose( ifile);
   return( rval);
}

#ifdef TEST_CODE

#define J2000 2451545.0
#define TEST_PREC (5.e-5 * (PI / 180.) / 3600.)
#define TEST_T       -.13400608189

void main( int argc, char **argv)
{
   FILE *ifile = fopen( "elp82.dat", "rb");
   double xyz[4], t;
   int i;
   const double prec = (argc == 1 ? 0 : atof( argv[1]));
   const double results_from_book[5][3] = {
             -361602.98536,    44996.99510,    -30696.65316,
             -363132.34248,    35863.65378,    -33196.00409,
             -371577.58161,    75271.14315,    -32227.94618,
             -373896.15893,   127406.79129,    -30037.79225,
             -346331.77361,   206365.40364,    -28502.11732 };
   if( !ifile)
      {
      printf( "elp82.dat not found\n");
      exit( -1);
      }

   for( i = 0; i < 5; i++)
      {
      t = 2469000.5 - 20000. * (double)i;
      compute_elp_xyz( ifile, (t - J2000) / 36525., prec, xyz);
      printf( "%.3lf: %15.5lf %15.5lf %15.5lf  %15.5lf\n", t,
                           xyz[0], xyz[1], xyz[2], xyz[3]);
      printf( "             %15.5lf %15.5lf %15.5lf\n",
                           results_from_book[i][0] - xyz[0],
                           results_from_book[i][1] - xyz[1],
                           results_from_book[i][2] - xyz[2] );
      }
   for( i = 2; i < argc; i++)
      {
      t = atof( argv[i]);
      compute_elp_xyz( ifile, (t - J2000) / 36525., 0., xyz);
      printf( "%.3lf: %15.5lf %15.5lf %15.5lf  %15.5lf\n", t,
                           xyz[0], xyz[1], xyz[2], xyz[3]);
      }
   fclose( ifile);
}
#endif
