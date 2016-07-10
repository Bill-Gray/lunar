/* colors.cpp: functions for color conversions

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

#include "colors.h"

/*
Code to do assorted rough transformations between color systems,
such as getting a B-V value given a V-I color.  I started out with
FORTRAN code supplied by Brian Skiff,  converted it to C,  and added
some inverse transformations not previously available.  Note that these
are _extremely_ rocky,  because they happen when one color doesn't
vary much with respect to another color.  They should not be taken
too seriously.  Brian comments:

   "High-order polynomials for transforming amongst colors are given
by Caldwell et al.:

1993SAAOC..15....1C
  CALDWELL J.A.R., COUSINS A.W.J., AHLERS C.C., VAN WAMELEN P., MARITZ E.J.
  South African Astron. Obs. Circ., 15, 1-29 (1993)
  Statistical relations between photometric colours of common types of stars in
    the UBV (RI)c, JHK and uvby systems."

   Also note that the V-I to V-R function is an exact inverse of the V-R to
V-I one,  and similarly for the B-V to V-I and V-I to B-V functions.  That
mathematical exactness is a result of using an inverse root-searching
algorithm,  though,  _not_ because the color conversion is so good.  Caldwell
gives polynomials for B-V to V-R and vice versa.  Feed a value through one
function and its "inverse" and you won't get exactly your original value.
The difference is usually about .005 mag or so,  but can be more.

   You'll see three separate test main( ) functions at the bottom,
originally used in testing these functions and now possibly helpful in
showing how they are used.  If you define LONEOS_PHOT (the last test
main( )),  you'll get a program to read loneos.phot and add colors
where no colors were previously available,  flagged so you can know
where they came from.
*/

#define COLOR_ORDER 13

static double compute_color_polynomial( const double ival,
              const double *coeffs, const double low_limit,
              const double high_limit)
{
   double rval = 0., power = 1.;
   int order = COLOR_ORDER;

   if( ival < low_limit || ival > high_limit)
      rval = 99.;
   else while( order--)
      {
      rval += (*coeffs++) * power;
      power *= ival - 1.;
      }
   return( rval + 1.);
}

static double compute_inverse_color_polynomial( const double ival,
              const double *coeffs, double ilow_limit, double ihigh_limit)
{
   double color_low, color_high, rval, color_rval;
   double high_limit = ihigh_limit;
   double low_limit = ilow_limit;
   int max_iterations = 100;

   color_low = compute_color_polynomial( ilow_limit, coeffs,
                        low_limit, high_limit);
   color_high = compute_color_polynomial( ihigh_limit, coeffs,
                        low_limit, high_limit);
   if( ival < color_low || ival > color_high)
      rval = 99.;
   else do
      {
      rval = (ival - color_low) / ( color_high - color_low);
      rval += rval * (rval -.5) * (rval - 1.);
      rval = low_limit + (high_limit - low_limit) * rval;
      color_rval = compute_color_polynomial( rval, coeffs,
                     ilow_limit, ihigh_limit);
      if( color_rval < ival)
         {
         low_limit = rval;
         color_low = color_rval;
         }
      else
         {
         high_limit = rval;
         color_high = color_rval;
         }
      max_iterations--;
      }
      while( high_limit - low_limit > 1.e-6 && max_iterations);
//    while( color_high - color_low > 1.e-6 && max_iterations);

   if( !max_iterations)       /* didn't converge on an answer */
      rval = 99.;
   return( rval);
}

static const double coeffs_vi_to_bv[COLOR_ORDER] = {
    -0.6865072E-01,  0.8837997E+00,   -0.3889774E+00,
    -0.4998126E-02,  0.3867544E+00,   -0.5422331E+00,
    -0.8926476E-01,  0.5194797E+00,   -0.2044681E+00,
    -0.1009025E+00,  0.9543256E-01,   -0.2567529E-01,
     0.2393742E-02 };

double v_minus_i_to_b_minus_v( const double v_minus_i)
{
   return( compute_color_polynomial(
                    v_minus_i, coeffs_vi_to_bv, -.23, 3.70));
}

double b_minus_v_to_v_minus_i( const double b_minus_v)
{
   return( compute_inverse_color_polynomial(
                    b_minus_v, coeffs_vi_to_bv, -.23, 3.70));
}

static const double coeffs_vi_to_vr[COLOR_ORDER] = {
   -0.4708373E+00,    0.5920728E+00,   -0.1095294E-01,
   -0.2281118E+00,   -0.9372892E-01,    0.1931393E+00,
    0.5077253E-01,   -0.9927284E-01,    0.8560631E-02,
    0.1922702E-01,   -0.7201880E-02,    0.7743020E-03, 0. };

double v_minus_i_to_v_minus_r( const double v_minus_i)
{
   return( compute_color_polynomial(
                    v_minus_i, coeffs_vi_to_vr, -.30, 4.00));
}

double v_minus_r_to_v_minus_i( const double v_minus_r)
{
   return( compute_inverse_color_polynomial(
                    v_minus_r, coeffs_vi_to_vr, -.30, 4.00));
}

double v_minus_r_to_b_minus_v( const double v_minus_r)
{
   static const double coeffs[COLOR_ORDER] = {
   0.4860429E+00,   0.6904008E+00,  -0.1229411E+01,   0.2990030E+01,
   0.7104513E+01,  -0.1637799E+02,  -0.2977123E+02,   0.4390751E+02,
   0.6145810E+02,  -0.5265358E+02,  -0.6135921E+02,   0.2297835E+02,
   0.2385013E+02};

   return( compute_color_polynomial( v_minus_r, coeffs, -.10, 1.75));
}

double b_minus_v_to_v_minus_r( const double b_minus_v)
{
   static const double coeffs[COLOR_ORDER] = {
  -0.4140951E+00,   0.7357165E+00,  -0.5242979E-01,  -0.6293304E+00,
   0.2332871E+01,   0.3812365E+01,  -0.5082941E+01,  -0.6520325E+01,
   0.4817797E+01,   0.5065505E+01,  -0.1706011E+01,  -0.1568243E+01, 0. };

   return( compute_color_polynomial( b_minus_v, coeffs, -.23, 1.95));
}

   /* This function derived from data on p 57,  _Intro & Guide to the Data_ */
double johnson_b_minus_v_from_tycho_b_minus_v( const double b_v_t)
{
   double delta = 0.;

   if( b_v_t < -.2 || b_v_t > 1.8)
      return( 99.);        /* no reasonable transformation possible */
   if( b_v_t < .1)
      delta = -.006 + .006 * (b_v_t + .2) / .3;
   else if( b_v_t < .5)
      delta = .046 * (b_v_t - .1) / .4;
   else if( b_v_t < 1.4)
      delta = .046 - .054 * (b_v_t - .5) / .9;
   else if( b_v_t < 1.8)
      delta = -.008 - .024 * (b_v_t - 1.4) / .4;
   return( .85 * b_v_t + delta);
}

/*
This function derived from data on p 57,  _Intro & Guide to the Data_.
It's probably better to use the tycho_to_johnson_colors( ) function,
in COLORS2.CPP,  instead.  I wrote the following some time ago,  and it's
somewhat obsolete now.
*/

double johnson_v_from_tycho_b_minus_v( const double b_v_t, const double tycho_v)
{
   double delta = 0.;

   if( b_v_t < -.2 || b_v_t > 1.8)
      return( 99.);        /* no reasonable transformation possible */
   if( b_v_t < .1)
      delta =  .014 - .014 * (b_v_t + .2) / .3;
   else if( b_v_t < .5)
      delta = -.005 * (b_v_t - .1) / .4;
   else if( b_v_t < 1.4)
      delta = -.005;
   else if( b_v_t < 1.8)
      delta = -.005 - .010 * (b_v_t - 1.4) / .4;
   return( tycho_v - .09 * b_v_t + delta);
}


#ifdef SIMPLE_TEST_PROGRAM

#include <stdio.h>
#include <stdlib.h>

int main( const int argc, const char **argv)
{
   if( argc != 2)
      {
      printf( "'colors' takes a color index as a command-line argument,\n");
      printf( "and outputs corresponding indices for other systems.\n");
      return( -1);
      }

   const double ival = atof( argv[1]);

   printf( "v_minus_i_to_b_minus_v( %f) = %f\n",
            ival, v_minus_i_to_b_minus_v( ival));
   printf( "v_minus_i_to_v_minus_r( %f) = %f\n",
            ival, v_minus_i_to_v_minus_r( ival));
   printf( "b_minus_v_to_v_minus_i( %f) = %f\n",
            ival, b_minus_v_to_v_minus_i( ival));
   printf( "v_minus_r_to_v_minus_i( %f) = %f\n",
            ival, v_minus_r_to_v_minus_i( ival));
   printf( "v_minus_r_to_b_minus_v( %f) = %f\n",
            ival, v_minus_r_to_b_minus_v( ival));
   printf( "b_minus_v_to_v_minus_r( %f) = %f\n",
            ival, b_minus_v_to_v_minus_r( ival));
   return( 0);
}
#endif

#ifdef GRAPHING_PROGRAM

#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <graph.h>

int main( const int argc, const char **argv)
{
   int function = atoi( argv[1]), i;
   char buff[200];
   FILE *ifile = fopen( "loneos.phot", "rb");
   double limits[12] = { -.4, 3.9, -.4, 2.,        /* V-I to B-V */
                         -.5, 4.2, -.4, 2.1,       /* V-I to V-R */
                         -.3, 2.,  -.4, 2.2 };     /* V-R to B-V */
   double *lim = limits + function * 4;
   int xcolumns[3] = { 75, 75, 69 };
   int ycolumns[3] = { 63, 69, 63 };

   _setvideomode( _VRES16COLOR);
   _setcolor( 2);
   for( i = 0; i < 640; i++)
      {
      double ocolor, icolor = lim[0] + (lim[1] - lim[0]) * (double)i / 640.;

      if( function == 0)
         ocolor = v_minus_i_to_b_minus_v( icolor);
      else if( function == 1)
         ocolor = v_minus_i_to_v_minus_r( icolor);
      else
         ocolor = v_minus_r_to_b_minus_v( icolor);
      _setpixel( (short)i, (short)( 480. * (ocolor-lim[3]) / (lim[2]-lim[3])));
      }

   _setcolor( 3);
   while( fgets( buff, sizeof( buff), ifile))
      {
      if( buff[xcolumns[function]] == '.' && buff[ycolumns[function]] == '.')
         {
         double icolor = atof( buff + xcolumns[function] - 2);
         double ocolor = atof( buff + ycolumns[function] - 2);
         int xpixel = (int)( 640. * (icolor-lim[0]) / (lim[1]-lim[0]));
         int ypixel = (int)( 480. * (ocolor-lim[3]) / (lim[2]-lim[3]));
         short ix, iy;
         short deltas[9 * 2] = { 0,0, 1,0, 0,1, -1,0, 0,-1,
                        1,1, 1,-1, -1,-1, -1,1 };

         for( i = 0; i < 18; i += 2)
            {
            ix = (short)( xpixel + deltas[i]);
            iy = (short)( ypixel + deltas[i + 1]);
            if( !_getpixel( ix, iy))
               i = 99;
            }
         if( i != 18)
            _setpixel( ix, iy);
         }
      buff[xcolumns[function]] = buff[ycolumns[function]] = ' ';
      }
   fclose( ifile);

   getch( );
   _setvideomode( _DEFAULTMODE);
   return( 0);
}

#endif

#ifdef LONEOS_PHOT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VI_OFFSET    84
#define VR_OFFSET    77
#define UB_OFFSET    70
#define BV_OFFSET    63
#define LINE_END     92

int main( const int argc, const char **argv)
{
   char buff[200];
   FILE *ifile = fopen( "loneos.pho", "rb"), *ofile;
   size_t i;

   if( !ifile)
      {
      printf( "Couldn't open 'loneos.pho'\n");
      exit( -1);
      }
   ofile = fopen( "loneos2.pho", "wb");
   while( fgets( buff, sizeof( buff), ifile))
      {
      if( buff[27] == '.' && buff[58] == '.' && buff[30] == ' ' &&
                                 strlen( buff) > 60)
         {
         int bv_flag, vr_flag, vi_flag;
//       int ub_flag;
         double b_minus_v = 99., v_minus_r = 99., v_minus_i = 99.;
//       double u_minus_b = 99.;

         for( i = 0; buff[i] >= ' '; i++)
            ;
         while( i < sizeof( buff))
            buff[i++] = '\0';
         bv_flag = (buff[BV_OFFSET + 2] == '.');
//       ub_flag = (buff[UB_OFFSET + 2] == '.');
         vr_flag = (buff[VR_OFFSET + 2] == '.');
         vi_flag = (buff[VI_OFFSET + 2] == '.');
         if( vi_flag)
            {
            v_minus_i = atof( buff + VI_OFFSET);
            if( !bv_flag)
               {
               b_minus_v = v_minus_i_to_b_minus_v( v_minus_i);
               bv_flag = 'i';
               }
            if( !vr_flag)
               {
               v_minus_r = v_minus_i_to_v_minus_r( v_minus_i);
               vr_flag = 'i';
               }
            }
         else
            if( vr_flag)
               {
               v_minus_r = atof( buff + VR_OFFSET);
               v_minus_i = v_minus_r_to_v_minus_i( v_minus_r);
               vi_flag = 'r';
               if( !bv_flag)
                  {
                  b_minus_v = v_minus_r_to_b_minus_v( v_minus_r);
                  bv_flag = 'r';
                  }
               }
            else     /* no V-R or V-I:  Only B-V is available */
               {
               b_minus_v = atof( buff + BV_OFFSET);
               v_minus_r = b_minus_v_to_v_minus_r( b_minus_v);
               v_minus_i = b_minus_v_to_v_minus_i( b_minus_v);
               vi_flag = vr_flag = 'b';
               }
         if( bv_flag > 1 && b_minus_v < 98.)
            {
            sprintf( buff + BV_OFFSET, "%5.2f", b_minus_v);
            buff[BV_OFFSET + 6] = bv_flag;
            }
         if( vr_flag > 1 && v_minus_r < 98.)
            {
            sprintf( buff + VR_OFFSET, "%5.2f", v_minus_r);
            buff[VR_OFFSET + 6] = vr_flag;
            }
         if( vi_flag > 1 && v_minus_i < 98.)
            {
            sprintf( buff + VI_OFFSET, "%5.2f", v_minus_i);
            buff[VI_OFFSET + 6] = vi_flag;
            }
         for( i = 0; i < LINE_END; i++)
            if( !buff[i])
               buff[i] = ' ';
         if( !buff[LINE_END])
            {
            for( i = LINE_END; i && buff[i - 1] == ' '; i--)
               ;
            buff[i] = '\0';
            }
         fprintf( ofile, "%s\n", buff);
         }
      else
         fprintf( ofile, "%s", buff);
      memset( buff, 0, sizeof( buff));
      }
   fclose( ifile);
   fclose( ofile);
   return( 0);
}
#endif
