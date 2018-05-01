/* de_plan.cpp: functions for planetary ephems from PS-1996

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "watdefs.h"
#include "lunar.h"
#include "get_bin.h"

#define POISSON struct poisson
#define POISSON_HEADER struct poisson_header

POISSON
   {
   double tzero, dt;
   int16_t n_blocks, total_fqs;    /* mx, imax are both = 2 in all cases... */
   double secular[12];
   double *fqs, *terms;
   int16_t nf[3];
   };

#pragma pack(2)

POISSON_HEADER
   {
   double tzero, dt;
   int16_t nf[3];
   int16_t n_blocks, total_fqs;
   };

#pragma pack()

/*
 In the PS_1996.DAT file,  the coefficients are stored in sets of six
doubles.  Usually,  they're small values that can be stored in a
single byte;  sometimes,  a int16_t or int32_t is needed.  Of course,
there are a few that really do require a full,  eight-byte double.
To handle this,  a set of six doubles is stored starting with a
int16_t integer;  it is parsed out as six pairs of bits.  Each pair
of bits tells you if a given coefficient is stored as a char,  int16_t,
int,  or full double.  Finally,  the number of bytes that were parsed
is returned.

   (2012 Jan 17) The function is used exactly once,  from within this
file,  and therefore should be 'static inline'.

   (2018 May 1) Puzzled out a longstanding MSVC optimization issue.
With optimizations turned on,  you _cannot_ simply extract a signed
32-bit integer;  it'll always be treated as unsigned.  Various
workarounds are, of course,  possible.  This was the simplest.  Not all
that efficient, but elderly MSVC isn't a high priority these days. */

static inline int unpack_six_doubles( double *ovals, const char *ibuff)
{
   int i;
   uint16_t flags = get16bits( ibuff);
   const char *iptr = ibuff + 2;

   for( i = 0; i < 6; i++, flags >>= 2)
      switch( flags & 3)
         {
         case 0:
            *ovals++ = (double)*iptr++;
            break;
         case 1:
            *ovals++ = (double)get16sbits( iptr);
            iptr += 2;
            break;
         case 2:
            *ovals++ = (double)get32sbits( iptr);
            iptr += 4;
#if defined(_MSC_VER) && _MSC_VER < 1900
            {
            const double two_31 = 2147483648.;

            if( ovals[-1] > two_31)
               ovals[-1] -= 2. * two_31;
            }
#endif
            break;
         case 3:
            *ovals++ = get_double( iptr);
            iptr += 8;
            break;
         }
   return( (int)( iptr - ibuff));
}

/* This function loads up the series for a given planet,  at a given time,
from a given file.  If you've already got PS_1996.DAT opened,  as Guide
usually does,  then you can just pass in the file pointer.  If not,  no
problem;  the file will be opened (and closed) for you.

   Next,  a buffer of ten int32_t integers is read in;  this provides the
offset within the file for the data concerning each planet.  Some checking
is done to figure out which 'block' covers the current jd;  of course,
it's possible that none will,  in which case we close up shop and return
a NULL.

   Otherwise,  we figure out how much memory is needed and allocate one
massive buffer to hold it all.  This may look a little odd,  but it
serves several purposes.  First,  you only have to check one memory
allocation for success or failure,  which simplifies the code a little.
Second,  freeing up the buffer when you're done with it is much easier
(see the following unload_ps1996_series() function).  I use this
"allocate one buffer" trick a lot.

   Next,  we read in the data.  The coefficients are read using the
above unpack_six_doubles() function.
*/

void * DLL_FUNC load_ps1996_series( FILE *ifile, double jd, int planet_no)
{
   int32_t offset;
   long jump = 0;
   int read_error = 0;
   int16_t block_sizes[30];
   int block, i, close_file = 0;
   POISSON_HEADER header;
   char *tbuff, *tptr;
   POISSON p;
   POISSON *rval;

   if( !ifile)
      {
      ifile = fopen( "ps_1996.dat", "rb");
      if( !ifile)
         return( NULL);
      close_file = 1;
      }

   fseek( ifile, (long)(planet_no - 1) * (long)sizeof( int32_t), SEEK_SET);
   if( !fread( &offset, 1, sizeof( int32_t), ifile))
      read_error = -1;
   fseek( ifile, offset, SEEK_SET);
   if( !fread( &header, sizeof( POISSON_HEADER), 1, ifile))
      read_error = -1;
   block = (int)floor( (jd - header.tzero) / header.dt);
   if( block < 0 || block >= header.n_blocks || read_error)
      {
      if( close_file)
         fclose( ifile);
      return( NULL);            /* outta bounds */
      }

   p.tzero = header.tzero + (double)block * header.dt;
   p.dt = header.dt;
   p.n_blocks = header.n_blocks;
   p.total_fqs = header.total_fqs;
   for( i = 0; i < 3; i++)
      p.nf[i] = header.nf[i];

   rval = (POISSON *)malloc( sizeof( POISSON) +
                            (size_t)header.total_fqs * 7 * sizeof( double));
   if( !rval)
      {
      if( close_file)
         fclose( ifile);
      return( NULL);
      }
   p.fqs = (double *)( rval + 1);
   p.terms = p.fqs + header.total_fqs;
   if( !fread( p.fqs, (size_t)header.total_fqs, sizeof( double), ifile))
      read_error = -3;
   if( !fread( block_sizes, (size_t)header.n_blocks, sizeof( int16_t), ifile))
      read_error = -4;
   for( i = 0; i < block; i++)
      jump += (long) block_sizes[i];
   fseek( ifile, jump, SEEK_CUR);
   tbuff = (char *)malloc( (size_t)block_sizes[block]);
   if( tbuff && !fread( tbuff, (size_t)block_sizes[block], 1, ifile))
      read_error = -5;
   if( !tbuff || read_error)
      {
      free( rval);
      if( close_file)
         fclose( ifile);
      if( tbuff)
         free( tbuff);
      return( NULL);
      }
   memcpy( p.secular, tbuff, 12 * sizeof( double));
   tptr = tbuff + 12 * sizeof( double);
   for( i = 0; i < p.total_fqs; i++)
      tptr += unpack_six_doubles( p.terms + i * 6, tptr);
   free( tbuff);
   if( close_file)
      fclose( ifile);
   memcpy( rval, &p, sizeof( POISSON));
   return( rval);
}

int DLL_FUNC unload_ps1996_series( void *p)
{
   free( p);
   return( 0);
}

/* I will have to mumble concerning much of what the following function
does.  There are parts I don't understand,  and parts that are
mathematically tricky enough that I don't want to spend much time
explaining them.  That doesn't leave much.  But here goes:

   The function first makes sure that the Poisson-series data pointer
you passed in,  iptr,  actually covers the jd you passed in;  if it
doesn't,  you get a return value of -1.  If all is well,  it goes
through and computes the Cartesian position (x, y, z) of the object,
in AU,  on the date in question,  in heliocentric J2000 coordinates.
x = state_vect[0], y = state_vect[1], z = state_vect[2].  If you pass
in a non-zero value for compute_velocity,  then it will also figure
out the velocity,  in AU/day,  as vx = state_vect[3],  vy = state_vect[4],
vz = state_vect[5].
*/

int DLL_FUNC get_ps1996_position( const double jd, const void *iptr,
                        double *state_vect, const int compute_velocity)
{
   const POISSON *p = (const POISSON *)iptr;
   const double x = 2. * (jd - p->tzero) / p->dt - 1.;
   const double fx = x * p->dt / 2.;
   const double *fq_ptr = p->fqs;
   const double *sec_ptr = p->secular;
   const double *term_ptr = p->terms;
   int i, j, m;
   double wx = 1.;
   double xpower[5];

   if( jd < p->tzero || jd > p->tzero + p->dt)
      return( -1);
   xpower[0] = xpower[1] = 1.;
   for( i = 2; i < 5; i++)
      xpower[i] = xpower[i - 1] * x;
   for( i = 0; i < 3; i++)       /* secular terms first: */
      {
      state_vect[i] = 0.;
      wx = 1.;
      for( j = 0; j < 4; j++)
         {
         if( compute_velocity)
            {
            if( j)
               state_vect[i + 3] += (double)j * (*sec_ptr) * xpower[j];
            else
               state_vect[i + 3] = 0.;
            }
         state_vect[i] += wx * (*sec_ptr++);
         wx *= x;
         }
      if( compute_velocity)
         state_vect[i + 3] *= 2. / p->dt;
      }

   wx = 1.;
   for( m = 0; m < 3; m++)
      {
      double new_sums[6];
      const double velocity_scale = (double)( m * 2) * xpower[m] / p->dt;

      for( i = 0; i < 6; i++)
         new_sums[i] = 0.;
      for( j = 0; j < p->nf[m]; j++)
         {
         const double amplitude = *fq_ptr++;
         const double f = amplitude * fx;
         const double cos_term = cos( f);
         const double sin_term = sin( f);

         for( i = 0; i < 3; i++, term_ptr += 2)
            {
            new_sums[i] += term_ptr[0] * cos_term + term_ptr[1] * sin_term;
            if( compute_velocity)
               new_sums[i + 3] += amplitude *
                          (term_ptr[1] * cos_term - term_ptr[0] * sin_term);
            }
         }

      for( i = 0; i < 3; i++)
         {
         state_vect[i] += new_sums[i] * wx;
         if( compute_velocity)
            {
            state_vect[i + 3] += new_sums[i + 3] * wx;
            if( m)
               state_vect[i + 3] +=
                               velocity_scale * new_sums[i];
            }
         }
      wx *= x;
      }

   for( i = 0; i < (compute_velocity ? 6 : 3); i++)
      state_vect[i] *= 1.e-10;         /* cvt to AU,  and to AU/day */
   return( 0);
}

#ifdef TEST_CODE

/* Example run,  giving heliocentric J2000 equatorial state
vector for Mars as of JD 2458239.5 = 2018 May 1 TDB :

../lunar/de_plan 4 2458239.5
Series loaded
-0.435280249 -1.303106084 -0.585949939  1.493616984
0.013913669 -0.002480622 -0.001513357
*/

int main( const int argc, const char **argv)
{
   double state_vect[6], r;
   const double t0 = atof( argv[2]);
   void *p;
   FILE *ifile;
   int rval;

   if( argc < 3)
      {
      printf( "You must supply a planet index (1=Mercury...8=Neptune)\n");
      printf( "and a JD as command-line arguments.\n");
      return( -3);
      }
   ifile = fopen( "ps_1996.dat", "rb");
   if( !ifile)
      {
      printf( "ps_1996.dat not loaded\n");
      return( -2);
      }
   p = load_ps1996_series( ifile, t0, atoi( argv[1]));
   fclose( ifile);
   if( !p)
      {
      printf( "Outside bounds of series\n");
      return( -1);
      }
   printf( "Series loaded\n");
   rval = get_ps1996_position( t0, p, state_vect, 1);
   unload_ps1996_series( p);
   if( rval)
      printf( "Error occurred: rval %d\n", rval);
   else
      {
      r = state_vect[0] * state_vect[0] + state_vect[1] * state_vect[1] +
                        state_vect[2] * state_vect[2];

      printf( "%.9lf %.9lf %.9lf  %.9lf\n", state_vect[0], state_vect[1],
            state_vect[2], sqrt( r));
      printf( "%.9lf %.9lf %.9lf\n", state_vect[3], state_vect[4],
            state_vect[5]);
      }
   return( rval);
}
#endif
