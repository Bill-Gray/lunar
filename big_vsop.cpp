/* big_vsop.cpp: functions for analytic (VSOP87) planetary ephems

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
#include <stdint.h>
#include "watdefs.h"
#include "lunar.h"

#define PI 3.141592653589793238462643383279502884197169399375105

/* Computes mean heliocentric ecliptic latitude/longitude of date
and radius for a planet,  using the full VSOP87 theory.

   BE WARNED that for most purposes,  rather than using VSOP (in
either its short or long forms),  it's well to use either the
PS-1996 series or DE ephemerides.  Source code for both is
available from

http://www.projectpluto.com/source.htm

   VSOP is based on the increasingly out-of-date DE-200 ephemerides
and will give you only eight planets (Mercury-Neptune).  JPL DE
will add in Pluto and the moon,  and is much faster to compute.
(VSOP,  ELP,  and PS-1996 all require one to sum up lots of trig
terms;  they express planetary/lunar positions as Poisson series.)
The only real down side to DE is file size (as much as half a GByte
if one wants to cover a six-thousand-year time span.)

   PS-1996 fits in a mere 132 KBytes,  is based on DE-4xx,  and
gives Mercury-Pluto for 1900 to 2100.  It's still a bit slow.

   VSOP-87 is not totally useless.  It can be packed into a small
file,  and you can sum up just a few terms to get an approximate
position.  (My numerical integration software uses this trick to
determine what planets are "significant perturbers" in a given
situation;  see sm_vsop.cpp in the Find_Orb source code for this.)

   In both 'big_vsop.bin' and 'vsop.bin',  for each of the eight
planets (Mercury through Neptune), data is provided in the
"usual" VSOP form,  as a Poisson series (a mix of a Fourier-like
series of trig terms and a Taylor-like power series.)  Thus,  for
example, the ecliptic latitude of a planet would be computed as

latitude = (sum of trig terms0)
         + (sum of trig terms1) * t
         + (sum of trig terms2) * t^2
         + (sum of trig terms3) * t^3
         + (sum of trig terms4) * t^4
         + (sum of trig terms5) * t^5

   ...with similar series given for ecliptic longitude and heliocentric
radius.

   't' = (jd - 2541545.0) / 365250,  the difference in millennia
between the time in question and 1.5 January 2000.  The 'trig terms'
are all of the form amplitude * cos( angle + rate * t).  Almost all of
'big_vsop.bin' consists of the values for 'amplitude',  'angle' and
'rate',  with a small header describing just where the 'amplitude',
'angle' and 'rate' data for a given series begins and ends.

   As you can see,  each value requires summing up six series,  then
multiplying by some integer power of t.  We've got six series per
value,  three values per planet,  and eight planets... therefore,
6 * 3 * 8 = 144 series.  Any given planet requires 18 series.

   The 'big_vsop' header could,  in theory,  give you the starting
and ending location of each series,  as short integer values.  You'd
then have 144 * 2,  or 288,  values in the header.  In reality,  all
I did was store the beginning of each series;  you figure out where
it ends by looking at the beginning of the next series.  That brings
us down to 145 values in the header (the last value giving you where
the final series actually ends.)

   Since each header entry is a short int,  we're looking at 290 bytes.

   In a possibly misguided effort to save space (well,  it _did_ make
lots of sense back in 1993!),  I store the header data needed for
one particular planet in RAM at a given time.  If someone asks for
a "new planet" (in big_vsop.cpp,  'if( curr_planet != planet)'),
then we dig through 'big_vsop.bin' for the required header data
and store it in a static array.  That requirement is for nineteen
short integers:  there are three values to be computed (lat/lon/r)
and six series for each of them (1, t, t^2, ...t^5),  and we need
to know where the last of them ends.  So the static array 'cache'
is dimensioned for 19 shorts.

   We then use that data stored in 'cache' to go forth and grab
VSOP data.  For longitude,  the coefficients for the (1, t, ...t^5)
series are stored at offsets indicated by cache[0], cache[1], ...
cache[5],  with the latter ending at cache[6] (that is,  the t^5
series would have cache[6]-cache[5] terms.)  For latitude,
the coefficients would be at offsets indicated by cache[6...11],
and for heliocentric radius,  cache[12...17],  with the last
t^5 term having cache[18]-cache[17] terms.

   The actual file offset,  in bytes,  is going to be 24 bytes
per term (each term consumes three double-precision floats) plus
the 290 bytes for the header.  That's why the actual 'fseek' call
in 'big_vsop.cpp' reads as

      fseek( ifile, 290L + (long)loc[0] * 24L, SEEK_SET);
*/

int DLL_FUNC calc_big_vsop_loc( FILE *ifile, const int planet,
                      double *ovals, double t, const double prec0)
{
   static int16_t cache[19];
   static int curr_planet = 99;
   int close_it = 0, value;

   ovals[0] = ovals[1] = ovals[2] = 0.;
   if( !planet)
      return( 0);       /* the sun */
   if( !ifile)
      {
      ifile = fopen( "big_vsop.bin", "rb");
      close_it = 1;
      }
   if( !ifile)
      return( -1);                              /* ...then give up. */
   if( curr_planet != planet)
      {                             /* reload the cache */
      fseek( ifile, (size_t)(planet - 1) * 6L * 3L * sizeof( int16_t), SEEK_SET);
      if( !fread( cache, 3 * 6 + 1, sizeof( int16_t), ifile))
         return( -2);
      curr_planet = planet;
      }

   t /= 10.;         /* convert to julian millennia */
   for( value = 0; value < 3; value++)
      {
      double sum, rval = 0., power = 1., prec = prec0;
      int16_t *loc = cache + value * 6;
      int i, j;

      fseek( ifile, 290L + (size_t)loc[0] * 24L, SEEK_SET);
      if( prec < 0.)
         prec = -prec;

      for( i = 6; i; i--, loc++)
         {
         double idata[3];

         sum = 0.;
         for( j = loc[1] - loc[0]; j; j--)
            {
            if( !fread( idata, 3, sizeof( double), ifile))
               return( -3);
            if( idata[0] > prec || idata[0] < -prec)
               {
               double argument = idata[1] + idata[2] * t;

               sum += idata[0] * cos( argument);
               }
            }
         rval += sum * power;
         power *= t;
         if( t)
            prec /= t;
         }
      ovals[value] = rval;
      }

   if( close_it)
      fclose( ifile);

   ovals[0] = fmod( ovals[0], 2. * PI);
   if( ovals[0] < 0.)
      ovals[0] += 2. * PI;
   return( 0);
}
