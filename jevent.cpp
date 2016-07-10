/* jevent.cpp: computes Galilean satellite events

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

/* The following code is an early and somewhat embarrassing effort.
Please do not judge my programming abilities by it.

   It computes the date/times of "classical" Jovian satellite events
(transits,  shadows,  occultations,  eclipses) over a given time
span,  and sorts them out to produce the sort of event list one
sees in _Sky and Telescope_ and similar publications.  It can
also store them in a binary form that is used within Guide,  so
that that program can show Galilean events without having to
do all the math itself.   */

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "watdefs.h"
#include "lunar.h"
#include "afuncs.h"
#include "date.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define JUPITER_DIAMETER_IN_KM (88000. * 1.609)
#define AU_PER_JRAD (JUPITER_DIAMETER_IN_KM / AU_IN_KM)
#define EVENT struct event
#define VSOP_CHUNK 2767U

EVENT
   {
   double t;
   unsigned sat, event_type;
   };

/* We don't actually explicitly use these four macros :
#define OCCULTA  0x00
#define TRANSIT  0x01
#define ECLIPSE  0x02
#define SHADOW   0x03
*/
#define IN_FRONT  0x01
#define FROM_SUN 0x02

#define EVENT_START   0x04
#define EVENT_UNSEEN  0x08

const double speeds[4] = {203.488955432, 101.374724550, 50.317609110, 21.571071314};
static char FAR *vsop_data;
static int quiet = 0;

static void show_event( FILE *ofile, EVENT *e)
{
   const char *event_str[4] = {"Occ", "Tra", "Ecl", "Sha"};
   char buff[80];
   const double ut_jd = e->t - td_minus_ut( e->t) / seconds_per_day;

   if( e->event_type & EVENT_UNSEEN)
      return;
   full_ctime( buff, ut_jd, FULL_CTIME_FORMAT_HH_MM);
   fprintf( ofile, "sat %u: %s %s: %s\n", e->sat, event_str[e->event_type & 3],
              (e->event_type & EVENT_START) ? "start" : "end  ", buff);
}

static unsigned find_events( unsigned sat_no, double t1, double t2, int viewpoint, EVENT *e)
{
   double t, lon_j, lat_j, rad_j, lon_e, lat_e, rad_e;
   double loc[9], tloc[18], *tptr, step, delta_lat, prev_delta_lat = 0.;
   double tc, lon_s, lat_s, rad_s, prev_delta = 0., delta;
   int i;
   unsigned rval = 0;

   t = t1;
   step = 10. / speeds[sat_no - 1];
   while( t < t2)
      {
      tc = (t - 2451545.) / 36525.;     /* re-cvt to julian centuries */

      lon_j = calc_vsop_loc( vsop_data, 5, 0, tc, 0.);
      lat_j = calc_vsop_loc( vsop_data, 5, 1, tc, 0.);
      rad_j = calc_vsop_loc( vsop_data, 5, 2, tc, 0.);

      lon_e = calc_vsop_loc( vsop_data, 3, 0, tc, 0.);
      lat_e = calc_vsop_loc( vsop_data, 3, 1, tc, 0.);
      rad_e = calc_vsop_loc( vsop_data, 3, 2, tc, 0.);

      loc[0] = rad_j * cos( lat_j) * cos( lon_j);
      loc[1] = rad_j * cos( lat_j) * sin( lon_j);
      loc[2] = rad_j * sin( lat_j);

      loc[6] = rad_e * cos( lat_e) * cos( lon_e);
      loc[7] = rad_e * cos( lat_e) * sin( lon_e);
      loc[8] = rad_e * sin( lat_e);

      calc_jsat_loc( t, tloc, 1 << (sat_no - 1), 0L);
      tptr = tloc + (sat_no - 1) * 3;
      loc[3] = loc[0] + tptr[0] * AU_PER_JRAD;
      loc[4] = loc[1] + tptr[1] * AU_PER_JRAD;
      loc[5] = loc[2] + tptr[2] * AU_PER_JRAD;

      if( viewpoint)
         {
         for( i = 0; i < 9; i++)
            loc[i] -= loc[6 + i % 3];
         lon_j = atan2( loc[1], loc[0]);
         lat_j = atan( loc[2] / sqrt( loc[0] * loc[0] + loc[1] * loc[1]));
         }
      rad_s = sqrt( loc[3] * loc[3] + loc[4] * loc[4] + loc[5] * loc[5]);
      lon_s = atan2( loc[4], loc[3]);
      while( lon_s - lon_j > PI)
         lon_s -= PI + PI;
      while( lon_s - lon_j <-PI)
         lon_s += PI + PI;
      lat_s = asin( loc[5] / rad_s);
      delta = (lon_s - lon_j) * rad_s / AU_PER_JRAD;
      delta_lat = (lat_s - lat_j) * rad_s / AU_PER_JRAD;
      delta_lat *= 1.071374;     /* stretch for jup's oblateness */
      if( delta * prev_delta < 0.)     /* zero crossed */
         {
         double t_crossing, diff, a, b, c, dx, dy, dist = 0.;

         dx = prev_delta - delta;
         dy = prev_delta_lat - delta_lat;
         a = dx * dx + dy * dy;              /* quadratic for intercept */
         b = 2. * (dx * delta + dy * delta_lat);
         c = delta_lat * delta_lat + delta * delta - 1.;
         diff = b * b - 4. * a * c;
         t_crossing = t + b * step / (2. * a);
         t = t_crossing + (180. / speeds[sat_no - 1]) * .9;
         if( diff > 0.)       /* if real solution, ie, sat doesn't miss */
            {
            diff = sqrt( diff) * step / (2. * a);
            diff = fabs( diff);
            for( i = 0; i < 3; i++)
               dist += (loc[i + 6] - loc[i]) * (loc[i + 6] - loc[i]);
            t_crossing += sqrt( dist) / AU_PER_DAY;
            e[0].t = t_crossing - diff;
            e[0].sat = sat_no;
            e[0].event_type = (prev_delta > 0.);
            if( !viewpoint)
               e[0].event_type |= FROM_SUN;
            e[1].t = t_crossing + diff;
            e[1].sat = sat_no;
            e[1].event_type = e[0].event_type;
            e[0].event_type |= EVENT_START;
            if( !quiet)
               {
               show_event( stdout, e);
               show_event( stdout, e + 1);
               }
            e += 2;
            rval += 2;
            }
         delta = 0.;
         }
      prev_delta = delta;
      prev_delta_lat = delta_lat;
      t += step;
      }
   return( rval);
}

int main( int argc, char **argv)
{
   unsigned i, j, k, julian = 0, gap;
   unsigned n_days = 30, sat_no = 15;
   unsigned max_events;
   unsigned n_events = 0, n_sun, n_earth;
   double t1, t2;
   long jd;
   EVENT *e;
   FILE *ofile = NULL, *data_file = NULL;
   FILE *vsop_file;
   char *vsop_tbuff;

   vsop_file = fopen( "vsop.bin", "rb");
   if( !vsop_file)
      {
      printf( "Couldn't open vsop.bin");
      return( -3);
      }
   vsop_tbuff = (char *)malloc( VSOP_CHUNK);
   vsop_data = (char FAR *)FMALLOC( VSOP_CHUNK * 22U);
   for( i = 0; i < 22; i++)
      {
      if( !fread( vsop_tbuff, VSOP_CHUNK, 1, vsop_file))
         {
         printf( "Couldn't read VSOP data\n");
         free( vsop_tbuff);
         free( vsop_data);
         return( -1);
         }
      memcpy( vsop_data + (unsigned)i * VSOP_CHUNK, vsop_tbuff, VSOP_CHUNK);
      }
   fclose( vsop_file);
   free( vsop_tbuff);

   if( argc < 4)
      {
      printf( "JEVENT calculates Jovian satellite events (transits, eclipses,\n");
      printf( "shadows,  occultations) as seen from Earth.  It requires, as a\n");
      printf( "minimum,  a day,  month,  and year.  Given,  say,  the command\n");
      printf( "\nJEVENT 18 4 1993\n\n");
      printf( "JEVENT will calculate all events from 18 Apr 1993 for the next\n");
      printf( "thirty days.  You can add on the following parameters:\n\n");
      printf( "   -j        Use Julian calendar\n");
      printf( "   -d(#)     Calculate for (#) days instead of 30\n");
      printf( "   -f(name)  Put results in ASCII file (name) as well as on screen\n");
      return( -2);
      }
   for( i = 0; i < (unsigned)argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'j': case 'J':
               julian = 1;
               break;
            case 'q': case 'Q':
               quiet = 1;
               break;
            case 's': case 'S':
               sat_no = (unsigned)atoi( argv[i] + 2);
               break;
            case 'd': case 'D':
               n_days = (unsigned)atoi( argv[i] + 2);
               break;
            case 'f': case 'F':
               ofile = fopen( argv[i] + 2, "wb");
               break;
            case 'r': case 'R':
               data_file = fopen( argv[i] + 2, "ab");
               printf( "Appending data to %s\n", argv[i] + 2);
               break;
            default:
               break;
            }
                     /* We seem to average about eight events/day,  before */
                     /* sorting and removing 'hidden' events.  But let's   */
                     /* allow a little margin :                            */
   max_events = n_days * 10 + 50;
   e = (EVENT *)calloc( max_events, sizeof( EVENT));
   if( !e)
      return( -1);
   jd = dmy_to_day( 0, atoi( argv[2]), atol( argv[3]), (int)julian);
   t1 = (double)jd - .5 + atof( argv[1]);
   t2 = t1 + (double)n_days;
   printf( "JD %f to %f\n", t1, t2);
   for( i = 0; i < 4; i++)
      if( sat_no & (1 << i))
         {
         printf( "Sat %u from sun\n", i + 1);
         n_sun = find_events( i + 1, t1 - 1., t2 + 1., 0, e + n_events);
         printf( "Sat %u from earth\n", i + 1);
         n_earth = find_events( i + 1, t1 - 1., t2 + 1., 1, e + n_events + n_sun);
         printf( "Finding hidden events\n");
         k = n_events + n_sun;
         for( j = n_events; j < n_events + n_sun; j++)
            while( k < n_events + n_sun + n_earth && e[k].t < e[j].t)
               {
               if( e[k + 1].t > e[j].t)
                  if( !(e[j].event_type & IN_FRONT))
                     e[j].event_type |= EVENT_UNSEEN;
               k += 2;
               }
         k = n_events;
         for( j = n_events + n_sun; j < n_events + n_sun + n_earth; j++)
            while( k < n_events + n_sun && e[k].t < e[j].t)
               {
               if( e[k + 1].t > e[j].t)
                  if( !(e[j].event_type & IN_FRONT))
                     e[j].event_type |= EVENT_UNSEEN;
               k += 2;
               }
         n_events += n_earth + n_sun;
         }
   printf( "Sorting %u events\n", n_events);
   for( gap = 1; gap < n_events / 3; gap = gap * 3 + 1)
      ;
   while( gap)
      {
      for( i = 0; i < gap; i++)
         for( j = i; j + gap < n_events; j += gap)
            if( e[j].t > e[j + gap].t)
               {
               EVENT temp;

               memcpy( &temp, e + j + gap, sizeof( EVENT));
               memcpy( e + j + gap, e + j, sizeof( EVENT));
               memcpy( e + j, &temp, sizeof( EVENT));
               if( j >= gap)
                  j -= gap + gap;
               }
      gap /= 3;
      }

   if( !quiet)
      {
      printf( "Final results:\n");
      for( i = 0; i < n_events; i++)
         if( e[i].t > t1 && e[i].t < t2)
            show_event( stdout, e + i);
      }

   if( ofile)
      {
      for( i = 0; i < n_events; i++)
         if( e[i].t > t1 && e[i].t < t2)
            show_event( ofile, e + i);
      fclose( ofile);
      }

   if( data_file)
      {
      for( i = 0; i < n_events; i++)
         if( e[i].t > t1 && e[i].t < t2)
            if( !( e[i].event_type & EVENT_UNSEEN))
               {
               char buff[5];
               int32_t tval;

                        /* store minutes from 2000.0 */

               tval = (int32_t)( (e[i].t - 2451545.0) * 1440.);
               memcpy( buff, &tval, sizeof( int32_t));
               buff[4] = (char)( (e[i].event_type & 7) | (e[i].sat << 3));
               fwrite( buff, 1, 5, data_file);
               }
      fclose( data_file);
      }
   return( 0);
}
