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

#include <stdio.h>
#include <stdlib.h>
#include "lun_tran.h"

int main( const int argc, const char **argv)
{
   double lon, lat;
   int day, zip_code, year, month;
   int time_zone, use_dst;
   char place_name[100];

   if( argc != 4)
      {
      printf( "'lun_test' requires,  as command-line arguments,  a year,\n");
      printf( "month, and US five-digit postal code.  For example,\n\n");
      printf( "./luntest 2015 2 04008\n\n");
      printf( "would produce a list of lunar transit times for 2015 February\n");
      printf( "for the town of Bowdoinham, Maine,  which has ZIP code 04008.\n\n");
      return( -2);
      }
   year = atoi( argv[1]);
   month = atoi( argv[2]);
   zip_code = atoi( argv[3]);
   if( get_zip_code_data( zip_code, &lat, &lon, &time_zone,
                           &use_dst, place_name))
      {
      printf( "Couldn't get data for ZIP code %05d\n", zip_code);
      return( -1);
      }

   printf( "Lat %f   Lon %f   Time zone %d  DST=%d  %s\n",
         lat, lon, time_zone, use_dst, place_name);
   for( day = 1; day <= 31; day++)
      {
      const double transit_time =
            get_lunar_transit_time( year, month, day,
                  lat, lon, time_zone, use_dst, 1);
      const double antitransit_time =
            get_lunar_transit_time( year, month, day,
                  lat, lon, time_zone, use_dst, 0);
      char transit_buff[6], antitransit_buff[6];

      format_hh_mm( transit_buff, transit_time);
      format_hh_mm( antitransit_buff, antitransit_time);

      printf( "%2d: %s %s\n", day,
            transit_buff, antitransit_buff);
      }
   return( 0);
}
