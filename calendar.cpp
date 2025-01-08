/* calendar.cpp: produces calendars in PostScript format
(probably not very useful to anybody but its author)

Copyright (C) 2013, Project Pluto

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

/* 2013 Jun 30:  added GPL.  Lunar phase reading assumed 32-bit
long integers;  replaced lots of "long"s with "int32_t"s. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"
#include "stringex.h"

#define LEGAL_SIZE (72 * 14)
#define US_LETTER_SIZE (72 * 11)

static int height = 72 * 17 / 2;
static int width = US_LETTER_SIZE;
         /* default is 8.5 x 11 inches */
static int xsize = 105;

#define X0 20
#define XEND (X0 + 7 * xsize)
#define TEXT_XOFFSET 3

#define Y0 20
#define YSIZE 104
#define YEND (Y0 + 5 * YSIZE)
#define TEXT_YOFFSET 15

#define DAY_OF_WEEK_HT 20
#define TOP_OF_DAYOFWK (Y0 + 5 * YSIZE + DAY_OF_WEEK_HT)

static const char *months[12] =  { "January", "February", "March",
               "April", "May", "June", "July", "August", "September",
               "October", "November", "December" };
static int show_jd_values = 0, dollhouse = 0, single_page = 0;

#define SHOW_JD    1
#define SHOW_MJD   2

      /* Kate prefers calendars without lunar phases or DST data: */
static int kate_style = 0;
      /* By default,  if the 24th and 31st land on Sundays or Mondays, */
      /* they share a cell;  similarly if the 23rd and 30th land on    */
      /* Sundays.  Use the '-w' switch,  and the 30th and/or 31st will */
      /* be moved to a cell of its own on the top line.                */
static int wraparound_style = 0;

#ifdef RESTORE_EASTER
void easter_date( const long year, int *month, int *day);
#endif

FILE *ofile = stdout;

static void show_month_text( const int month, const int year)
{
   char buff[80];

   snprintf_err( buff, sizeof( buff), "%s %d", months[month - 1], year);
   if( show_jd_values)
      {
      const long jd = dmy_to_day( 0, month, (long)year, CALENDAR_JULIAN_GREGORIAN) - 1;

      if( show_jd_values == SHOW_JD)
         snprintf_append( buff, sizeof( buff), " (JD %ld.5)", jd);
      else
         snprintf_append( buff, sizeof( buff), " (MJD %ld)", jd - 2400000L);
      }

   fprintf( ofile, "%d %d moveto (%s) show\n",
            X0 + 7 * xsize / 2 - (int)strlen( buff) * 20 / 3,
            TOP_OF_DAYOFWK + 5, buff);
}

static int get_phase_data( const double t, double *t_phases)
{
#ifdef _WIN32
   const char *phase_file_name = "z:\\vsop\\phases.dat";
#else
   const char *phase_file_name = "/home/phred/guide9/vsop/phases.dat";
#endif
   FILE *phase_file;
   int rval = -1;

   phase_file = fopen( phase_file_name, "rb");
   if( !phase_file)
      phase_file = fopen( "phases.dat", "rb");
   if( phase_file)
      {
      int32_t dt[12];
      long offset;
      double k;
      int i;
      int32_t phase0;

      k = floor( (t - 2451550.09765) / 29.530588853 - .5);
      if( !fread( &phase0, 1, sizeof( int32_t), phase_file))
         fprintf( stderr, "Error reading phase file\n");
      offset = (4L * ((long)k - phase0) + 1L) * (long)sizeof( int32_t);
      if( !fseek( phase_file, offset, SEEK_SET))
         {
         if( !fread( dt, 12, sizeof( int32_t), phase_file))
            {
            fprintf( stderr, "Read error on phase file\n");
            fprintf( stderr, "phase0 = %ld; offset %ld\n", (long)phase0, offset);
            exit( -1);
            }
         rval = 0;

         for( i = 0; i < 12; i++, k += .25)
            {
            const double t_phase1 = ( 2451550.09765 + 29.530588853 * k);

            t_phases[i] = t_phase1
                           + ((double)dt[i] - td_minus_utc( t_phase1)) / seconds_per_day;
            }
         }
      fclose( phase_file);
      }
   return( rval);
}

static void show_small_month( int month, int year, const int x0, int y0)
{
   long jd1, jd2;
   int n_days, starting_loc, i, j, k, day;
   char buff[100];
   double t_phases[12];

   while( month <= 0)
      {
      month += 12;
      year--;
      }
   while( month >= 13)
      {
      month -= 12;
      year++;
      }
   jd1 = dmy_to_day( 0, month, (long)year, CALENDAR_JULIAN_GREGORIAN);
   jd2 = dmy_to_day( 0, month + 1, (long)year, CALENDAR_JULIAN_GREGORIAN);
   n_days = (int)( jd2 - jd1);
   starting_loc = (int)( (jd1 + 2L) % 7L);
   strcpy( buff, months[month - 1]);
   snprintf_append( buff, sizeof( buff), " %d", year);
   y0 -= 12;
   if( get_phase_data( jd1 - 10, t_phases))
      for( i = 0; i < 12; i++)
         t_phases[i] = 0.;
   fprintf( ofile, "%d %d moveto (%s) show\n",
            x0 + (28 - (int)strlen( buff)) * TEXT_XOFFSET / 2, y0, buff);
   for( i = 0; i < 6; i++)
      {
      y0 -= 12;
      memset( buff, ' ', 20);
      for( j = 0; j < 7; j++)
         {
         day = i * 7 + j - starting_loc + 1;
         if( day >= 1 && day <= n_days)
            {
            char *tptr = buff + j * 3;
            const char *phase_text = "NM 1Q FM 3Q";

            tptr[0] = (char)( '0' + day / 10);
            tptr[1] = (char)( '0' + day % 10);
            if( day < 10)
               tptr[0] = ' ';
            for( k = 0; k < 12; k++)
               if( (long)(t_phases[k] + 0.5) == jd1 + day)
                  memcpy( tptr, phase_text + (k % 4) * 3, 2);
            }
         }
      for( j = 20; j && buff[j - 1] == ' '; j--)
         ;
      buff[j] = buff[20] ='\0';
      if( j)
         fprintf( ofile, "%d %d moveto (%s) show\n", x0, y0, buff);
      }
   if( show_jd_values)
      {
      jd1--;
      fprintf( ofile, "%d %d moveto ", x0, y0 - 12);
      if( show_jd_values == SHOW_JD)
         fprintf( ofile, "((JD %ld.5)) show\n", jd1);
      else
         fprintf( ofile, "((MJD %ld)) show\n", jd1 - 2400000L);
      }
}

static void make_postscript_substitutions( char *buff)
{
   while( *buff)
      {
      if( *buff == '(' || *buff == ')')
         {
         memmove( buff + 4, buff + 1, strlen( buff));
         buff[3] = (*buff == '(' ? '0' : '1');
         *buff++ = '\\';
         *buff++ = '0';
         *buff++ = '5';
         }
      buff++;
      }
}

   /* Dates will be drawn from 'dates.txt' for all years,  to cover */
   /* recurring events;  and from 'date####.txt' for events in a    */
   /* specific year.                                                */

#define MAX_DATES_PER_MONTH    300

static char **grab_dates( const int month, const int year)
{
   char **rval = (char **)calloc( MAX_DATES_PER_MONTH, sizeof( char *));
   char file_for_year[40];
   unsigned n_found = 0, i, j, pass;
   long jd;
   double t_phases[12];

   if( !rval)
      return( rval);
   snprintf_err( file_for_year, sizeof( file_for_year), "date%d.txt", year);
   for( pass = 0; pass < (kate_style ? 1u : 2u); pass++)
      {
      const char *date_filename = (kate_style ? "katedate.txt" : "dates.txt");
      FILE *dates_file = fopen( pass ? file_for_year : date_filename, "rb");

      if( dates_file)
         {
         char buff[80];

         while( fgets( buff, sizeof( buff), dates_file))
            if( !memcmp( months[month - 1], buff, 3))
               {
               unsigned loc = 7, n_days = 1;
               unsigned days[10];

               days[0] = (unsigned)atoi( buff + 3);
               if( buff[6] == '/')     /* See 'dates.txt' for an explanation */
                  {
                  while( dmy_to_day( (int)days[0] + 1, month, year,
                              CALENDAR_JULIAN_GREGORIAN) % 7 != buff[7] - '0')
                     (days[0])++;
                  loc = 9;
                  }
               while( buff[loc - 1] == ',')
                  {
                  int bytes_scanned;

                  sscanf( buff + loc, "%u%n", &days[n_days], &bytes_scanned);
                  loc += (unsigned)bytes_scanned + 1;
                  n_days++;
                  }
               make_postscript_substitutions( buff);
               for( i = 0; buff[i] >= ' '; i++)
                  ;
               buff[i] = '\0';
               for( j = 0; j < n_days; j++)
                  {
                  rval[n_found] = (char *)malloc( i - 4);
                  rval[n_found][0] = (char)days[j];
                  strcpy( rval[n_found] + 1, buff + loc);
                  n_found++;
                  assert( n_found < MAX_DATES_PER_MONTH);
                  }
               }
         fclose( dates_file);
         }
      }
            /* Reverse the order of the input: */
   if( n_found)
      for( i = 0, j = n_found - 1; i < j; i++, j--)
         {
         char *tptr = rval[i];

         rval[i] = rval[j];
         rval[j] = tptr;
         }

#ifdef RESTORE_EASTER
   for( i = 0; i < 3; i++)
#else
   for( i = 0; i < 2; i++)
#endif
      {
      static const char *names[3] = { "Rosh Hashanah", "Yom Kippur", "Easter" };
      int tmonth, day;

#ifdef RESTORE_EASTER
      if( i == 2)        /* easter */
         easter_date( (long)year, &tmonth, &day);
      else            /* Rosh Hashanah/Yom Kippur */
#endif
         {
         long tyear;

         jd = dmy_to_day( (i ? 10 : 1), 1,
                             (long)year + 3761L, CALENDAR_HEBREW);
         day_to_dmy( jd, &day, &tmonth, &tyear, CALENDAR_JULIAN_GREGORIAN);
         }
      if( tmonth == month)
         {
         rval[n_found] = (char *)malloc( strlen( names[i]) + 2);
         if( rval[n_found])
            {
            rval[n_found][0] = (char)day;
            strcpy( rval[n_found] + 1, names[i]);
            n_found++;
            }
         }
      }

#ifdef OBSOLETE
            /* 2006 and before:  DST change times are first Sunday in */
            /* April,  last Sun in October.  2007 and later:  second  */
            /* Sunday in March,  first Sunday in November.            */
   if( year < 2007 && (month == 4 || month == 10))
      {
      const unsigned day = year + year / 4 - year / 100 + year / 400;

      if( month == 4)
         dst_day = 7 - (day + 5) % 7;
      else
         dst_day = 31 - (day + 2) % 7;
      }
   if( year > 2006 && (month == 3 || month == 11))  /* incl start/end of D"S"T */
      {
      const unsigned day = year + year / 4 - year / 100 + year / 400;

      if( month == 3)
         dst_day = 14 - (day + 2) % 7;
      else
         dst_day = 7 - (day + 2) % 7;
      }
   if( dst_day && !kate_style)
      {
      rval[n_found] = (char *)malloc( 20);
      rval[n_found][0] = (char)dst_day;
      strcpy( rval[n_found] + 1, (month < 6 ? "D'S'T begins" : "D'S'T ends"));
      n_found++;
      }
#endif   // #ifdef OBSOLETE

   jd = dmy_to_day( 0, month, year, CALENDAR_JULIAN_GREGORIAN);
   if( !get_phase_data( jd - 10, t_phases) && !kate_style)
      for( i = 0; i < 12; i++)
         {
         unsigned day;
         const long minutes = (long)( (t_phases[i] + .5 - jd)
                                 * (double)minutes_per_day);

         day = (unsigned)( minutes / minutes_per_day);
         if( day > 0 && day < 32)
            {
            const char *phases[4] = { "*New moon", "*First quarter",
                       "*Full moon", "*Last quarter" };
            const char *phasestr = phases[i & 3];
            char *tptr;

            rval[n_found] = tptr = (char *)malloc( strlen( phasestr) + 20);
            if( tptr)
               {
               *tptr++ = (char)day;
               strcpy( tptr, phasestr);
               n_found++;
               }
            }
         }
   return( rval);
}

static const char * const trailer_data[51] = {
                "%%Page: 1 1",
                "%%PageOrientation: Landscape",
                "gsave height width translate 180 rotate calendar grestore",
                "",
                "%%Page: 2 2",
                "%%PageOrientation: Portrait",
                "  gsave",
                "  /doublescale 1.34 def",
                "  90 rotate 0 -450 doublescale mul translate",
                "  doublescale doublescale scale",
                "  calendar grestore",
                "",
                "%%Page: 3 3",
                "%%PageOrientation: Portrait",
                "  gsave",
                "  /doublescale 1.34 def",
                "  90 rotate 0 -850 doublescale mul translate",
                "  doublescale doublescale scale",
                "  calendar grestore",
                "",
                "%%Page: 4 4",
                "%%PageOrientation: Landscape",
                "  gsave",
                "  /quadscale 2 def",
                "  quadscale quadscale scale",
                "  calendar grestore",
                "",
                "%%Page: 5 5",
                "%%PageOrientation: Landscape",
                "  gsave",
                "  /quadscale 2 def",
                "  0 -370 quadscale mul translate",
                "  quadscale quadscale scale",
                "  calendar grestore",
                "",
                "%%Page: 6 6",
                "%%PageOrientation: Landscape",
                "  gsave",
                "  /quadscale 2 def",
                "  -290 quadscale mul -370 quadscale mul translate",
                "  quadscale quadscale scale",
                "  calendar grestore",
                "",
                "%%Page: 7 7",
                "%%PageOrientation: Landscape",
                "  gsave",
                "  /quadscale 2 def",
                "  -290 quadscale mul 0 quadscale mul translate",
                "  quadscale quadscale scale",
                "  calendar grestore",
                NULL };


static const char *dollhouse_trailer_data[] = {
                "%%Page: 1 1",
                "%%PageOrientation: Landscape",
                "  gsave",
                "  /shrinkscale .23 def",
                "  shrinkscale shrinkscale scale",
                "  130 2450 translate calendar_1",
                "  800 0 translate calendar_2",
                "  800 0 translate calendar_3",
                "  -1600 -690 translate calendar_4",
                "  800 0 translate calendar_5",
                "  800 0 translate calendar_6",
                "  -1600 -690 translate calendar_7",
                "  800 0 translate calendar_8",
                "  800 0 translate calendar_9",
                "  -1600 -690 translate calendar_10",
                "  800 0 translate calendar_11",
                "  800 0 translate calendar_12",
                " grestore showpage",
                NULL };


#define FONT_UNSET -1
#define FONT_ITALIC 1
#define FONT_PLAIN  2
#define FONT_BOLD   3

static int calendar( const int month, const int year)
{
   int i, j;
   const long jd1 = dmy_to_day( 0, month, (long)year, CALENDAR_JULIAN_GREGORIAN);
   const long jd2 = dmy_to_day( 0, month + 1, (long)year, CALENDAR_JULIAN_GREGORIAN);
   const int n_days = (int)( jd2 - jd1), starting_loc = (int)( (jd1 + 1L) % 7L);
   char buff[100];
   char lines_used[35], phases_shown[35];
   char **dates = grab_dates( month, year);

   if( month == 1 || !dollhouse)
      {
      fprintf( ofile, "%%!PS-Adobe-2.0\n");
      fprintf( ofile, "%%%%Pages: %d\n", (dollhouse || single_page ? 1 : 7));
      fprintf( ofile, "%%%%PageOrder: Ascend\n");
      fprintf( ofile, "%%%%DocumentMedia: Default 612 %d 0 () ()\n", width);
      fprintf( ofile, "%%%%Creator: calendar.cpp\n");
      fprintf( ofile, "%%%%Copyright: none\n");
      fprintf( ofile, "%%%%Title: Calendar for %s %d\n", months[month - 1], year);
      fprintf( ofile, "%%%%Version: none\n");
      fprintf( ofile, "%%%%DocumentData: Clean7Bit\n");
      fprintf( ofile, "%%%%EndComments\n");
      fprintf( ofile, "%%%%BeginDefaults\n");
      fprintf( ofile, "%%%%PageResources: font Times-Roman\n");
      fprintf( ofile, "%%%%PageResources: font Times-Italic\n");
      fprintf( ofile, "%%%%PageResources: font Courier-Bold\n");
      fprintf( ofile, "%%%%EndDefaults\n");
      fprintf( ofile, "\n/width %d def\n", width);
      fprintf( ofile, "/height %d def\n\n", height);
      fprintf( ofile, "/blood {\n");
      fprintf( ofile, "/aa 3 def /bb 2.5 def\n");
      fprintf( ofile, "aa 0 rmoveto bb 0 rlineto 0 aa rlineto aa 0 rlineto 0 bb rlineto\n");
      fprintf( ofile, "0 aa sub 0 rlineto 0 aa rlineto 0 bb sub 0 rlineto\n");
      fprintf( ofile, "0 0 aa sub rlineto 0 aa sub 0 rlineto\n");
      fprintf( ofile, "0 0 bb sub rlineto aa 0 rlineto 0 0 aa sub rlineto aa 3 mul 0 rmoveto\n");
      fprintf( ofile, "} def\n\n");
      }
   if( dollhouse)
      fprintf( ofile, "/calendar_%d {\n", month);
   else
      {
      fprintf( ofile, "/calendar {\n");
      fprintf( ofile, "0 %d translate -90 rotate\n", width);
      }

   for( i = 0; i <= 5; i++)             /* horizontal lines separating weeks */
      fprintf( ofile, "%d %d moveto %d %d lineto\n", X0, Y0 + i * YSIZE,
                                           XEND, Y0 + i * YSIZE);
   fprintf( ofile, "%d %d moveto %d %d lineto\n", X0, TOP_OF_DAYOFWK,
                                        XEND, TOP_OF_DAYOFWK);

   for( i = 0; i <= 7; i++)       /* vertical lines */
      fprintf( ofile, "%d %d moveto %d %d lineto\n", X0 + i * xsize, Y0,
                                             X0 + i * xsize, YEND);
   fprintf( ofile, "/defaultfontsize  { 14 scalefont } def\n");
   fprintf( ofile, "/itemfontsize  { 12 scalefont } def\n");
   fprintf( ofile, "/Times-Roman findfont defaultfontsize setfont\n");
   for( i = 0; i < 35; i++)
      lines_used[i] = phases_shown[i] = 0;

   for( i = n_days; i >= 1; i--)
      {
      int xcell = (i + starting_loc) % 7, curr_font = FONT_UNSET;
      int ycell = 4 - (i + starting_loc) / 7, double_cell = 0, x0, y0;
      int cell_idx;

      if( starting_loc == 6)        /* first of month is on a sunday */
         ycell++;
      if( !wraparound_style && !ycell
                      && i + 7 <= n_days)     /* two days in one cell */
         {
         snprintf_err( buff, 6, "%d/%d", i, i + 7);
         double_cell = 1;
         }
      else
         snprintf_err( buff, 3, "%d", i);
      if( wraparound_style && ycell < 0)    /* this is a wraparound case */
         ycell = 4;
      x0 = X0 + xcell * xsize;
      y0 = Y0 + (ycell + 1) * YSIZE;
      if( ycell >= 0)
         {
         const int xboxsize = TEXT_XOFFSET + 8 * (int)strlen( buff);

         fprintf( ofile, "%d %d moveto (%s) show\n",
                     x0 + TEXT_XOFFSET, y0 - TEXT_YOFFSET, buff);
         fprintf( ofile, "%d %d moveto %d %d rlineto %d %d rlineto\n",
               x0, y0 - (TEXT_YOFFSET + 2),
               xboxsize, 0, 0, TEXT_YOFFSET + 2);
         }
      else
         {
         double_cell = 1;
         ycell = 0;
         y0 = Y0 + YSIZE;
         }
      cell_idx = (4 - ycell) * 7 + xcell;
      for( j = 0; dates[j]; j++)
         if( dates[j][0] == i)
            {
            char *text_to_show = dates[j] + 1;
            const int is_lunar_phase = (*text_to_show == '*');
            int font_to_use = (is_lunar_phase ? FONT_ITALIC : FONT_PLAIN);
            int xloc = X0 + xcell * xsize + TEXT_XOFFSET +
                           ((is_lunar_phase && !double_cell) ? 24 : 0);
            int yloc = (is_lunar_phase ?
                         y0 - TEXT_YOFFSET * (1 + double_cell)
                                        - phases_shown[cell_idx] * 11 :
                         y0 - YSIZE + 3 + lines_used[cell_idx] * 11);

            if( *text_to_show == '^')
               font_to_use = FONT_BOLD;
            if( curr_font != font_to_use)
               {
               curr_font = font_to_use;
               if( curr_font == FONT_ITALIC)
                  fprintf( ofile, "/Times-Italic findfont 14 scalefont setfont\n");
               if( curr_font == FONT_PLAIN)
                  fprintf( ofile, "/Times-Roman findfont itemfontsize setfont\n");
               if( curr_font == FONT_BOLD)
                  fprintf( ofile, "/Times-Bold findfont itemfontsize setfont\n");
               }
            if( curr_font != FONT_PLAIN)
               text_to_show++;    /* first char is a font indicator */

            fprintf( ofile, "%d %d moveto ", xloc, yloc);
            if( double_cell && *text_to_show)
               fprintf( ofile, "(\050%d\051 ) show ", i);
            if( *text_to_show == '+' && text_to_show[1] == ' ')
               {
               fprintf( ofile, "blood ");
               text_to_show += 2;
               }
            if( *text_to_show == '\\' && text_to_show[1] == 'p')
               {
               fprintf( ofile, "/Symbol findfont itemfontsize setfont (p) show\n");
               fprintf( ofile, "/Times-Roman findfont itemfontsize setfont\n");
               text_to_show += 2;
               }
            fprintf( ofile, "(%s) show\n", text_to_show);
            if( is_lunar_phase)
               phases_shown[cell_idx]++;
            else
               lines_used[cell_idx]++;
            }
      if( !lines_used[cell_idx])         /* mark it as being a "used" cell */
         lines_used[cell_idx] = 1;
      if( curr_font != FONT_UNSET)
         fprintf( ofile, "/Times-Roman findfont defaultfontsize setfont\n");
      }

   for( i = 0; i < 7; i++)
      {
      const char *day_of_week_text[7] = { "Sunday", "Monday", "Tuesday",
                           "Wednesday", "Thursday", "Friday", "Saturday" };

      fprintf( ofile, "%d %d moveto (%s) show\n",
                              X0 + i * xsize + xsize / 2 -
                              (int)strlen( day_of_week_text[i]) * 4,
                              Y0 + 5 * YSIZE + 5,
                              day_of_week_text[i]);
      }

   fprintf( ofile, "/Times-Roman findfont 24 scalefont setfont\n");

   show_month_text( month, year);

   fprintf( ofile, "/Courier-Bold findfont %d scalefont setfont\n",
            (width == LEGAL_SIZE ? 10 : 8));

               /* show two preceding & two (or more) following months: */
   for( i = j = 0; i < 35; i++)
      if( !lines_used[i])
         {
         show_small_month( month + j - 2, year,
                     X0 + (i % 7) * xsize + 1, Y0 + YSIZE * (5 - i / 7));
         j++;
         if( j == 2)    /* avoid duplicating the 'main' month */
            j = 3;
         }
   fprintf( ofile, "stroke %s } def\n", (dollhouse ? "" : "showpage"));
   if( !dollhouse)
      for( i = 0; (single_page ? i < 3 : trailer_data[i] != NULL); i++)
         fprintf( ofile, "%s\n", trailer_data[i]);
   return( 0);
}

static void error_exit( const int error_code)
{
   fprintf( stdout, "'calendar' needs a year and month on the command line.\n"
                    "It will produce a PostScript (R) calendar.  Given only\n"
                    "the year,  it will produce a twelve-month calendar.\n\n"
                    "Options are :\n\n"
                    "   -d  'dollhouse' style (small,  full year on one page)\n"
                    "   -s  single page,  don't create large versions\n"
                    "   -j  show JD values\n"
                    "   -m  show MJD values\n"
                    "   -l  create for US legal-size paper\n"
                    "   -w  wraparound style;  no 'divided' dates\n"
                    "   -o(filename)  specify output file;  default is stdout\n");
   exit( error_code);
}

int main( const int argc, const char **argv)
{
   int i, month = 0, year;

   if( argc == 1)
      error_exit( -2);
   year = atoi( argv[1]);
   month = atoi( argv[2]);
   if( month > year)       /* allow for possible reversed entry */
      {
      month = year;
      year = atoi( argv[2]);
      }
   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'd':
               dollhouse = 1;
               break;
            case 'j':
               show_jd_values = SHOW_JD;
               break;
            case 's':
               single_page = 1;
               break;
            case 'k':
               kate_style = 1;
               break;
            case 'l':            /* legal size */
               width = LEGAL_SIZE;
               fprintf( stderr, "Use Atril;  specify US legal-size\n");
               break;
            case 'm':
               show_jd_values = SHOW_MJD;
               break;
            case 'o':
               ofile = fopen( argv[i] + 2, "wb");
               assert( ofile);
               break;
            case 'w':
               wraparound_style = 1;
               break;
            default:
               printf( "Unknown option '%s'\n", argv[i]);
               return( -1);
            }
   if( !dollhouse && (month < 1 || month > 12))
      {
      fprintf( stderr, "Can't do year %d, month %d\n", year, month);
      error_exit( -1);
      }
   xsize = (width - 60) / 7;
   if( !dollhouse)
      return( calendar( month, year));
   else
      {
      for( i = 1; i <= 12; i++)
         calendar( i, year);
      for( i = 0; dollhouse_trailer_data[i]; i++)
         fprintf( ofile, "%s\n", dollhouse_trailer_data[i]);
      }
   return( 0);
}
