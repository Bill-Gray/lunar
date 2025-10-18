#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "watdefs.h"
#include "date.h"
#include "lunar.h"

/* Small program-let to generate the ASCII calendar with lunar phases
shown in interactive Find_Orb when you hit 'c'.  (See the file
'calendar.txt' in the Find_Orb repository.)  Compile with

gcc -Wextra -Wall -O3 -pedantic text_ca2.c -o text_ca2 liblunar.a -lm

   See 'boxize.c' in the 'junk' repo for code to turn the output from
+ - | dividers into Unicode box characters.

   make_base_line( ) will produce one of the following four lines,
for line = 0, 1, 2, 3,  except with the width depending on the value
of 'cols'.

0: (all spaces)
1: +---------+---------+---------+---------+---------+---------+---------+
2: |  Sunday |  Monday | Tuesday |Wednesday| Thursday|  Friday | Saturday|
3: |         |         |         |         |         |         |         |
*/

static void make_base_line( char *txt, const int cols, const int line)
{
   int i;

   memset( txt, (line != 1 ? ' ' : '-'), cols * 7);
   txt[cols * 7 + 1] = '\0';
   if( !line)
      return;
   for( i = 0; i < 8; i++)
      txt[i * cols] = (line != 1 ? '|' : '+');
   if( line == 2)
      for( i = 0; i < 7; i++)
         {
         const char *days[7] = { "Sunday", "Monday", "Tuesday", "Wednesday",
                        "Thursday", "Friday", "Saturday" };
         int len = (int)strlen( days[i]);

         if( len > cols - 1)
            len = cols - 1;
         memcpy( txt + i * cols + (cols - len + 1) / 2, days[i], len);
         }
}

static int lines_per_week = 4;

static void show_calendar( const long year, const int month, const int cols)
{
   int i, day, n_output_lines = 4 + 5 * lines_per_week;
   const long jd0 = dmy_to_day( 0, month, year, CALENDAR_GREGORIAN);
   const int loc0 = (int)( (jd0 + 2) % 7) - 1;
   char **text;
   const char *month_names[12] = { "January", "February", "March",
            "April", "May", "June", "July", "August", "September",
            "October", "November", "December" };

   text = (char **)calloc( n_output_lines, sizeof( char *));
   text[0] = (char *)malloc( n_output_lines * (7 * cols + 2));
   for( i = 1; i < n_output_lines; i++)
      text[i] = text[i - 1] + 7 * cols + 2;
   for( i = 0; i < 3; i++)
      make_base_line( text[i], cols, i);
   for( i = 0; i <= 5 * lines_per_week; i++)
      make_base_line( text[i + 3], cols, (i % lines_per_week ? 3 : 1));
   for( day = days_in_month( month, year, CALENDAR_GREGORIAN); day; day--)
      {
      const unsigned col = cols * ((unsigned)(day + loc0) % 7u) + 2u;
      unsigned line = (unsigned)((day + loc0) / 7u) % 5u;
      const long t2k = (double)jd0 + (double)day - 2451545.;
      int phase;

      line = line * (unsigned)lines_per_week + 4u;
      if( day >= 10)
         text[line][col] = '0' + day / 10u;
      text[line][col + 1] = '0' + day % 10u;
      for( phase = 0; phase < 4; phase++)
         if( (long)( find_nearest_lunar_phase_time( phase, (long double)t2k) + .5) == t2k)
            {
            const char *phases[4] = { "New Moon", "1st Qtr", "Full Moon", "3rd Qtr" };

            memcpy( text[line + 1] + col - 1, phases[phase], strlen( phases[phase]));
            }
      }
   sprintf( text[0] + (cols * 7 - strlen( month_names[month - 1])) / 2 - 2,
                             "%s %ld", month_names[month - 1], year);
   for( i = 0; i < 4 + 5 * lines_per_week; i++)
      printf( "%s\n", text[i]);
   free( text[0]);
   free( text);
   printf( "\n");
}

static void usage_exit( void)
{
   fputs( "usage: text_ca2 year n_years [-c cols] [-l lines]\n", stderr);
   exit( -1);
}

int main( const int argc, const char **argv)
{
   long year = (argc > 1 ? atoi( argv[1]) : 0);
   long n_years = (argc > 2 ? atoi( argv[2]) : 0);
   int cols = 11;

   if( !n_years || !year)
      usage_exit( );
   for( int i = 3; i < argc - 1; i++)
      if( argv[i][0] == '-')
         {
         switch( argv[i][1])
            {
            case 'l' :
               lines_per_week = atoi( argv[i + 1]);
               break;
            case 'c' :
               cols = atoi( argv[i + 1]);
               break;
            default:
               fprintf( stderr, "Unrecognized argument '%s'\n", argv[i]);
               usage_exit( );
               break;
            }
         i++;
         }
      else
         usage_exit( );
   while( n_years--)
      {
      int month;

      for( month = 1; month <= 12; month++)
         show_calendar( year, month, cols);
      year++;
      }
   return( 0);
}
