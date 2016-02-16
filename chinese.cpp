#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define YEAR_DATA struct year_data

#pragma pack(1)

YEAR_DATA
   {
   unsigned short mask;
   unsigned char intercalary_month, offset;
   } *ydata;

#pragma pack()

short n_years, year0;
unsigned short bitmask = 1u;
long prev_jd;

static void dump_data( char *text[], const int n_found)
{
   int i, month = 12, year = atoi( text[0] + 35), intercalary_month = 0;

   if( text[0][8] != 'z')
      printf( "NOT SOLSTICE\n");
   if( n_found < 24 || n_found > 25)
      {
      printf( "n_found = %d\n", n_found);
      exit( 0);
      }
   if( !year0)               /* Remember,  we get the last month of the */
      year0 = (short)( year + 1);    /* previous year... skip it */

   printf( "\n");
   for( i = 0; i < n_found; i++)
      {
      if( text[i][8] == ' ')        /* it's a month */
         {
         const long jd = atol( text[i]);
         const int prev_month = (month + 10) % 12 + 1;
         const int prev_year = year - (prev_month == 12);

         if( jd - prev_jd == 30L)
            if( prev_year >= year0 && prev_year < year0 + n_years)
               ydata[prev_year - year0].mask |= bitmask;
         bitmask <<= 1;

         if( month == 1 && intercalary_month != 1)
            {           /* it's the for-real first month of the year: */
            const long offset = jd - (year * 365L + year / 4L + 757861L);

            bitmask = 1;
            if( year >= year0 && year < year0 + n_years)
               ydata[year - year0].offset = (unsigned char)offset;
            }

         if( n_found == 25 && !intercalary_month && text[i + 1][8] == ' ')
            {
            sprintf( text[i] + 31, "%4di%5d\n", prev_month, prev_year);
            intercalary_month = prev_month;
            if( prev_year >= year0 && prev_year < year0 + n_years)
               ydata[prev_year - year0].intercalary_month =
                                                 (unsigned char)prev_month;
            }
         else           /* it was just a normal month... move on by one */
            {
            sprintf( text[i] + 31, "%4d%6d\n", month, year);
            month++;
            if( month == 13)
               {
               month = 1;
               year++;
               }
            }

         prev_jd = jd;
         }
      printf( "%s", text[i]);
      }
}

int main( const int argc, const char **argv)
{
   char *buff, *text[30], ibuff[60];
   int i, n_found = 0;
   FILE *ifile = fopen( argv[1], "rb");

   if( !ifile)
      {
      printf( "%s not opened\n", argv[1]);
      return( -1);
      }
   if( argc == 3)
      {
      n_years = (short)atoi( argv[2]);
      ydata = (YEAR_DATA *)calloc( n_years, sizeof( YEAR_DATA));
      printf( "Looking for %d years of data\n", n_years);
      }
   buff = (char *)malloc( 30 * 60);
   for( i = 0; i < 30; i++)
      text[i] = buff + i * 60;
   while( fgets( ibuff, 60, ifile))
      {
      strcpy( text[n_found], ibuff);
      if( ibuff[33] == '1' && ibuff[34] == '1')
         {
         if( n_found == 24 || n_found == 25)
            dump_data( text, n_found);
         n_found = 0;
         }
      strcpy( text[n_found++], ibuff);
      memset( ibuff, 0, 60);
      if( n_found > 25)
         {
         printf( "ERROR in input file\n");
         for( i = 0; i < n_found; i++)
            printf( "%s", text[i]);
         exit( 0);
         }
      }
   if( ydata)
      {
      FILE *ofile = fopen( "chinese.dat", "wb");

      printf( "Writing %d years,  year0 = %d\n", n_years, year0);
      fwrite( &n_years, sizeof( short), 1, ofile);
      fwrite( &year0, sizeof( short), 1, ofile);
      for( i = 0; i < n_years; i++)
         {
         unsigned long oval = (unsigned long)ydata[i].mask +
                (unsigned long)ydata[i].intercalary_month * 8192L +
                (unsigned long)ydata[i].offset * 8192L * 14L;
         fwrite( &oval, 3, 1, ofile);

         printf( "  Year %4d: offset %d, intercalary %2d, mask %04x\n",
            i + year0,
            ydata[i].offset,
            ydata[i].intercalary_month,
            ydata[i].mask);
         }
      fclose( ofile);
      }
   return( 0);
}
