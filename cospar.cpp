/* cospar.cpp: functions for planet/satellite orientations

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
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define FREE_INTERNAL_DATA       -99999
#define NO_OBJECT_SELECTED       -99999
#define PLANET_SELECTED          -99998

const char *cospar_filename = "cospar.txt";

/* Mostly complete code to compute COSPAR planetary,  satellite,  and
asteroid orientations using data extracted from 'cospar.txt'.  Thus far,
it can take the Guide-style object number and a JD and compute the pole
position and rotation angle Omega,  including the periodic and linear
and quadratic terms.  It still needs some of the logic from the original
'cospar.cpp' to convert this to a matrix and to handle odd cases such as
the Earth, which has a rotation matrix based on a separate precession
formula, and to set an identity matrix for unknown objects and such.
'cospar.txt' may be extended to include asteroids... someday. */

static int get_cospar_data_from_text_file( int object_number,
         const int system_number, const double jde,
         double *pole_ra, double *pole_dec, double *omega,
         bool *is_retrograde)
{
   const double J2000 = 2451545.0;        /* JD 2451545.0 = 1.5 Jan 2000 */
   const double d = (jde - J2000);
   const double t_cen = d / 36525.;
   static char **cospar_text = NULL;
   char buff[300];
   char planet = 0;
   int line, angular_coeffs_line = 0;
   int i, err = 0, done = 0, got_omega = 0;
   int curr_obj_from_file = NO_OBJECT_SELECTED;

   if( object_number == FREE_INTERNAL_DATA)
      {
      if( cospar_text)
         free( cospar_text);
      cospar_text = NULL;
      return( 0);
      }

   if( !cospar_text)       /* must load file */
      {
      FILE *ifile = fopen( cospar_filename, "rb");
      int pass;

      if( !ifile)
         return( -1);
                  /* make two passes through file:  one to count lines */
                  /* and file size, another to load the file           */
      for( pass = 0; pass < 2; pass++)
         {
         size_t bytes_read = 0;

         fseek( ifile, 0L, SEEK_SET);
         line = 0;
         while( fgets( buff, sizeof( buff), ifile) && memcmp( buff, "END", 3))
            {
            for( i = 0; buff[i] >= ' ' && buff[i] != '#'; i++)
               ;
            if( i)      /* yes,  it's a for-real line */
               {
               buff[i] = '\0';
                        /* remove redundant spaces: */
               if( memcmp( buff, "Remap:", 6))
                  for( i = 0; buff[i]; i++)
                     if( buff[i] == ' ')
                        {
                        memmove( buff + i, buff + i + 1, strlen( buff + i));
                        i--;
                        }
               if( pass)
                  {
                  strcpy( cospar_text[line], buff);
                  cospar_text[line + 1] = cospar_text[line] + i + 1;
                  }
               else  /* just counting bytes and lines */
                  bytes_read += (size_t)i + 1;
               line++;
               }
            }
         if( !pass)     /* we've counted lines & bytes; now alloc memory */
            {
            cospar_text = (char **)malloc( (size_t)(line + 1) * sizeof( char *)
                        + bytes_read);
            cospar_text[0] = (char *)(cospar_text + line + 1);
            }
         }
      cospar_text[line] = NULL;
      fclose( ifile);
      if( !pole_ra)        /* just loading coefficients;  not actually */
         return( 0);       /* computing orientations quite yet (see    */
      }                    /* load_cospar_file( ) below)               */
   *is_retrograde = false;
   for( line = 0; cospar_text[line] && !done && !err; line++)
      {
      char *tptr = cospar_text[line];

      if( *tptr == 'R')        /* "Remap:" */
         if( object_number > 9 && object_number < 1000) /* remap, e.g, 10 to 3001 */
             {                                          /* or 28 to 4001.  Only   */
             int loc = 6, bytes_read, idx1, idx2;   /* needed for idxs 10 to 999. */

             while( sscanf( tptr + loc, "%d %d%n", &idx1, &idx2, &bytes_read) == 2)
                {
                if( object_number == idx2)
                   object_number = idx1;
                loc += bytes_read;
                }
             }
      if( *tptr == 'P')   /* "Planet: "*/
         {
         planet = tptr[7];
         if( curr_obj_from_file == object_number)
            done = 1;
         curr_obj_from_file = PLANET_SELECTED;
         }
      else if( *tptr == 'O')    /* "Object:" */
         {
         if( curr_obj_from_file == object_number)
            done = 1;
         curr_obj_from_file = atoi( tptr + 4);
         }
      else if( planet && curr_obj_from_file == PLANET_SELECTED
               && tptr[0] == planet && tptr[1] == '1'
               && tptr[2] == '=')
         angular_coeffs_line = line - 1;
      else if( curr_obj_from_file == object_number)
         {
         double *oval = NULL;

         if( *tptr == 'a')    /* "a0=" */
            oval = pole_ra;
         else if( *tptr == 'd')   /* "d0=" */
            oval = pole_dec;
         else if( !got_omega && *tptr == 'W')
            if( tptr[1] == (char)(system_number + '0') || tptr[1] == '=')
               {
               got_omega = 1;
               oval = omega;
               }
         if( oval)
            {
            while( *tptr != '=')
               tptr++;
            tptr++;
            *oval = atof( tptr);
            if( *tptr == '-')     /* skip leading neg sign */
               tptr++;
            while( *tptr)
               if( *tptr != '+' && *tptr != '-')
                  tptr++;        /* just skip on over... */
               else
                  {
                  double coeff;
                  int number_length;

                  sscanf( tptr, "%lf%n", &coeff, &number_length);
                  tptr += number_length;
                  if( *tptr == 'd')
                     {
                     if( tptr[1] == '2')
                        coeff *= d;
                     else if( coeff < 0. && oval == omega)
                        *is_retrograde = true;
                     coeff *= d;
                     }
                  else if( *tptr == 'T')
                     {
                     if( tptr[1] == '2')
                        coeff *= t_cen;
                     else if( coeff < 0. && oval == omega)
                        *is_retrograde = true;
                     coeff *= t_cen;
                     }
                  else
                     {
                     int idx, multiplier = 1;
                     double angle;
                     double linear, constant_term;
                     char d_or_T, *ang_ptr;

                     if( tptr[3] == planet)
                        idx = atoi( tptr + 4);
                     else
                        {
                        multiplier = atoi( tptr + 3);
                        idx = atoi( tptr + 5);
                        if( tptr[4] != planet)
                           err = -4;
                        }
                     ang_ptr = cospar_text[angular_coeffs_line + idx] + 1;
                     while( ang_ptr[-1] != '=')
                        ang_ptr++;

                     sscanf( ang_ptr, "%lf%lf%c", &constant_term, &linear, &d_or_T);
                     angle = constant_term + linear *
                                   (d_or_T == 'd' ? d : t_cen);

                     angle *= (double)multiplier * PI / 180.;
                     if( !multiplier)
                        err = -5;
                     else if( angle == 0.)
                        err = -3;
                     else if( *tptr == 's')     /* sine term */
                        coeff *= sin( angle);
                     else if( *tptr == 'c')     /* cosine term */
                        coeff *= cos( angle);
                     else
                        err = -2;
                     }
                  *oval += coeff;
                  }
            }
         }
      }
   if( !err)
      if( !done)    /* never did find the object... fill with  */
         {          /* semi-random values and signal an error: */
         *pole_ra = *pole_dec = (double)( object_number * 20);
         if( omega)
            *omega = d * 360. / 1.3;   /* rotation once every 1.3 days */
         err = -1;
         }
#ifdef TEST_MAIN
   if( err && err != -1)
      printf( "ERROR %d: %s\n", err, buff);
#endif
   return( err);
}

/* Some programs (e.g.,  Find_Orb) put 'cospar.txt' in some directory
other than the working one.  They can use the following function to
unload the COSPAR data (if any) and then direct that it be loaded from
a specified file.          */

int DLL_FUNC load_cospar_file( const char *filename)
{
   const char *temp_name = cospar_filename;
   int rval, pass;

   cospar_filename = filename;
   for( pass = 0; pass < (filename ? 2 : 1); pass++)
      rval = get_cospar_data_from_text_file( (pass ? 0 : FREE_INTERNAL_DATA),
                    0, 0., NULL, NULL, NULL, NULL);
   cospar_filename = temp_name;
   return( rval);
}

/* The returned matrix contains three J2000 equatorial unit vectors :

   matrix[0, 1, 2] = vector in direction of lat=0, lon=0;
   matrix[3, 4, 5] = vector in direction of lat=0, lon=90;
   matrix[6, 7, 8] = vector in direction of planet's north pole (lat=90)

   For a prograde object,  this will form a right-handed coordinate system.
For a retrograde object,  the IAU's boneheaded convention results in a
left-handed system (see the 'if( is_retrograde)' code that flips the
middle of the above three vectors).  */

#ifdef TEST_MAIN
int DLL_FUNC calc_planet_orientation2( const int planet_no, const int system_no,
                  const double jd, double *matrix)
#else
int DLL_FUNC calc_planet_orientation( const int planet_no, const int system_no,
                  const double jd, double *matrix)
#endif
{
   static int prev_planet_no = NO_OBJECT_SELECTED;
   static int prev_system_no = NO_OBJECT_SELECTED, prev_rval = 0;
   static double prev_jd = -1.;
   static double prev_matrix[9];
   int i, rval;
   bool is_retrograde;
   double pole_ra, pole_dec, omega;

   if( planet_no == prev_planet_no && system_no == prev_system_no
                           && jd == prev_jd)
      {
      memcpy( matrix, prev_matrix, 9 * sizeof( double));
      return( prev_rval);
      }

   prev_planet_no = planet_no;
   prev_system_no = system_no;
   prev_jd = jd;

   if( planet_no == 3)        /* handle earth with "normal" precession: */
      {
      const double J2000 = 2451545.;   /* 1.5 Jan 2000 = JD 2451545 */
      const double t_cen = (jd - J2000) / 36525.;

//    setup_precession( matrix, 2000., 2000. + t_cen * 100.);
      setup_precession_with_nutation( matrix, 2000. + t_cen * 100.);
      for( i = 3; i < 6; i++)
         matrix[i] = -matrix[i];
      spin_matrix( matrix, matrix + 3, green_sidereal_time( jd));
      memcpy( prev_matrix, matrix, 9 * sizeof( double));
      prev_rval = 0;
      return( 0);
      }

         /* For everybody else,  we use TD.  Only the earth uses UT. */
         /* (This correction added 5 Nov 98,  after G Seronik pointed */
         /* out an error in the Saturn central meridian.)             */

   rval = get_cospar_data_from_text_file( planet_no, system_no,
              jd + td_minus_ut( jd) / seconds_per_day,
              &pole_ra, &pole_dec, &omega, &is_retrograde);
   pole_ra *= PI / 180.;
   pole_dec *= PI / 180.;
   polar3_to_cartesian( matrix, pole_ra - PI / 2., 0.);
   polar3_to_cartesian( matrix + 3, pole_ra - PI, PI / 2. - pole_dec);
   polar3_to_cartesian( matrix + 6, pole_ra, pole_dec);

   spin_matrix( matrix, matrix + 3, omega * PI / 180. + PI);
   if( is_retrograde)
      for( i = 3; i < 6; i++)
         matrix[i] *= -1.;
   memcpy( prev_matrix, matrix, 9 * sizeof( double));
   prev_rval = rval;
   return( rval);
}


#ifdef TEST_MAIN
void main( int argc, char **argv)
{
   double pole_ra, pole_dec, omega;
   const int planet_number = atoi( argv[1]);
   const double jde = atof( argv[2]);
   const int system_number = (argc > 3 ? atoi( argv[3]) : 0);
   int i, j;

   for( i = (planet_number == -1 ? 0 : planet_number);
        i < (planet_number == -1 ? 100 : planet_number + 1); i++)
      {
      int err = get_cospar_data_from_text_file( i, system_number, jde,
                      &pole_ra, &pole_dec, &omega);

      printf( "Planet %d\n", i);
      if( !err)
         {
         double old_mat[9], new_mat[9], delta = 0.;

         printf( "   pole RA: %lf\n", pole_ra);
         printf( "   pole dec %lf\n", pole_dec);
         printf( "   Omega    %lf (%lf)\n", omega, fmod( omega, 360.));
         if( !calc_planet_orientation2( i, system_number, jde, new_mat))
            if( !calc_planet_orientation( i, system_number, jde, old_mat))
               {
               for( j = 0; j < 9; j += 3)
                  printf( "%10.6lf %10.6lf %10.6lf    %10.6lf %10.6lf %10.6lf\n",
                           new_mat[j], new_mat[j + 1], new_mat[j + 2],
                           old_mat[j], old_mat[j + 1], old_mat[j + 2]);
               for( j = 0; j < 9; j++)
                  delta += (new_mat[j] - old_mat[j]) * (new_mat[j] - old_mat[j]);
               printf( "Diff: %lf\n", sqrt( delta));
               }
         }
      else
         printf( "   Error %d\n", err);
      }
}
#endif
