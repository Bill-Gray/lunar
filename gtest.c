#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Code created largely out of curiosity about the effects of
roundoff errors.  Greenwich sidereal time is computed as a linear
function of UT1 time since 2000 Jan 1 = JD 2451545.0 :

(1) gst = omega0 + rate * (jd - 2451545.0) + quadratic and cubic terms

   The official IAU 1976 expression in seconds of time is :

GMST= 67310.548 + (3155760000. + 8640184.812866)*T_cen
    = 67310.548 + (86400. + 8640184.812866 / 36525) * (jd - J2000)
                     + quadratic and cubic terms

   (the quad/cubic terms are used below,  but they're irrelevant
for this bit of discussion.)

rate = 360 + excess

   where excess is 0.98564... degrees,  or 236.5553... seconds.  That
is to say,  each UT day,  the earth turns 360 degrees plus (almost)
one degree.  That 'excess' accumulates to 360 degrees over the course
of a year, so it's roughly 360/365.2475 degrees.

   We're interested in 'gst' in degrees,  and therefore not interested
in multiples of 360 degrees.  If we used the above equation,  we might
follow it with

(2) gst = fmod( gst, 360.);

   As we get further from jd = 2451545.0,  we'll see roundoff error
as a result.  (1) will produce a large number,  and (2) will effectively
subtract a nearly equal large number.  In _Astronomical Algorithms_,
Meeus suggests avoiding this with

(3a) ijd = (int)jd;    (select nearest,  or at least nearby,  integer JD)
(3b) gst = omega0 + rate * (jd - ijd) + excess * (ijd - 2451545.0) + terms

   and then,  as before,  applying (2) to discard unneeded full rotations.
Mathematically,  (3b) produces the same result as (1),  but with fewer
excess rotations.  (The value is roughly 365.25 times smaller as jd
tends to +/- inf.)  So the roundoff error is drastically decreased.

   (I should perhaps mention that Meeus did a lot of work in an era
where floating-point math was often done with 32-bit quantities,
giving about six digits of precision.  Which made tricks of this sort
especially important.)

   But for large/small enough values of jd,  excess * (jd - 2451545.0)
can still be a large number.  So I decided to try to make use of the
fact that in 1461 days (four years of 365 days plus a leap day),
the earth makes about 1465 rotations,  plus an 'excess_1461' angle
of a mere 0.030802 degrees :

(4) gst = omega0 + rate * (jd - ijd)
         + excess *      ((ijd - 2451545) % 1461)
         + excess_1461 * ((ijd - 2451545) / 1461) + terms

   (in which the divisions and modulus 1461 are done in integer math.)
This results in a still greater decrease in the magnitude of gst,
causing still less roundoff after (2) is applied.

   The whole point is somewhat theoretical here.  Our knowledge of the
earth's rotation,  and its unpredictable irregular nature,  are greater
sources of error.  But it _would_ be relevant in similar situations
where lots of rotations accumulate and one then wants to know the
remaining angle.  Or any linear function where your interest is in a
fractional remainder. */

const double pi =
   3.1415926535897932384626433832795028841971693993751058209749445923;

const double omega0 = 280.460618375;
const double excess =
   0.9856473662863335614875655943417750399269906456764772986;
const double quad_term =
   0.0003879333333333333333333333333333333333333333333333333333;
const double cubic_term =
   -2.58333333333e-8;

double green_sidereal_time( double jd_ut)
{
   double t_cen, rval;

   jd_ut -= 2451545.0;        /* set relative to 2000.0 */
   t_cen = jd_ut / 36525.;    /* convert to julian centuries */
   rval = omega0 + (360 + excess) * jd_ut
            + t_cen * t_cen * (quad_term + t_cen * cubic_term);

         /* See p 84,  in Meeus:  the following should get apparent */
         /* Greenwich sidereal time: */
   return( rval * pi / 180.);
}

double green_sidereal_time_better( double jd_ut)
{
   double t_cen, rval, base_t;

   jd_ut -= 2451545.0;        /* set relative to 2000.0 */
   t_cen = jd_ut / 36525.;    /* convert to julian centuries */
   base_t = floor( jd_ut);
   jd_ut -= base_t;
   rval = omega0 + (360. + excess) * jd_ut + excess * base_t
            + t_cen * t_cen * (quad_term + t_cen * cubic_term);

         /* See p 84,  in Meeus:  the following should get apparent */
         /* Greenwich sidereal time: */
   return( rval * pi / 180.);
}

double green_sidereal_time_best( double jd_ut)
{
   double t_cen, rval, base_t;
   long idays;
   const double excess_1461 =
      0.030802144333333333333333333333333333333333333333333333333333333333;

   jd_ut -= 2451545.0;        /* set relative to 2000.0 */
   t_cen = jd_ut / 36525.;    /* convert to julian centuries */
   base_t = floor( jd_ut);
   jd_ut -= base_t;
   idays = (long)base_t;
   printf( "idays %ld jd_ut %f base_t %f\n", idays, jd_ut, base_t);
   rval = omega0 + (360. + excess) * jd_ut
            + (double)( idays / 1461L) * excess_1461
            + (double)( idays % 1461L) * excess
            + t_cen * t_cen * (quad_term + t_cen * cubic_term);

         /* See p 84,  in Meeus:  the following should get apparent */
         /* Greenwich sidereal time: */
   return( rval * pi / 180.);
}

static double fix_ang( double ang)
{
   ang = fmod( ang, 360.);
   if( ang < 0.)
      ang += 360.;
   return( ang);
}

int main( const int argc, const char **argv)
{
   const double jd = atof( argv[1]);
   double v1 = green_sidereal_time( jd) * 180. / pi;
   double v2 = green_sidereal_time_better( jd) * 180. / pi;
   double v3 = green_sidereal_time_best( jd) * 180. / pi;

   printf( "v1 %30.15f %25.15f (simple)\n", v1, fix_ang( v1));
   printf( "v2 %30.15f %25.15f (Meeus) \n", v2, fix_ang( v2));
   printf( "v3 %30.15f %25.15f ('best')\n", v3, fix_ang( v3));
   return( 0);
}
