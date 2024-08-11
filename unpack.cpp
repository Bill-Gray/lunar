#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "watdefs.h"
#include "mpc_func.h"
#include "stringex.h"

/* "Mutant hex" is frequently used by MPC.  It uses the usual hex digits
   0123456789ABCDEF for numbers 0 to 15,  followed by G...Z for 16...35
   and a...z for 36...61,  to encode numbers in base 62.  */

int mutant_hex_char_to_int( const char c)
{
   int rval = -1;

   if( c >= 'a')
      {
      if( c <= 'z')
         rval = (int)c - 'a' + 36;
      }
   else if( c >= 'A')
      {
      if( c <= 'Z')
         rval = (int)c - 'A' + 10;
      }
   else if( c >= '0' && c <= '9')
      rval = (int)c - '0';
   return( rval);
}

char int_to_mutant_hex_char( const int ival)
{
   int rval;

   if( ival < 0 || ival > 61)
      rval = '\0';
   else if( ival < 10)
      rval = '0';
   else if( ival < 36)
      rval = 'A' - 10;
   else
      rval = 'a' - 36;
   return( rval ? (char)( rval + ival) : '\0');
}

int get_mutant_hex_value( const char *buff, size_t n_digits)
{
   int rval = 0;

   while( n_digits--)
      {
      const int digit = mutant_hex_char_to_int( *buff++);

      if( digit == -1)
         return( -1);
      rval = rval * 62 + digit;
      }
   return( rval);
}

/* If this returns a non-zero result,  'value' was too large to fit into
the allocated space (i.e.,  value >= 62 ^ n_digits).                  */

int encode_value_in_mutant_hex( char *buff, size_t n_digits, int value)
{
   buff += n_digits - 1;
   while( n_digits--)
      {
      *buff-- = int_to_mutant_hex_char( value % 62);
      value /= 62;
      }
   return( value);
}

/* This will unpack a packed designation such as 'K04J42X' into
'2004 JX42'. Returns 0 if it's a packed desig,  non-zero otherwise. Call
with obuff == NULL just to find out if ibuff is actually a packed desig. */

static int unpack_provisional_packed_desig( char *obuff, const size_t obuff_size,
                                        const char *ibuff)
{
   int rval = 0;

   if( *ibuff ==  '_' && ibuff[2] >= 'A' && ibuff[2] <= 'Z')
      {
      const int year = mutant_hex_char_to_int( ibuff[1]);
      const int desig_num = get_mutant_hex_value( ibuff + 3, 4);

      if( year < 0 || desig_num < 0)
         rval = -1;
      else
         {
         char letter = 'A' + desig_num % 25;

         if( letter >= 'I')
            letter++;
         snprintf_err( obuff, 14, "%d %c%c%d", year + 2000, ibuff[2],
                           letter, 620 + desig_num / 25);
         }
      }
   else if( *ibuff >= 'G' && *ibuff <= 'K' && isdigit( ibuff[1])
            && isdigit( ibuff[2]) && isupper( ibuff[3])
            && isdigit( ibuff[5]) && (isalpha( ibuff[6]) || ibuff[6] == '0'))
      {
      int output_no = mutant_hex_char_to_int( ibuff[4]);

      if( output_no == -1)
         rval = -2;

      if( obuff)
         {
         obuff[0] = ((*ibuff >= 'K') ? '2' : '1');          /* millennium */
         obuff[1] = (char)( '0' + ((*ibuff - 'A') % 10));   /* century   */
         obuff[2] = ibuff[1];                               /* decade   */
         obuff[3] = ibuff[2];                               /* year    */
         obuff[4] = ' ';
         obuff[5] = ibuff[3];       /* half-month designator */
         if( isupper( ibuff[6]))    /* asteroid second letter */
            {
            obuff[6] = ibuff[6];
            obuff[7] = '\0';
            }
         else
            obuff[6] = '\0';
         output_no = output_no * 10 + ibuff[5] - '0';
         if( output_no)
            snprintf_append( obuff, obuff_size, "%d", output_no);
         if( islower( ibuff[6]))    /* comet fragment letter */
            snprintf_append( obuff, obuff_size, "-%c", ibuff[6] + 'A' - 'a');
         }
      }
   else
      {
      if( obuff)
         *obuff = '\0';
      rval = -1;
      }
   return( rval);
}

/* Looks for designations of the form YYYY-NNN(letter).  */

static int is_artsat_desig( const char *desig)
{
   size_t slen;

   for( slen = strlen( desig); slen > 8; desig++, slen--)
      if( desig[4] == '-' && atoi( desig) > 1956 && atoi( desig) < 2100
               && isdigit( desig[5]) && isdigit( desig[6])
               && isdigit( desig[7]) && isupper( desig[8]))
         return( 1);
   return( 0);
}

const char *planet_names_in_english[] = { "Venus", "Earth", "Mars", "Jupiter",
                  "Saturn", "Uranus", "Neptune", "Pluto" };

/* In an MPC astrometric report line,  the name can be stored in assorted
   highly scrambled and non-intuitive ways.  Those ways _are_ documented
   on the MPC Web site,  which is probably the best place to look to
   understand why the following bizarre code does what it does:

   https://www.minorplanetcenter.org/iau/info/PackedDes.html      */

int unpack_mpc_desig( char *obuff, const char *packed)
{
   int rval = OBJ_DESIG_OTHER;
   size_t i;
   int mask = 1, digit_mask = 0, space_mask = 0, digit;
   char provisional_desig[40];
   const size_t obuff_size = 16;

   if( *packed == '$')         /* Find_Orb extension to allow storing of */
      {                        /* an unpacked name,  up to 11 chars */
      if( obuff)
         strlcpy_err( obuff, packed + 1, obuff_size);
      return( OBJ_DESIG_ASTEROID_NUMBERED);     /* fix later... will be a project! */
      }
   for( i = 0; i < 12; i++, mask <<= 1)
      if( isdigit( packed[i]))
         digit_mask |= mask;
      else if( packed[i] == ' ')
         space_mask |= mask;

   if( packed[4] == 'S')   /* Possible natural satellite */
      {
      if( strchr( "MVEJSUNP", *packed) && (digit_mask & 0xe) == 0xe
                  && (space_mask & 0xfe0) == 0xfe0)
         {
         if( obuff)
            {
            const char *roman_digits[10] = { "", "I", "II", "III", "IV",
                     "V", "VI", "VII", "VIII", "IX" };
            const char *roman_tens[10] = { "", "X", "XX", "XXX", "XL",
                     "L", "LX", "LXX", "LXXX", "XC" };
            const char *roman_hundreds[10] = { "", "C", "CC", "CCC", "CD",
                     "D", "DC", "DCC", "DCCC", "CM" };
            const int obj_number = atoi( packed + 1);

            for( i = 0; i < 8; i++)
               if( planet_names_in_english[i][0] == *packed)
                  {
                  strlcpy( obuff, planet_names_in_english[i], obuff_size);
                  strlcat( obuff, " ", obuff_size);
                  }
            if( obj_number / 100)
               strlcat_err( obuff, roman_hundreds[obj_number / 100], obuff_size);
            if( (obj_number / 10) % 10)
               strlcat_err( obuff, roman_tens[(obj_number / 10) % 10], obuff_size);
            if( obj_number % 10)
               strlcat_err( obuff, roman_digits[obj_number % 10], obuff_size);
            }
         rval = OBJ_DESIG_NATSAT_NUMBERED;
         }
      else if( strchr( "MVEJSUNP", packed[8])
               && (digit_mask & 0xcc0) == 0xcc0 && space_mask == 0xf
               && (digit = mutant_hex_char_to_int( packed[9])) >= 0
               && packed[11] == '0' && packed[5] >= 'H' && packed[5] <= 'Z')
         {
         if( obuff)
            {
            snprintf_err( obuff, 7, "S/%d%c%c", 20 + packed[5] - 'K',
                                     packed[6], packed[7]);
            obuff[6] = ' ';
            obuff[7] = packed[8];       /* planet identifier */
            obuff[8] = ' ';
            i = 9;
            if( digit > 10)
               obuff[i++] = '0' + digit / 10;
            if( digit)
               obuff[i++] = '0' + digit % 10;
            obuff[i++] = packed[10];
            obuff[i] = '\0';
            }
         rval = OBJ_DESIG_NATSAT_PROVISIONAL;
         }
      }
   unpack_provisional_packed_desig( provisional_desig, sizeof( provisional_desig),
                                 packed + 5);

            /* Check for numbered asteroids (620000) or greater.  See MPEC
            2019-O55 for explanation of this 'extended numbering scheme'. */
   if( packed[0] == '~' && ((space_mask & 0x3e0) == 0x3e0 || !(space_mask & 0x3e0)))
      {
      const int num = get_mutant_hex_value( packed + 1, 4);

      if( num >= 0)    /* yes,  it's a valid 'extended' numbered object */
         {
         if( obuff)
            {
            snprintf_err( obuff, 11, "(%d)", num + 620000);
            if( *provisional_desig)
               snprintf_append( obuff, 24, " %s", provisional_desig);
            }
         return( OBJ_DESIG_ASTEROID_NUMBERED);
         }
      }
            /* For numbered asteroids or comets,  we require either that
            columns 6-12 be blank (usually the case) or have a valid
            provisional designation.  And,  of course,  the packed desig must
            start with an alphanumeric and be followed by three digits.
            Four,  if it's a numbered asteroid. */
   if( isalnum( packed[0]) && (digit_mask & 0xe) == 0xe
                  && ((space_mask & 0x3e0) == 0x3e0 || *provisional_desig))
      {
      if( isdigit( packed[4]))
         {                                /* it's a numbered asteroid */
         int number = mutant_hex_char_to_int( *packed);

         if( number >= 0)
            {
            number = number * 10000L + atol( packed + 1);
            rval = OBJ_DESIG_ASTEROID_NUMBERED;
            if( obuff)
               {
               snprintf_err( obuff, 13, "(%d)", number);
                  /* Desig may be,  e.g., "U4330K06SL7X" : both the number */
                  /* and the provisional desig.  We'd like to decipher this */
                  /* as (for the example) "304330 = 2006 SX217".            */
               if( *provisional_desig)
                  snprintf_append( obuff, 23, " = %s", provisional_desig);
               }
            }
         }
      else if( strchr( "PCDXAI", packed[4]))   /* it's a numbered comet */
         {
         if( obuff)
            {
            int j = 0;
            const char suffix_char = packed[11];

            i = 0;
            while( packed[i] == '0')         /* skip leading zeroes */
               i++;
            while( packed[i] >= '0' && packed[i] <= '9')   /* read digits... */
               obuff[j++] = packed[i++];
            obuff[j++] = packed[4];
            if( suffix_char >= 'a' && suffix_char <= 'z')
               {
               const char extra_suffix_char = packed[10];

               obuff[j++] = '-';
                  /* possibly _two_ suffix letters... so far,  only the  */
                  /* components of 73P/Schwassmann-Wachmann 3 have this: */
               if( extra_suffix_char >= 'a' && extra_suffix_char <= 'z')
                  obuff[j++] = extra_suffix_char + 'A' - 'a';
               obuff[j++] = suffix_char + 'A' - 'a';
               }
            obuff[j] = '\0';
            }
         rval = OBJ_DESIG_COMET_NUMBERED;
         }
      }

   if( rval == OBJ_DESIG_OTHER && !memcmp( packed, "    ", 4)
            && strchr( " PCDXA", packed[4]))
      {
      if( obuff && packed[4] != ' ')
         {
         *obuff++ = packed[4];
         *obuff++ = '/';
         }
      if( !unpack_provisional_packed_desig( obuff, obuff_size, packed + 5))
         rval = ((packed[4] == ' ' || packed[4] == 'A') ?
                                      OBJ_DESIG_ASTEROID_PROVISIONAL
                                    : OBJ_DESIG_COMET_PROVISIONAL);
      }
   if( rval == OBJ_DESIG_OTHER && space_mask == 0x1f
            && (digit_mask & 0xf00) == 0xf00 && packed[7] == 'S')
      {
      if( (packed[5] == 'P' && packed[6] == 'L') ||
          (packed[5] == 'T' && packed[6] >= '1' && packed[6] <= '3'))
         {
         if( obuff)
            {
            memcpy( obuff, packed + 8, 4);
            obuff[4] = ' ';
            obuff[5] = packed[5];
            obuff[6] = '-';
            obuff[7] = packed[6];
            obuff[8] = '\0';
            }
         rval = OBJ_DESIG_ASTEROID_PROVISIONAL;
         }
      }

   if( rval == OBJ_DESIG_OTHER && obuff)    /* store the name "as is",   */
      {             /* assuming no encoding (except skip leading spaces) */
      for( i = 0; i < 12 && packed[i] == ' '; i++)
         ;
      memcpy( obuff, packed + i, 12 - i);
      i = 12 - i;
      while( i && obuff[i - 1] <= ' ')
         i--;
      obuff[i] = '\0';
      }

   if( rval == OBJ_DESIG_OTHER && is_artsat_desig( packed))
      rval = OBJ_DESIG_ARTSAT;
   return( rval);
}

/* The above function will work if the designation in question is
aligned in the correct columns.  The following function will generally
unpack designations even if they're misaligned. */

int unpack_unaligned_mpc_desig( char *obuff, const char *packed)
{
   size_t len = 0, column = 0;
   char aligned[13];

   while( *packed == ' ')
      packed++;
   while( packed[len] > ' ' && packed[len] < 127)
      len++;
   if( len > 12)
      {
      memcpy( obuff, packed, 12);
      return( OBJ_DESIG_OTHER);
      }
   memset( aligned, ' ', 12);
   aligned[12] = '\0';
   if( len == 7 || len == 8)
      column = (size_t)( 12 - len);
   memcpy( aligned + column, packed, len);
   return( unpack_mpc_desig( obuff, aligned));
}
