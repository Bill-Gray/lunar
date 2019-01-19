#include <stdio.h>
#include <string.h>

/* See 'ades2mpc.cpp'.  This takes a list of ADES tags,  in the order
specified in the documentation,  and produces arrays suitable for use
in 'ades2mpc.cpp'.  It also sorts the tags;  this will allow for a
binary search,  gaining 0.001% (probably) in speed.  */

#define INTENTIONALLY_UNUSED_PARAMETER( param) (void)(param)

int main( const int intentionally_unused_argc,
          const char **intentionally_unused_argv)
{
   static const char *tags[] = {
         "permID", "provID", "artSat", "trkSub", "obsID", "trkID",
         "mode", "stn", "trx", "rcv", "sys", "ctr", "pos1", "pos2", "pos3",
         "posCov11", "posCov12", "posCov13",
         "posCov22", "posCov23", "posCov33",
         "prog", "obsTime", "ra", "dec", "raStar", "decStar", "obsCenter",
         "deltaRA", "deltaDec", "dist", "pa", "rmsRA", "rmsDec", "rmsDist",
         "rmsPA", "rmsCorr", "delay", "rmsDelay", "doppler", "rmsDoppler",
         "astCat", "mag", "rmsMag", "band", "photCat", "photAp",
         "nucMag", "logSNR", "seeing", "exp",
                     /* end p6, start p7 */
         "rmsFit", "nStars", "com",
         "frq", "ref", "disc", "subFmt", "subFrm",
                     /* end p7/start p8 */
         "precTime", "precRA", "precDec", "uncTime", "notes", "remarks",
         "deprecated", "localUse", "orbProd", "orbID", "resRA", "resDec",
         "selAst", "sigRA", "sigDec", "sigCorr", "sigTime", "biasRA",
         "biasDec", "biasTime", "photProd", "resMag", "selPhot",
                     /* end p9/start p10 */
         "sigMag", "biasMag", "photMod", "resDelay", "selDelay", "sigDelay",
         "resDoppler", "selDoppler", "sigDoppler", "observatory", "submitter",
         "observers", "measurers", "telescope", "software", "coinvestigators",
         "collaborators", "fundingSource",
                     /* end p10/start p11 */
         "comment", "optical", "offset", "occultation", "radar", "obsContext",
         "obsData", "obsBlock", "opticalResidual", "radarResidual", "ades",
                     /* start p 18: complex types */
         "MPCID", "OpticalID", "RadarID", "RadarValue", "Precision", "Location",
         "Photometry", "OffsetVal", "OpticalRes", "OpticalResMag",
         "OpticalResiduals", "RadarResiduals",
                     /* start p 25 */
         "mpcCode", "institution", "design", "aperture", "detector", "fRatio",
         "filter", "arrayScale", "pixelScale", "astrometry", "fitOrder",
         "photometry", "objectDetection", "line",
         "name",
                     /* added Sep 2018 */
         "rmsTime",
         NULL };
   size_t i, j;

   INTENTIONALLY_UNUSED_PARAMETER( intentionally_unused_argv);
   INTENTIONALLY_UNUSED_PARAMETER( intentionally_unused_argc);

   for( i = 0; tags[i]; i++)
      for( j = i + 1; tags[j]; j++)
         if( strcmp( tags[i], tags[j]) > 0)
            {
            const char *tptr = tags[j];

            tags[j] = tags[i];
            tags[i] = tptr;
            }

   j = 60;
   printf( "   static const char *tags[] = {");
   for( i = 0; tags[i]; i++)
      {
      int len = strlen( tags[i]) + 4;

      if( j + len > 60)
         {
         printf( "\n       ");
         j = 0;
         }
      printf( "\"%s\", ", tags[i]);
      j += len;
      }

   printf( "\n   NULL };\n\n");
   for( i = 0; tags[i]; i++)
      printf( "#define ADES_%-27s%4d\n", tags[i], (int)i + 1);
   return( 0);
}
