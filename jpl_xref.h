/* The following table is going to need occasional fixes,  mostly
as new spacecraft are launched.  Cas = Cassini,  SoO = Solar Orbiter,
etc. are _not_ official MPC codes;  if they ever get them,  that
will also cause updates to be made.  The list does not reflect
every spacecraft listed on Horizons;  it lists those that have
either gotten astrometry or for which I compute TLEs.  */

typedef struct
{
   const char mpc_code[4];
   int jpl_desig, norad_number;
   const char *intl_desig, *name;
} jpl_xref_t;

/* MPC code    JPL#   NORAD#  Intl desig    Common name */
static const jpl_xref_t jpl_xrefs[] = {
   { "249",      -21, 23726, "1995-065A",   "SOHO" },
   { "Jui",      -28, 56176, "2023-053A",   "JUICE" },
   { "Ha2",      -37, 40319, "2014-076A",   "Hayabusa 2" },
   { "250",      -48, 20580, "1990-037B",   "Hubble Space Telescope" },
   { "Luc",      -49, 49328, "2021-093A",   "Lucy" },
   { "OsR",      -64, 41757, "2016-055A",   "OSIRIS-REx" },
   { "Cas",      -82, 25008, "1997-061A",   "Cassini" },
   { "245",      -79, 27871, "2003-038A",   "Spitzer Space Telescope" },
   { "C57",      -95, 43435, "2018-038A",   "TESS" },
   { "PSP",      -96, 43592, "2018-065A",   "Parker Solar Probe" },
   { "C54",      -98, 28928, "2006-001A",   "New Horizons" },
   { "Equ",     -101, 79970, "2022-156ZZZ", "EQUULEUS" },
   { "SoO",     -144, 45167, "2020-010A",   "Solar Orbiter" },
   { "Cha",     -151, 25867, "1999-040B",   "Chandra X-ray Observatory" },
   { "Cdr",     -158, 57320, "2023-098A",   "Chandrayaan-3" },
   { "C51",     -163, 36119, "2009-071A",   "WISE" },
   { "LFL",     -164, 54697, "2022-168B",   "Lunar Flashlight" },
   { "274",     -170, 50463, "2021-130A",   "James Webb Space Telescope" },
   { "C55",     -227, 34380, "2009-011A",   "Kepler" },
   { "C49",     -234, 29510, "2006-047A",   "STEREO-A" },
   { "C50",     -235, 29511, "2006-047B",   "STEREO-B" },
   { "Sli",     -240, 57803, "2023-137D",   "SLIM" },
   { "Psy",     -255, 58049, "2023-157A",   "Psyche" },
   { "273",     -680, 57209, "2023-092A",   "Euclid" },
   { "C52",  -128485, 28485, "2004-047A",   "Swift" },
   { "C53",  -139089, 39089, "2013-009D",   "NEOSSat" },
   { "258",  -139479, 39479, "2013-074A",   "Gaia" },
   { "C56",  -141043, 41043, "2015-070A",   "LISA Pathfinder" },
   { "C59",  -148840, 48841, "2021-050B",   "Yangwang 1" },
   { "Tia", -9901491, 45935, "2020-049A",   "Tianwen-1" } };
/* MPC code    JPL#   NORAD#  Intl desig    Common name */
