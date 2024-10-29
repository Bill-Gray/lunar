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
   { "Cap",    -1176, 52914, "2022-070A",   "CAPSTONE" },
   { "Cas",      -82, 25008, "1997-061A",   "Cassini" },
   { "Cha",     -151, 25867, "1999-040B",   "Chandra X-ray Observatory" },
   { "Cdr",     -158, 57320, "2023-098A",   "Chandrayaan-3" },
   { "Equ",     -101, 79970, "2022-156ZZZ", "EQUULEUS" },
   { "273",     -680, 57209, "2023-092A",   "Euclid" },
   { "EuC",     -159, 61507, "2024-182A",   "Europa Clipper" },
   { "258",  -139479, 39479, "2013-074A",   "Gaia" },
   { "Goe",  -160133, 60133, "2024-119A",   "GOES-19" },
   { "Ha2",      -37, 40319, "2014-076A",   "Hayabusa 2" },
   { "Her",      -91, 61449, "2024-180A",   "Hera", },
   { "250",      -48, 20580, "1990-037B",   "Hubble Space Telescope" },
   { "IM1",     -229, 58963, "2024-030A",   "IM-1 (Odysseus)" },
   { "274",     -170, 50463, "2021-130A",   "James Webb Space Telescope" },
   { "Jui",      -28, 56176, "2023-053A",   "JUICE" },
   { "C55",     -227, 34380, "2009-011A",   "Kepler" },
   { "C56",  -141043, 41043, "2015-070A",   "LISA Pathfinder" },
   { "Luc",      -49, 49328, "2021-093A",   "Lucy" },
   { "LFL",     -164, 54697, "2022-168B",   "Lunar Flashlight" },
   { "C53",  -139089, 39089, "2013-009D",   "NEOSSat" },
   { "C58",      -33, 99999, "2025-000ZZZ", "NEO Surveyor" },
   { "C54",      -98, 28928, "2006-001A",   "New Horizons" },
   { "OsR",      -64, 41757, "2016-055A",   "OSIRIS-REx" },
   { "PSP",      -96, 43592, "2018-065A",   "Parker Solar Probe" },
   { "Prg",     -244, 58751, "2024-006A",   "Peregrine" },
   { "Psy",     -255, 58049, "2023-157A",   "Psyche" },
   { "Sli",     -240, 57803, "2023-137D",   "SLIM" },
   { "249",      -21, 23726, "1995-065A",   "SOHO" },
   { "SoO",     -144, 45167, "2020-010A",   "Solar Orbiter" },
   { "245",      -79, 27871, "2003-038A",   "Spitzer Space Telescope" },
   { "C49",     -234, 29510, "2006-047A",   "STEREO-A" },
   { "C50",     -235, 29511, "2006-047B",   "STEREO-B" },
   { "C52",  -128485, 28485, "2004-047A",   "Swift" },
   { "C57",      -95, 43435, "2018-038A",   "TESS" },
   { "Tia", -9901491, 45935, "2020-049A",   "Tianwen-1" },
   { "C51",     -163, 36119, "2009-071A",   "WISE" },
   { "C59",  -148840, 48841, "2021-050B",   "Yangwang 1" } };
/* MPC code    JPL#   NORAD#  Intl desig    Common name */
