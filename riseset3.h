#define PLANET_DATA struct planet_data

PLANET_DATA
   {
   double ecliptic_loc[3], equatorial_loc[3], altaz_loc[3];
   double r, ecliptic_lon, ecliptic_lat, jd;
   double hour_angle;
   };

int fill_planet_data( PLANET_DATA *pdata, const int planet_no, const double jd,
                  const double observer_lat, const double observer_lon,
                  const char *vsop_data);
double look_for_rise_set( const int planet_no,
                  const double jd0, const double jd1,
                  const double observer_lat, const double observer_lon,
                  const char *vsop_data, int *is_setting);
char *load_file_into_memory( const char *filename, size_t *filesize);
