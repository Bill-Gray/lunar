double get_lunar_transit_time( const int year, const int month,
   const int day, const double latitude, const double longitude,
   const int time_zone, const int dst, const int real_transit);
void format_hh_mm( char *buff, const double time);
int get_zip_code_data( const int zip_code, double *latitude,
      double *longitude, int *time_zone, int *use_dst,
      char *place_name);
