#ifndef MPC_FUNC_H_INCLUDED
#define MPC_FUNC_H_INCLUDED

bool is_valid_mpc_code( const char *mpc_code);        /* mpc_fmt.cpp */
double extract_date_from_mpc_report( const char *buff, unsigned *format);
int get_ra_dec_from_mpc_report( const char *ibuff,    /* mpc_fmt.cpp */
                       int *ra_format, double *ra, double *ra_precision,
                       int *dec_format, double *dec, double *dec_precision);

char net_name_to_byte_code( const char *net_name);
const char *byte_code_to_net_name( const char byte_code);

typedef struct
{
   double lat, lon;     /* in radians */
   double alt;          /* in metres */
   double rho_cos_phi, rho_sin_phi;    /* in planet radii */
   int planet, prec1, prec2, format;
   const char *name;
   char code[5];
} mpc_code_t;

#define MPC_CODE_PARALLAXES         1
#define MPC_CODE_LAT_LON_ALT        2
#define MPC_CODE_SATELLITE          3

int get_mpc_code_info( mpc_code_t *cinfo, const char *buff);
double point_to_ellipse( const double a, const double b,
                         const double x, const double y, double *dist);

#endif
