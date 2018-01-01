void avoid_runaway_process( const int max_time_to_run);   /* cgi_func.c */
int get_urlencoded_form_data( const char **idata,       /* cgi_func.c */
                              char *field, const size_t max_field,
                              char *buff, const size_t max_buff);
int get_multipart_form_data( const char *boundary, char *field,
                char *buff, char *filename, const size_t max_len);
