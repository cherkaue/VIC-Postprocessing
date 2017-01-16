// Prototype stats_1d functions
extern double   get_mean(double *,int,double);
extern double   get_var(double *,int,double,double);
extern double   get_stdev(double *,int,double,double);
extern double   get_skew(double *,int,double,double,double);
extern double   get_sum(double *,int,double);
extern double   get_min(double *,int,double);
extern double   get_max(double *,int,double);
extern int      get_count_over_thres(double *, double, int, double);
extern int      get_count_under_thres(double *, double, int, double);
extern int      get_first_day_over_thres(double *, double, int, double);
extern int      get_first_day_under_thres(double *, double, int, double);
extern int      get_last_day_over_thres(double *, double, int, double);
extern int      get_last_day_under_thres(double *, double, int, double);
extern void     get_3LN(double,double,double,double *,double *,double *);
extern double   get_rk(int,int,double *,double *,double,double);
