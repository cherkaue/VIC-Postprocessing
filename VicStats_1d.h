#ifndef DATEFUNCS
#include <DateFuncs.h>
#endif

#ifndef VICSTATS
#define VICSTATS

// Prototype stats_1d functions
extern double   get_mean(double *,int,double);
extern double   get_var(double *,int,double,double);
extern double   get_stdev(double *,int,double,double);
extern double   get_skew(double *,int,double,double,double);
extern double   get_sum(double *,int,double);
extern double   get_min(double *,int,double);
extern double   get_max(double *,int,double);
extern int      get_count_over_thres(double *, double, int, double);
extern double   get_average_days_over_thres(double *, double, int, double);
extern double   get_average_value_over_thres(double *, double, int, double);
extern int      get_consecutive_days_over_thres(double *, double, int, double);
extern int      get_count_under_thres(double *, double, int, double);
extern double   get_average_days_under_thres(double *, double, int, double);
extern double   get_average_value_under_thres(double *, double, int, double);
extern int      get_consecutive_days_under_thres(double *, double, int, double);
extern int      get_first_day_over_thres(double *, double, int, double, DATE_STRUCT);
extern int      get_first_day_over_thres_from_middle(double *, double, int, double, DATE_STRUCT);
extern int      get_first_day_under_thres(double *, double, int, double, DATE_STRUCT);
extern int      get_first_day_under_thres_from_middle(double *, double, int, double, DATE_STRUCT);
extern int      get_last_day_over_thres(double *, double, int, double, DATE_STRUCT);
extern int      get_last_day_over_thres_from_middle(double *, double, int, double, DATE_STRUCT);
extern int      get_last_day_under_thres(double *, double, int, double, DATE_STRUCT);
extern int      get_last_day_under_thres_from_middle(double *, double, int, double, DATE_STRUCT);
extern int      get_count_times_thres_crossed(double *, double, int, double);
extern int      get_modified_growing_degree_day(double *, double *, int, double);
extern double   get_RB_index(double *, int, double);
extern double   get_TQ_mean(double *, int, double);
extern double   get_seven_day_low_value(double *, int, double);
extern double   get_quantile(double *, double, int, double);
extern void     get_3LN(double,double,double,double *,double *,double *);
extern double   get_rk(int,int,double *,double *,double,double);
extern void     swap(double *, double *);
extern int      partition (double *, int, int);
extern void     quickSort(double *, int, int);

#endif
