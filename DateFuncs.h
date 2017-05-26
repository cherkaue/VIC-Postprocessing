// define leap year
#ifndef _LEAPYR
#define LEAPYR(y) (!((y)%400) || (!((y)%4) && ((y)%100)))
#endif

// define calendar variables
static char *SeasonNames[] = { "SPR", "SUM", "AUT", "WIN" };
static char *MonthNames[] = { "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };
static int   DaysInMonth[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
static int   SeasonNumbers[] = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2 };

// date structure
typedef struct {
  int day;
  int month;
  int year;
  int season;
  int doy;
  int rec;
  double hour;
  double juldate;
} DATE_STRUCT;

// date subroutines
DATE_STRUCT get_juldate( DATE_STRUCT self );
DATE_STRUCT get_season( DATE_STRUCT self );
DATE_STRUCT get_next_hour(DATE_STRUCT self);
DATE_STRUCT get_last_hour( DATE_STRUCT self );
DATE_STRUCT get_next_day( DATE_STRUCT self );
DATE_STRUCT get_last_day( DATE_STRUCT self );
DATE_STRUCT get_next_week( DATE_STRUCT self );
DATE_STRUCT get_last_week( DATE_STRUCT self );
DATE_STRUCT get_next_month( DATE_STRUCT self );
DATE_STRUCT get_last_month( DATE_STRUCT self );
DATE_STRUCT get_next_season( DATE_STRUCT self );
DATE_STRUCT get_last_season( DATE_STRUCT self );
DATE_STRUCT get_next_year( DATE_STRUCT self );
DATE_STRUCT get_last_year( DATE_STRUCT self );
DATE_STRUCT get_next_month_start( DATE_STRUCT self );
DATE_STRUCT get_next_season_start( DATE_STRUCT self );
double calc_juldate( int, int, int, double );
