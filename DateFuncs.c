#include <VicStuff.h>

double calc_juldate(int year,
		    int month,
		    int day,
		    double time)
/***********************************************************************
  juldate                     Keith Cherkauer             June 16, 1999

  This subroutine computes the julian date for the given decimal hour, 
  day, month, and year.  Julian dates from this routine start on 
  January 1, 1900.  Subtract 29220 to start on January 1, 1980, and 
  32873 to start on January 1, 1990.
***********************************************************************/
{
  double date, calc1, calc2, calc3, calc4;
  double hour, min, sec;

  hour = floor(time);
  time = time - (hour);
  min  = floor(time * 60.);
  time = time - (min / 60.);
  sec  = floor(time * 3600);
  time = time - (min / 3600);
  hour = hour + (min/60) + (sec/3600);

  calc1 = (367 * (double)year);
  calc2 = abs(floor(((double)month + 9)/12)); 
  calc2 = abs(floor(7*((double)year + calc2)/4));
  calc3 = abs(floor(275 * (double)month/9));
  calc4 = ((double)day + (hour/24.0) - 694006);
  
  date = calc1 - calc2 + calc3 + calc4;

  return(date);

}

int calc_doy( DATE_STRUCT self ) {
  // Routine to calculate the doy of the year for the given year, month, day and hour
  double juldate;
  juldate = calc_juldate( self.year, 1, 1, 0 );
  return ( (int)(self.juldate-juldate+1) );
}

DATE_STRUCT get_juldate( DATE_STRUCT self ) {
  // Routine to calculate julian date for the given year, month, day and hour
  self.juldate = calc_juldate( self.year, self.month, self.day, self.hour );
  return ( self );
}

DATE_STRUCT get_season( DATE_STRUCT self ) {
  // Routine to determine the season for the given date
  int    idx;

  idx = 0;
  while ( self.month != SeasonNumbers[idx] && idx < 12 ) idx++;
  if ( idx >= 12 ) {
    fprintf( stderr, "ERROR: invalid month, unable to determine season.\n" );
    exit (FAIL);
  }
  self.season = (int)( idx / 3 );
  return ( self );
}

DATE_STRUCT get_next_hour(DATE_STRUCT self) {
  // Routine returns a date, one hour after the current date
  if ( LEAPYR(self.year) ) DaysInMonth[1] = 29;
  else DaysInMonth[1] = 28;
  self.hour += 1.;
  if ( self.hour >= 24 ) {
    self.hour -= 24;
    self.day  += 1;
    if ( self.day > DaysInMonth[self.month-1] ) {
      self.day = 1;
      self.month = self.month + 1;
      if ( self.month > 12 ) {
	self.month = 1;
	self.year += 1;
      }
    }
  }
  return (self);
}

DATE_STRUCT get_last_hour( DATE_STRUCT self ) {
  // Routine returns a date, one hour before the current date
  if ( LEAPYR(self.year) ) DaysInMonth[1] = 29;
  else DaysInMonth[1] = 28;
  self.hour -= 1.;
  if ( self.hour < 0 ) {
    self.hour += 24;
    self.day  -= 1;
    if ( self.day < 1 ) {
      self.month -= 1;
      if ( self.month < 1 ) {
	self.month = 12;
	self.year -= 1;
	self.day =  DaysInMonth[self.month-1];
      }
    }
  }
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_next_day( DATE_STRUCT self ) {
  // Routine returns a date, one day after the current date
  if ( LEAPYR(self.year) ) DaysInMonth[1] = 29;
  else DaysInMonth[1] = 28;
  self.day += 1;
  if ( self.day > DaysInMonth[self.month-1] ) {
    self.day = 1;
    self.month += 1;
    if ( self.month > 12 ) {
      self.month = 1;
      self.year += 1;
    }
  }
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_last_day( DATE_STRUCT self ) {
  // Routine returns a date, one day before the current date
  self.day -= 1;
  if ( self.day < 1 ) {
    self.month -= 1;
    if ( self.month < 1 ) {
      self.month = 12;
      self.year -= 1;
      if ( LEAPYR(self.year) ) DaysInMonth[1] = 29;
      else DaysInMonth[1] = 28;
      self.day = DaysInMonth[self.month-1];
    }
  }
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_next_week( DATE_STRUCT self ) {
  // Routine returns a date, one week after the current date
  if ( LEAPYR(self.year) ) DaysInMonth[1] = 29;
  else DaysInMonth[1] = 28;
  self.day += 7;
  if ( self.day > DaysInMonth[self.month-1] ) {
    self.day = 1;
    self.month += 1;
    if ( self.month > 12 ) {
      self.month = 1;
      self.year += 1;
    }
  }
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_last_week( DATE_STRUCT self ) {
  // Routine returns a date, one week before the current date
  self.day -= 7;
  if ( self.day < 1 ) {
    self.month -= 1;
    if ( self.month < 1 ) {
      self.month = 12;
      self.year -= 1;
      if ( LEAPYR(self.year) ) DaysInMonth[1] = 29;
      else DaysInMonth[1] = 28;
      self.day = DaysInMonth[self.month-1];
    }
  }
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_next_month( DATE_STRUCT self ) {
  // Routine returns a date, one month after the current date"
  self.month += 1;
  if ( self.month > 12 ) {
    self.month = 1;
    self.year += 1;
  }
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_last_month( DATE_STRUCT self ) {
  // Routine returns a date, one month before the current date
  self.month -= 1;
  if ( self.month < 1 ) {
    self.month = 12;
    self.year -= 1;
  }
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_next_season( DATE_STRUCT self ) {
  // Routine returns the first day of the next season
  if ( self.season == 3 && self.month > 10 )
    self.year += 1;
  self.season ++;
  if ( self.season >= 4 ) self.season = 0;
  self.month = SeasonNumbers[self.season*3];
  self.day = 1;
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_last_season( DATE_STRUCT self ) {
  // Routine returns the first day of the previous season
  if ( self.season == 0 )
    self.year -= 1;
  self.season --;
  if ( self.season < 0 ) self.season = 3;
  self.month = SeasonNumbers[self.season*3];
  self.day = 1;
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_next_year( DATE_STRUCT self ) {
  // Routine returns a date, one year after the current date
  self.year += 1;
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_last_year( DATE_STRUCT self ) {
  // Routine returns a date, one year before the current date
  self.year -= 1;
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_next_month_start( DATE_STRUCT self ) {
  // Routine returns the month starting on or after the given date
  if ( self.day == 1 && self.hour == 0 ) return ( self );
  self.month += 1;
  if ( self.month > 12 ) {
    self.month = 1;
    self.year += 1;
  }
  self.day = 1;
  self.hour = 0;
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT get_next_season_start( DATE_STRUCT self ) {
  // Routine returns the season starting on or after the selected date
  int    idx;

  for ( idx = 0; idx < 4 ; idx ++ )
    if ( self.month == SeasonNumbers[idx*3] && self.day == 1 && self.hour == 0 )
      return ( self );
  if ( self.season == 3 && self.month > 10 )
    self.year += 1;
  self.season ++;
  if ( self.season >= 4 ) self.season = 0;
  self.month = SeasonNumbers[self.season*3];
  self.day = 1;
  self = get_juldate( self );
  return ( self );
}

DATE_STRUCT copy_date( DATE_STRUCT self ) {
  // creates copy of given date structure
  DATE_STRUCT new_date;

  new_date.day     = self.day;
  new_date.month   = self.month;
  new_date.year    = self.year;
  new_date.season  = self.season;
  new_date.doy     = self.doy;
  new_date.rec     = self.rec;
  new_date.hour    = self.hour;
  new_date.juldate = self.juldate;

  return ( new_date );

}

