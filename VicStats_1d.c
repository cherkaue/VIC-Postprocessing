#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <VicStats_1d.h>

/*********************************************************************
	stats.c		Keith Cherkauer		February 25, 1996

  This file contains routine for calculating statistics (mean, standard
  deviation, ...) for data in single dimension arrays of length N.

*********************************************************************/

double get_mean(double *values, int N, double NoData) {

  int index;
  int cnt;
  double mean=0.0;

  cnt = 0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      mean+=values[index];
      cnt++;
    }
  }

  if(cnt>0) mean /= (double)cnt;
  else mean = NoData;

  return(mean);

}

double get_var(double *values, int N, double mean, double NoData) {

  int index;
  int cnt;
  double stdev=0.0;

  cnt = 0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      stdev+=pow(values[index]-mean,2.0);
      cnt++;
    }
  }

  if(cnt>0) stdev = stdev/(double)(cnt-1);
  else stdev = NoData;

  return(stdev);

}

double get_stdev(double *values, int N, double mean, double NoData) {

  int index;
  int cnt;
  double stdev=0.0;

  cnt = 0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      stdev+=pow(values[index]-mean,2.0);
      cnt++;
    }
  }

  if(cnt>0) stdev = sqrt((stdev)/(double)(cnt-1));
  else stdev = NoData;

  return(stdev);

}

double get_skew(double *values, int N, double mean, double stdev, double NoData) {

  int index;
  int cnt;
  double skew=0.0;

  cnt = 0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      skew+=pow(values[index]-mean,3.0);
      cnt++;
    }
  }

  if(cnt>0) skew = (skew*(double)cnt)/(pow(stdev,3.0)*(double)(cnt-1)
                   *(double)(cnt-2));
  else skew = NoData;

  return(skew);

}

double get_sum(double *values, int N, double NoData) {

  int index;
  double sum;

  sum=0.0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      sum += values[index];
    }
  }
  
  return (sum);

}

double get_min(double *values, int N, double NoData) {

  int index;
  double min;

  min=values[0];
  for(index=1;index<N;index++)
    if(values[index]!=NoData)
      if(min>values[index]) min=values[index];
  
  return (min);

}

double get_max(double *values, int N, double NoData) {

  int index;
  double max;

  max=values[0];
  for(index=1;index<N;index++)
    if(values[index]!=NoData)
      if(max<values[index]) max=values[index];
  
  return (max);

}

int get_count_over_thres(double *values, double thres, int N, double NoData) {
  // computes the number of records over a threshold

  int index;
  int count;

  count=0;
  for(index=0;index<N;index++)
    if(values[index]!=NoData)
      if(values[index]>thres) count++;
  
  return (count);

}

double get_average_days_over_thres(double *values, double thres, int N, double NoData) {
  // computes the average number of days over a threshold

  int index;
  int sum;
  int count;
  int over;

  sum=0;
  count=0;
  over=0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      if(values[index]>thres) {
	sum++;
	over=1;
      }
      else {
	if ( over ) {
	  count++;
	  over=0;
	}
      }
    }
  }
  
  return ((double)sum/(double)count);

}

double get_average_value_over_thres(double *values, double thres, int N, double NoData) {
  // computes the average value of events over a threshold

  int index;
  int count;
  double sum;

  sum=0;
  count=0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      if(values[index]>thres) {
	sum+=(values[index]-thres);
	count++;
      }
    }
  }
  
  if ( count == 0 ) return (NoData);
  else return (sum/(double)count);

}

int get_consecutive_days_over_thres(double *values, double thres, int N, double NoData) {
  // computes the number of records under a threshold between first and last occurance identifeid from the middle of the record (e.g., growing sesason between first and last frost)

  int first;
  int last;

  for(first=N/2;first>=0;first--)
    if(values[first]!=NoData)
      if(values[first]<thres) break;
  
  for(last=N/2;last<N;last++)
    if(values[last]!=NoData)
      if(values[last]<thres) break;
  
  return (last-first);

}

int get_count_under_thres(double *values, double thres, int N, double NoData) {
  // computes the number of records under a threshold

  int index;
  int count;

  count=0;
  for(index=0;index<N;index++)
    if(values[index]!=NoData)
      if(values[index]<thres) count++;
  
  return (count);

}

double get_average_days_under_thres(double *values, double thres, int N, double NoData) {
  // computes the average number of days under a threshold

  int index;
  int sum;
  int count;
  int under;

  sum=0;
  count=0;
  under=0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      if(values[index]<thres) {
	sum++;
	under=1;
      }
      else {
	if ( under ) {
	  count++;
	  under=0;
	}
      }
    }
  }
  
  return ((double)sum/(double)count);

}

double get_average_value_under_thres(double *values, double thres, int N, double NoData) {
  // computes the average value of events under a threshold

  int index;
  double sum;
  int count;

  sum=0;
  count=0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      if(values[index]<thres) {
	sum+=thres-values[index];
	count++;
      }
    }
  }

  if ( count == 0 ) return (NoData);
  else return (sum/(double)count);

}

int get_consecutive_days_under_thres(double *values, double thres, int N, double NoData) {
  // computes the number of records under a threshold between first and last occurance identifeid from the middle of the record (e.g., growing sesason between first and last frost)

  int first;
  int last;

  for(first=N/2;first>=0;first--)
    if(values[first]!=NoData)
      if(values[first]>thres) break;
  
  for(last=N/2;last<N;last++)
    if(values[last]!=NoData)
      if(values[last]>thres) break;
  
  return (last-first);

}

int get_first_day_over_thres(double *values, double thres, int N, double NoData, DATE_STRUCT firstdate) {
  // find the first day value exceeds a threshold

  int index;

  for(index=0;index<N;index++) {
    if(values[index]!=NoData)
      if(values[index]>thres) break;
    firstdate = get_next_day( firstdate );
  }

  if ( index == N ) return( (int)NoData );
  else return ( calc_doy(firstdate) );

}

int get_first_day_over_thres_from_middle(double *values, double thres, int N, double NoData, DATE_STRUCT firstdate) {
  // find the first day value exceeds a threshold but start in the middle of the record (e.g. find first frost in calendar year)

  int index, idx;

  for(index=N/2;index<N;index++) {
    if(values[index]!=NoData)
      if(values[index]>thres) break;
  }

  if ( index == N ) return( (int)NoData );
  else {
    for(idx=0;idx<index;idx++) firstdate = get_next_day( firstdate );
    return ( calc_doy(firstdate) );
  }

}

int get_first_day_under_thres(double *values, double thres, int N, double NoData, DATE_STRUCT firstdate) {
  // find the first day value is below a threshold

  int index;

  for(index=0;index<N;index++) {
    if(values[index]!=NoData)
      if(values[index]<thres) break;
    firstdate = get_next_day( firstdate );
  }
  
  if ( index == N ) return( (int)NoData );
  else return ( calc_doy(firstdate) );

}

int get_first_day_under_thres_from_middle(double *values, double thres, int N, double NoData, DATE_STRUCT firstdate) {
  // find the first day value is below a threshold but start in the middle of the record (e.g. find first frost in calendar year)

  int index, idx;

  for(index=N/2;index<N;index++) {
    if(values[index]!=NoData) 
      if(values[index]<thres) break;
  }
  
  if ( index == N ) return( (int)NoData );
  else {
    for(idx=0;idx<index;idx++) firstdate = get_next_day( firstdate );
    return ( calc_doy(firstdate) );
  }

}

int get_last_day_over_thres(double *values, double thres, int N, double NoData, DATE_STRUCT firstdate) {
  // find the last day value is over a threshold

  int index;
  int last;
  DATE_STRUCT tmpdate;

  last = -99;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData)
      if(values[index]>thres) {
	last = index;
	tmpdate = copy_date( firstdate );
      }
    firstdate = get_next_day( firstdate );
  }
  
  if ( index == N ) return( (int)NoData );
  else return ( calc_doy(firstdate) );

}

int get_last_day_over_thres_from_middle(double *values, double thres, int N, double NoData, DATE_STRUCT firstdate) {
  // find the last day value is over a threshold but work back from the middle of the record (e.g. find last frost in calendar year)

  int index, idx;

  for(index=N/2;index>=0;index--) {
    if(values[index]!=NoData)
      if(values[index]>thres) break;
  }
  
  if ( index < 0 ) return( (int)NoData );
  else {
    for(idx=0;idx<index;idx++) firstdate = get_next_day( firstdate );
    return ( calc_doy(firstdate) );
  }

}

int get_last_day_under_thres(double *values, double thres, int N, double NoData, DATE_STRUCT firstdate) {
  // find the last day value is over a threshold

  int index;
  int last;
  DATE_STRUCT tmpdate;

  last = -99;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData)
      if(values[index]<thres) { 
	last = index;
	tmpdate = copy_date( firstdate );
      }
    firstdate = get_next_day( firstdate );
  }
  
  if ( index == N ) return( (int)NoData );
  else return ( calc_doy(firstdate) );

}

int get_last_day_under_thres_from_middle(double *values, double thres, int N, double NoData, DATE_STRUCT firstdate) {
  // find the last day value is over a threshold but work back from the middle of the record (e.g. find last frost in calendar year)

  int index, idx;

  for(index=N/2;index>=0;index--) {
    if(values[index]!=NoData)
      if(values[index]<thres) break;
  }
  
  if ( index < 0 ) return( (int)NoData );
  else {
    for(idx=0;idx<index;idx++) 
      firstdate = get_next_day( firstdate );
    return ( calc_doy(firstdate) );
  }

}

int get_count_times_thres_crossed(double *values, double thres, int N, double NoData) {
  // find the number of times a threshold is crossed

  int index;
  int count;
  char HIGH = 0;

  if ( values[0] > thres ) HIGH = 1;
  for(index=0;index<N;index++)
    if(values[index]!=NoData)
      if(values[index]>thres && !HIGH) {
	count++;
	HIGH=1;
      }
      else if(values[index]<thres && HIGH) {
	count++;
	HIGH=0;
      }

  return (count);

}

int get_modified_growing_degree_day(double *TMIN, double *TMAX, int N, double NoData) {
  // find the seasonal growing degree day count based on MRCC definition

  int first;
  int last;
  int index;
  int count;
  int GDD;

  // find start of growing season searching back from middle of record
  for(first=N/2;first>=0;first--)
    if(TMIN[first]!=NoData)
      if(TMIN[first]<0) break;

  // find end of growing season searching back from middle of record
  for(last=N/2;last<N;last++)
    if(TMIN[last]!=NoData)
      if(TMIN[last]<0) break;

  // compute growing degree day from start of growing season
  for (index=first;index<last;index++) {
    if(TMIN[index]!=NoData && TMAX[index]!=NoData)
      if(TMIN[index]<10) TMIN[index]=10;
  }

  return (count);

}

double get_RB_index(double *values, int N, double NoData) {
  // find the Richards-Baker flashiness index

  int index;
  int start;
  double sum, diff;

  for(start=0;start<N;start++) if (values[start]!=NoData) break;
 
  if ( start >= N-1 ) return( NoData ); // no valid data

  for(index=start+1;index<N;index++) {
    if(values[index]!=NoData) {
      sum += values[index];
      diff += fabs( values[index] - values[index-1] );
    }
  }

  return (diff/sum);

}

double get_TQ_mean(double *values, int N, double NoData) {
  // find the TQ mean of total runoff

  int index;
  int count;
  double mean;

  mean = get_mean( values, N, NoData );

  count=0;
  for(index=0;index<N;index++)
    if(values[index]!=NoData)
      if(values[index]>mean) count++;

  return ((double)count/(double)N);

}

double get_seven_day_low_value(double *values, int N, double NoData) {
  // get the lowest value out of a 7-day moving average.

  int index;
  int count;
  int Nval;
  double *avg_values;
  double min_value;

  min_value = 1.e99;
  avg_values = (double *)calloc(sizeof(double),(N-7));
  for(index=0;index<N-7;index++) {
    Nval=0;
    for(count=0;count<7;count++) {
      if(values[index+count]!=NoData) {
	avg_values[index] += values[index+count];
	Nval++;
      }
    }
    if(Nval>0) {
      avg_values[index] /= (double)Nval;
      if ( avg_values[index] < min_value ) min_value = avg_values[index];
    }
    else avg_values[index] = NoData;
  }

  return (min_value);

}

double get_quantile(double *values, double quantile, int N, double NoData) {
  // find the value associated with the provided quantile when values are ranked
  // low to high.

  int index;
  int start;
  double PP;
  
  // sort data from smallest to largest
  quickSort(values, 0, N-1);

  // skip over no data values, which are now lowest values
  for(start=0;start<N;start++) if (values[start]!=NoData) break;
 
  // find desired quantile
  for ( index=start; index<N; index++ ) {
    PP = ( 1. - (double)(index-start) / (double)(N-start+1) ) * 100.;
    if ( PP <= quantile ) break;
  }
  
  return (values[index]);

}

void get_3LN(double mean,double stdev,double skew,double *tau,double *muy,double *sigmay) {

  double GAMMA,phi,epsilon;

  GAMMA = skew*sqrt(4.0*pow(skew,2.0));
  phi = pow(1.0+0.5*(pow(skew,2.0)+GAMMA),1.0/3.0)+pow(1.0+0.5*(pow(skew,2.0)-GAMMA),1.0/3.0)-1;
  epsilon = sqrt(pow(stdev,2.0)/phi/(phi-1));
  *tau = mean-epsilon*sqrt(phi);
  *muy = log(epsilon);
  *sigmay = sqrt(log(phi));

  return;

}

double get_rk(int N, int k, double *season1, double *season2, double mean1, double mean2) {
/**********************************************************************
  This function computes the lag (k) correlation coefficient between
  two seasons of normalized streamflow data.
**********************************************************************/
  double rk=0.0,sum1=0.0,sum2=0.0;
  int r;

  for(r=0;r<N-k;r++) {
    rk += (season1[r+k] - mean1)*(season2[r] - mean2);
    sum1 += pow(season1[r+k] - mean1,2.0);
    sum2 += pow(season2[r] - mean2,2.0);
  }

  sum1 = sqrt(sum1/(double)(N-k));
  sum2 = sqrt(sum2/(double)(N-k));

  rk = rk/(double)(N-k)/(sum1*sum2);

  return (rk);

}

/* C implementation QuickSort from http://www.geeksforgeeks.org/quick-sort/ */
 
// A utility function to swap two elements
void swap(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}
 
/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partition (double *arr, int low, int high)
{
    double pivot = arr[high];    // pivot
    int i = (low - 1);  // Index of smaller element
 
    for (int j = low; j <= high- 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[j] <= pivot)
        {
            i++;    // increment index of smaller element
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}
 
/* The main function that implements QuickSort
 arr[] --> Array to be sorted,
  low  --> Starting index,
  high  --> Ending index */
void quickSort(double *arr, int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition(arr, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}
 
