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

int get_first_day_over_thres(double *values, double thres, int N, double NoData) {
  // find the first day value exceeds a threshold

  int index;

  for(index=0;index<N;index++)
    if(values[index]!=NoData)
      if(values[index]>thres) break;
  
  return (index);

}

int get_first_day_under_thres(double *values, double thres, int N, double NoData) {
  // find the first day value is below a threshold

  int index;

  for(index=0;index<N;index++)
    if(values[index]!=NoData)
      if(values[index]<thres) break;
  
  return (index);

}

int get_last_day_over_thres(double *values, double thres, int N, double NoData) {
  // find the last day value is over a threshold

  int index;
  int last;

  last = -99;
  for(index=0;index<N;index++)
    if(values[index]!=NoData)
      if(values[index]>thres) last = index;
  
  return (last);

}

int get_last_day_under_thres(double *values, double thres, int N, double NoData) {
  // find the last day value is over a threshold

  int index;
  int last;

  last = -99;
  for(index=0;index<N;index++)
    if(values[index]!=NoData)
      if(values[index]<thres) last = index;
  
  return (last);

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

