#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <zlib.h>

/** User defined variables **/
#define MAXLINE 256 // maximum number of characters on a forcing file line
#define Nplaces  4         /* number of decimal places to use in file names */
#define NODATA_val   -9999
#define MAX_VARS 25
/** End User Defined Variables **/

#ifndef _LEAPYR
#define LEAPYR(y) (!((y)%400) || (!((y)%4) && ((y)%100)))
#endif

int read_arcinfo_grid(char *, double **, double **, int *, 
		      int *, double *, double *, double *, int *, double **);
int write_arcinfo_grid(char *, double *, double *, double *, int,
		       int, double, double, double, int, int);
double get_mean(double *, int, double);
double get_stdev(double *, int, double, double);
double get_sum(double *, int, double);
double get_min(double *, int, double);
double get_max(double *, int, double);
double juldate( int, int, int, double );

#define FALSE 0
#define TRUE  !FALSE

int main(int argc, char *argv[]) {
/**********************************************************************
  compute_average_grid_forcings.c  Keith Cherkauer    April 14, 1999

  This program was written to compute the average basin forcing date
  for the period of interest.  Basin forcing data should be stored in
  gridded short int binary files.  Active grid cells are determined 
  from the given ARC/INFO mask file.  All output statistics are written
  to ARC/INFO ASCII grid files.

  Modifications:
  01/12/2011 Modified to increase flexability, including allowing the 
    user to define parameters and file format at the command line.  And
    allowing for the processing of ASCII column files.            KAC

**********************************************************************/

  gzFile    *fin;

  char       tmpstr[512];
  char       junk[30];
  char       forcing_prefix[1024];
  char       mask_file_name[1024];
  char       output_path[1024];
  char       ASCII;
  char       LineStore[MAXLINE];
  char      *VarNameList[MAX_VARS];
  char       VarStatType[MAX_VARS], VarSignList[MAX_VARS];
  char      *TmpStr;
  int        VarColNum[MAX_VARS], NumCols, Nparam, varcnt, col, TAVG[3];
  int       *ColParam;
  double     VarMultList[MAX_VARS];
  double    *lat;
  double    *lng;
  double    *tmp;
  double   **values;
  double   **max, **min, **mean, **stdev;
  double    *tmp_ptr;
  double     ll_lat;
  double     ll_lng;
  double     cellsize;
  double     tmp_max, tmp_min;
  int        ncols;
  int        nrows;
  int        NODATA;
  int        Ncells;
  int        i;
  int        cell, rec, param;
  int        curryear;
  int        Ndays, Nyears;
  int        filestartdate[4], startdate[4], enddate[4];
  int        NumDays, SkipBytes;
  double     filestartjday, startjday, endjday, currjday, nextjday;
  short int  si_val;
  unsigned short int usi_val;
  char       DONE;

  // Check usage
  if ( argc != 10 ) {

    fprintf( stderr, "\nUsage: %s <Forcing Prefix> <Mask File> <Output Path> <File Start Date> <Average Start Date> <Average End Date> <ASCII|BINARY> <Num Columns> <Var1/Col/<S|A>[/Mult/<S|U>],<Var2/Col/<S|A>[/Mult/<S|U>],....<VarN/Col/<S|A>[/Mult/<S|U>]>\n", argv[0] );
    fprintf( stderr, "\n\tThis program was written to compute the average basin forcing date \n\tfor the period of interest.  Basin forcing data should be stored in \n\tgridded short int binary files.  Active grid cells are determined \n\tfrom the given ARC/INFO mask file.  All output statistics are written \n\tto ARC/INFO ASCII grid files.  Number and type of columns is \n\tdefined in the source code and changes must be recompiled.\n" );
    fprintf( stderr, "\n\t<Forcing Prefix> prefix, including directory, of the forcing files \n\t\t(e.g. $PATH/data_ ),\n" );
    fprintf( stderr, "\t<Mask File> ArcInfo ASCII grid file with active cells indicated as \n\t\tdata,\n" );
    fprintf( stderr, "\t<Output Path> path to which output grids will be written,\n" );
    fprintf( stderr, "\t<Start Date> and <End Date> are in the form MMDDYYYY (MM = month, \n\t\tDD = day, YYYY = year) and computations start at midnight \n\t\tthe morning of that day.\n");
    fprintf( stderr, "\t<ASCII|BINARY> indicates whether the raw file is ASCII columns, or \n\t\tshort integer binary,\n" );
    fprintf( stderr, "\t<Var1/Col/<S|A>[/Mult/<S|U>] each variable must provide a name <VarX>, a \n\t\tcolumn number <Col>, and whether the variable should be <S>ummed \n\t\tor <A>veraged.  If data is binary you must also include the \n\t\tmultiplier (use 1, for no change to the value) and whether the \n\t\tshort int number is <S>igned or <U>nsigned.  If the variables \n\t\tTMAX and TMIN are defined, TAVG will also be output.\n\n" );
    exit(-1);

  }

  // Handle command line arguments
  strcpy( forcing_prefix, argv[1] );
  strcpy( mask_file_name, argv[2] );
  strcpy( output_path, argv[3] );
  
  filestartdate[1] = (int)(atof(argv[4])/1000000);
  filestartdate[2] = (int)(atof(argv[4])/10000) - filestartdate[1]*100;
  filestartdate[0] = (int)(atof(argv[4])) - filestartdate[1]*1000000 
    - filestartdate[2]*10000;
  filestartdate[3] = 0;
  filestartjday    = juldate(filestartdate[0],filestartdate[1],filestartdate[2],filestartdate[3]);
  startdate[1] = (int)(atof(argv[5])/1000000);
  startdate[2] = (int)(atof(argv[5])/10000) - startdate[1]*100;
  startdate[0] = (int)(atof(argv[5])) - startdate[1]*1000000 
    - startdate[2]*10000;
  startdate[3] = 0;
  startjday    = juldate(startdate[0],startdate[1],startdate[2],startdate[3]);
  enddate[1]   = (int)(atof(argv[6])/1000000);
  enddate[2]   = (int)(atof(argv[6])/10000) - enddate[1]*100;
  enddate[0]   = (int)(atof(argv[6])) - enddate[1]*1000000 
    - enddate[2]*10000;
  enddate[3]   = 0;
  endjday      = juldate(enddate[0],enddate[1],enddate[2],enddate[3]);
  NumDays      = (int)(endjday-startjday);
  if ( strcasecmp( argv[7], "ASCII" ) == 0 ) ASCII = TRUE;
  else ASCII = FALSE;
  NumCols = atoi( argv[8] );
  ColParam = (int *)calloc(sizeof(int),NumCols);
  // process variable definition string starting with comma separations
  varcnt = 0;
  VarNameList[varcnt] = strtok( argv[9], "," );
  while ( VarNameList[varcnt] != NULL ) {
    varcnt++;
    VarNameList[varcnt] = strtok( '\0', "," );
  }
  Nparam = varcnt;
  TAVG[0] = 0;
  // process each variable string as separated with slashes
  for ( varcnt = 0; varcnt < Nparam; varcnt++ ) {
    VarNameList[varcnt] = strtok( VarNameList[varcnt], "/" );
    VarColNum[varcnt] = atoi(strtok( '\0', "/" )) - 1;
    VarStatType[varcnt] = strtok( '\0', "/" )[0];
    if ( ! ASCII ) {
      VarMultList[varcnt] = atof(strtok( '\0', "/" ));
      VarSignList[varcnt] = strtok( '\0', "/" )[0];
    } else {
      VarMultList[varcnt] = 1.0;
      VarSignList[varcnt] = 'S';
    }
    ColParam[VarColNum[varcnt]] = varcnt+1; // link column with parameter
    // determine if TAVG can be calculated
    if ( strcasecmp( VarNameList[varcnt], "TMAX" ) == 0 ) {
      TAVG[0] ++;
      TAVG[1] = varcnt;
    } else if ( strcasecmp( VarNameList[varcnt], "TMIN" ) == 0 ) {
      TAVG[0] ++;
      TAVG[2] = varcnt;
    }
    VarStatType[Nparam] = 'A';
    printf( "%s\t%i\t%c\t", VarNameList[varcnt], VarColNum[varcnt], VarStatType[varcnt] );
    if ( ASCII ) printf("\n");
    else printf( "%f\t%c\n", VarMultList[varcnt], VarSignList[varcnt] );
  }
  if ( TAVG[0] == 2 ) {
    printf( "TAVG will be calculated.\n" );
    strcpy( VarNameList[Nparam], "PRCP" );
    VarColNum[Nparam] = 0;
    VarMultList[Nparam] = 1.0;
    VarSignList[Nparam] = 'A';
  }

  // set number of bytes to skip at the start of the forcing file
  currjday = filestartjday;
  if ( ASCII ) SkipBytes = ( startjday - currjday ); // lines not bytes
  else SkipBytes = Nparam * sizeof(short int) * ( startjday - currjday );

  // Allocate time series arrays
  values = (double **)calloc(Nparam+1,sizeof(double *));
  for ( i = 0; i <= Nparam; i++ ) 
    values[i] = (double *)calloc(NumDays,sizeof(double));

  // Read mask file
  Ncells = read_arcinfo_grid(mask_file_name,&lat,&lng,&ncols,&nrows,
			     &ll_lat,&ll_lng,&cellsize,&NODATA,&tmp);
  fprintf(stdout,"Processing %i Cells.\n",Ncells);

  // Allocate statistic arrays
  mean  = (double **)calloc(Nparam+1,sizeof(double *));
  stdev = (double **)calloc(Nparam+1,sizeof(double *));
  max   = (double **)calloc(Nparam+1,sizeof(double *));
  min   = (double **)calloc(Nparam+1,sizeof(double *));
  for ( i = 0; i <= Nparam; i++ ) {
    mean[i]  = (double *)calloc(Ncells,sizeof(double));
    stdev[i] = (double *)calloc(Ncells,sizeof(double));
    max[i]   = (double *)calloc(Ncells,sizeof(double));
    min[i]   = (double *)calloc(Ncells,sizeof(double));
  }

  sprintf(junk, "%%s%%.%ilf_%%.%ilf", Nplaces, Nplaces);

  /** open each file and compute stats **/
  for ( cell = 0; cell < Ncells; cell++ ) {

    sprintf(tmpstr,junk,forcing_prefix,lat[cell],lng[cell]);

    fprintf(stdout,"Processing %s\n",tmpstr);

    /** Open forcing file **/
    if((fin=gzopen(tmpstr,"rb"))==NULL) {
      if((fin=gzopen(strcat(tmpstr,".gz"),"rb"))==NULL) {
	fprintf(stderr,"ERROR: Unable to open forcing file %s\n",tmpstr);
	exit(0);
      }
    }

    // Reinitialize average temperature array
    if ( TAVG[0] == 2 )
      for ( rec = 0; rec < NumDays; rec++ ) values[Nparam][rec] = 0;

    /** Read forcing data for current cell **/
    if ( ASCII ) {
      // read data from ASCII column file
      // skip start of file to get to requested time period
      for ( i = 0; i < SkipBytes; i++ )
	if ( gzgets(fin,LineStore,MAXLINE) == Z_NULL ) {
	  fprintf(stderr, "ERROR: End of forcing file %s found unexpectedly!\n", tmpstr );
	  return(-1);
	}

      // process all days in given file
      for(rec = 0; rec < NumDays; rec++ ) {
	if ( gzgets(fin,LineStore,MAXLINE) == Z_NULL ) {
	  fprintf(stderr, "ERROR: End of forcing file %s found unexpectedly!\n", tmpstr );
	  return(-1);
	}
	for ( col = 0; col < NumCols; col++ ) {
	  if ( col == 0 ) TmpStr = strtok( LineStore, " \t" );
	  else TmpStr = strtok( '\0', " \t" );
	  if ( ColParam[col] != 0 ) {
	    values[ColParam[col]-1][rec] = atof(TmpStr);
	  }
	}
	if ( TAVG[0] == 2 )
	  // compute average temperature value
	  values[Nparam][rec] = ( values[TAVG[1]][rec] + values[TAVG[2]][rec] ) / 2.;
      }
    } else {
      // read data from a short int binary forcing file

      // skip start of file to get to requested time period
      gzseek(fin,SkipBytes,SEEK_SET);

      // process all days in given file
      for(rec = 0; rec < NumDays; rec++ ) {
	for ( col = 0; col < NumCols; col++ ) {
	  if ( ColParam[col] != 0 ) {
	    if ( VarSignList[ColParam[col]-1] == 'S' ) {
	      // read signed short int value
	      gzread(fin,&si_val,1*sizeof(short int));
	      values[ColParam[col]-1][rec] = (double)si_val / VarMultList[ColParam[col]-1];
	    }
	    else {
	      // read unsigned short int
	      gzread(fin,&usi_val,1*sizeof(unsigned short int));
	      values[ColParam[col]-1][rec] = (double)usi_val / VarMultList[ColParam[col]-1];
	    }
	  } else {
	    // read signed short int value by default
	    gzread(fin,&si_val,1*sizeof(short int));
	    values[ColParam[col]-1][rec] = (double)si_val / VarMultList[ColParam[col]-1];
	  }
	}
	if ( TAVG[0] == 2 )
	  // compute average temperature value
	  values[Nparam][rec] = ( values[TAVG[1]][rec] + values[TAVG[2]][rec] ) / 2.;
      }
    }

    curryear = startdate[0];
    nextjday = startjday;
    
    DONE     = FALSE;
    Nyears   = 0;
    rec      = 0;

    // Compute annual statistics
    while (!DONE) {

      currjday = nextjday;
      nextjday = juldate(curryear+1,startdate[1],startdate[2],startdate[3]);
      Ndays = (int)( nextjday - currjday + 0.5 );
      
      if ( rec+Ndays <= NumDays ) {
	for ( param = 0; param <= Nparam; param++ ) {
	  
	  // Compute statistics for parameter data
	  tmp_ptr = &values[param][rec];
	  if( VarStatType[param] == 'S' )
	    // annual sum of values
	    mean[param][cell] += get_sum(tmp_ptr,Ndays,NODATA_val);
	  else
	    // annual average values
	    mean[param][cell] += get_mean(tmp_ptr,Ndays,NODATA_val);
	  // standard deviation of values
	  stdev[param][cell] += get_stdev(tmp_ptr,Ndays,
					  get_mean(tmp_ptr,
						   Ndays,NODATA),NODATA_val);
	  // find maximum value
	  tmp_max = get_max(tmp_ptr,Ndays,NODATA_val);
	  if(tmp_max > max[param][cell] || rec==0) max[param][cell] = tmp_max;
	  // find minimum values
	  tmp_min = get_min(tmp_ptr,Ndays,NODATA_val);
	  if(tmp_min < min[param][cell] || rec==0) min[param][cell] = tmp_min;
	  
	}
 
	rec += Ndays;
	curryear ++;
	Nyears ++;
	
      }
      else {
	// Finished computing annual statistics, average all years
	DONE = TRUE;
	printf("  --> completed %i years\n",Nyears);
	for(param=0;param<=Nparam;param++) {
	  mean[param][cell] /= (double)Nyears;
	  stdev[param][cell] /= (double)Nyears;
	}
      }
    }
    gzclose(fin);
  }
  
  /** write statistics for entire basin **/
  for(param=0;param<=Nparam;param++) {
    sprintf(tmpstr,"%s/%s_mean.asc",output_path,VarNameList[param]);
    cell = write_arcinfo_grid(tmpstr,mean[param],lat,lng,ncols,nrows,
			      ll_lat,ll_lng,cellsize,NODATA,Ncells);
    sprintf(tmpstr,"%s/%s_stdev.asc",output_path,VarNameList[param]);
    cell = write_arcinfo_grid(tmpstr,stdev[param],lat,lng,ncols,nrows,
			      ll_lat,ll_lng,cellsize,NODATA,Ncells);
    sprintf(tmpstr,"%s/%s_max.asc",output_path,VarNameList[param]);
    cell = write_arcinfo_grid(tmpstr,max[param],lat,lng,ncols,nrows,
			      ll_lat,ll_lng,cellsize,NODATA,Ncells);
    sprintf(tmpstr,"%s/%s_min.asc",output_path,VarNameList[param]);
    cell = write_arcinfo_grid(tmpstr,min[param],lat,lng,ncols,nrows,
			      ll_lat,ll_lng,cellsize,NODATA,Ncells);
  }

  return(0);

}

int read_arcinfo_grid(char    *filename,
		      double **lat,
		      double **lng,
		      int     *ncols,
		      int     *nrows,
		      double  *ll_lat,
		      double  *ll_lng,
		      double  *cellsize,
		      int     *NODATA,
		      double **values) {
/**********************************************************************
  read_arcinfo_info           Keith Cherkauer           May 5, 1998

  This subroutine reads data from an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

**********************************************************************/

  FILE   *farc;
  int    i, j;
  double tmp_lat;
  double tmp_lng;
  int    cell;
  int    Ncells;
  int    tmpvalue;

  if((farc=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable of open ARCINFO grid %s for reading.\n",filename);
    exit(0);
  }

  /***** Read ARC/INFO Header *****/
  fscanf(farc,"%*s %i",&ncols[0]);
  fscanf(farc,"%*s %i",&nrows[0]);
  fscanf(farc,"%*s %lf",&ll_lng[0]);
  fscanf(farc,"%*s %lf",&ll_lat[0]);
  fscanf(farc,"%*s %lf",&cellsize[0]);
  fscanf(farc,"%*s %i",&NODATA[0]);

  /***** Allocate Latitude and Longitude Arrays for maximum size *****/
  Ncells    = ncols[0]*nrows[0];
  lat[0]    = (double *)calloc(Ncells,sizeof(double));
  lng[0]    = (double *)calloc(Ncells,sizeof(double));
  values[0] = (double *)   calloc(Ncells,sizeof(double));

  /***** Check for Valid Location *****/
  cell = 0;
  for(j=0;j<nrows[0];j++) {
    tmp_lat = ll_lat[0]+(double)(nrows[0]-j-0.5)*cellsize[0];
    for(i=0;i<ncols[0];i++) {
      tmp_lng = ll_lng[0] + (double)(i+0.5)*cellsize[0];
      fscanf(farc,"%i",&tmpvalue);
      if(tmpvalue!=NODATA[0]) {
	lat[0][cell]  = tmp_lat;
	lng[0][cell]  = tmp_lng;
	values[0][cell] = tmpvalue;
	cell++;
      }
    }
  }
  fclose(farc);
  Ncells = cell;

  return Ncells;

}

int write_arcinfo_grid(char   *filename,
		       double *values,
		       double *lat,
		       double *lng,
		       int     ncols,
		       int     nrows,
		       double  ll_lat,
		       double  ll_lng,
		       double  cellsize,
		       int     NODATA,
		       int     Ncells) {
/**********************************************************************
  write_arcinfo_info       Keith Cherkauer           April 14, 1999

  This subroutine reads data from an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

**********************************************************************/

  FILE   *farc;
  int    i, j;
  double tmp_lat;
  double tmp_lng;
  int    cell;

  if((farc=fopen(filename,"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable of open ARCINFO grid %s for writing.\n",filename);
    exit(0);
  }

  /***** Write ARC/INFO Header *****/
  fprintf(farc,"ncols\t%i\n",ncols);
  fprintf(farc,"nrows\t%i\n",nrows);
  fprintf(farc,"xllcorner\t%lf\n",ll_lng);
  fprintf(farc,"yllcorner\t%lf\n",ll_lat);
  fprintf(farc,"cellsize\t%lf\n",cellsize);
  fprintf(farc,"NODATA_value\t%i\n",NODATA);

  /***** Check for Valid Location *****/
  cell = 0;
  for(j=0;j<nrows;j++) {
    tmp_lat = ll_lat+(double)(nrows-j-0.5)*cellsize;
    for(i=0;i<ncols;i++) {
      tmp_lng = ll_lng + (double)(i+0.5)*cellsize;
      if(lat[cell] == tmp_lat && lng[cell] == tmp_lng) {
	fprintf(farc,"%f",values[cell]);
	cell++;
      }
      else
	fprintf(farc,"%i",NODATA);
      if(i<ncols-1) fprintf(farc,"\t");
      else fprintf(farc,"\n");
    }
  }
  fclose(farc);
  if(Ncells != cell) {
    fprintf(stderr,"ERROR: number of cells printed does not equal number of cells defined.\n");
    exit(0);
  }

  return cell;

}

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

double juldate(int year,
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

