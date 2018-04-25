#include <VicStuff.h>

// get information from the staistics control file
int GetStatInfo( char *, StatInfoStruct *, char **, int, int *, int *, int *, int *, int *, int *, int *, int *, int *, char, double *, double *, int, int, double, double, double, int );
// get special information for calculating supplemental output variables
int GetPenmanInfo( PenInfoStruct *, char **, int );
int GetTotalRunoffInfo( PenInfoStruct *, char **, int );
int GetTotalSoilMoistureInfo( PenInfoStruct *, char **, int, int );
int GetModifiedGrowingDegreeDayInfo( PenInfoStruct *, StatInfoStruct *, char **, int );
int GetPlantingDayInfo( PenInfoStruct *, StatInfoStruct *, char **, int );
int GetChillingHoursInfo( PenInfoStruct *, StatInfoStruct *, char **, int );
int GetBudFreezeInfo( PenInfoStruct *, StatInfoStruct *, char **, int );
int GetDrainModWorkingDaysInfo( PenInfoStruct *, StatInfoStruct *, char **, int );
int GetDaysSuitable4FieldWorkInfo( PenInfoStruct *, StatInfoStruct *, char **, int );
// Functions for computing supplemental output variables
int CalcModifiedGrowingDegreeDays( double *, double *, double *, DATE_STRUCT, int );
int CalcPlantingDay( double *, double *, double *, double, DATE_STRUCT, int );
int CalcChillingHours( double *, double *, double *, DATE_STRUCT, int );
int CalcBudFreeze( double *, double *, double *, DATE_STRUCT, int );
int CalcDrainModWorkingDays( double *, double *, double *, double, DATE_STRUCT, int );
int CalcDaysSuitable4FieldWork( double *, double *, double *, double *, DATE_STRUCT, int );
// File handling routines
char CheckExistingOutputFiles ( char, char, char *, StatInfoStruct, char ***, 
				DATE_STRUCT, int, int, int, int, int, 
				double, double, double, int, int, char );
char CheckExistingArcInfoGridFiles ( char, char *, StatInfoStruct, char ***, 
				     DATE_STRUCT, int, int, int *, int *, 
				     double *, double *, double *, int * );
char CheckExistingXyzFiles ( char, char *, StatInfoStruct, char ***, DATE_STRUCT,
			     int, int, int *, int *, double *, double *, double *,
			     int * );
int read_arcinfo_grid(char *, double **, double **, int *, int *, double *, 
		      double *, double *, int *, double ***);
int write_arcinfo_grid(char *, double **, double *, double *, int, int, double, 
		       double, double, int, int);
int write_blank_arcinfo_grid(char *, int, int, double, double, double, int);
int WriteOutputFiles ( char, char *, double **, double *, double *,
		       int, int, double, double, double, int, int, char );
int WriteXyzFiles( char *, double **, double *, double *, int, int,
		   double, double, double, int, int );
int ReadXyzFiles(char *, double **, double **, int *, int *, double *, 
		      double *, double *, int *, double ***);
int ReadSpatialDataFiles ( char, char *, double **, double *, double *, int,
			   int, double, double, double, int);

int main(int argc, char *argv[]) {
/************************************************************************
  get_arc_mean_grid.c     Keith Cherkauer          April 16, 2001

  This program was written to extract summary statistics from VIC model
  daily output files.  The program requires that you give it a file 
  extension list for the output files to be included in the grid,
  the directory which contains the output, the file prefix of the output
  files (so the program can be run on fluxes, snow, snow_band, and fdepth
  files), the time step and date range for the statistics being computed,
  and a file that defines the columns to be assessed and the statistics
  to be computed on them.  Multiple statistics can be computed on any one 
  column with a single pass of this program.

  Based on VicBinary2ArcStatGrids, conversion compelted January 13, 2017.
  This code now works with ASCII and Binary VIC model output files, as long
  as output files are written with the output file header implemented in 
  VIC 4.1.0 r4 and available in all later versions of VIC 4.1.x when set
  in the global file.

  Program can output ESRI (ArcInfo) ASCII raster files and XYZ files.

  Modifications:
  2017-Mar-27 Modified to skip over cell files that fall outside of the 
    previously defined extent of the output ArcInfo raster file.  
    Previously this threw an error and cause the program to end, even when 
    valid cell files were found within the extent of the raster domain.  KAC
  2017-Apr-18 Updated calculation of number of seasons to output by 
    looping with get_next_season rather than trying to estimate the
    number of times 4 periods occurs within a year.  Produces more 
    reliable output from Seasonal statistics, as previously it was 
    crashing when it could not open files past the end of the data.  KAC
  2017-05-17 Fixed accounting so that when end of analysis period 
    equals the end of the file, the program will not crash with a 
    Segementation Fault.                                             KAC
  2017-05-23 Added additional metrics to help with Agricultural metrics
    for the IN-CCIA report.                                         KAC
  2017-05-29 Added statistic type "P" to define a fixed period within
    an annual time step (such as a growing season) over which to 
    compute the requested metric.  Period is defined with starting and
    ending month and day combinations, and data outside this range are
    set to NoData values before calculating annual metrics.      KAC
  2017-05-31 Added Chilling Hours and modified Growing Degree Days as 
    variables that are computed before final statistics for each period.
    This is different than PE and TOTAL_RUNOFF, which are computed as 
    the variables are read from the input files.  These variables requried 
    access to multiple columns of input data, not all of which were 
    necessarily requested for processing, so tose variables had to be
    read in and then calculations were made once the entire period was
    ready for statistical analysis.                      KAC
  2017-05-31 Modified first and last date functions to provide the day
    of year (doy) rather than the number of days since the start of the 
    analysis period.  This should provide results that are easier to 
    assess and explain.                                    KAC
  2017-05-31 Completed transition to supporting both ArcGIS grid and 
    XYZ file output formats.  This had been left incomplete.  KAC
  2017-06-12 Added the ability to define an input file rather than a 
    constant value for the threshold.    KAC
  2017-Jun-16 Added OUT_TOTAL_SOIL_MOIST to provide a sum of soil
    moisture in all layers.  Currently only works with three layer 
    simulations.  Probably also need to add OUT_TOTAL_SOIL_ICE.  KAC
  2017-Jun-16 Added new statistics to compute the average value of
    events over or under a threshold AVG(event - threshold).    KAC
  2017-Jun-20 Additional variable calculations have been added, and
    TOTAL_SOIL_MOISTURE has been updated to work with any number of
    soil layers, as defined in the VIC model output header.   KAC

  NOTE: As of 2017-Jun-12 - Support for XYZ files has not been checked.  KAC
  NOTE: As of 2017-Jun-12 - Support for no Overwrite is not completed.  KAC
  NOTE: Should be able to double check that a file has been given for 
    OUT_PLANT_DAY and report that the special case has been violated.
    Should also be able to catch a blank line in the stat control file,
    and gracefully skip it rather than throwing a segementation fault.  KAC

***********************************************************************/

  StatInfoStruct StatInfo;

  // variables that must be adjusted when adding a new calculated column
  PenInfoStruct PenInfo, TRoffInfo, TSMInfo, MGDDInfo, PDayInfo, ChillHrInfo;
  PenInfoStruct BudFreezeInfo, DModWDaysInfo, DSFWInfo;
  int     CalcPE, CalcTR, CalcTSM, CalcMGDD, CalcPD, CalcCH, CalcBF;
  int     CalcDMWD, CalcDSFW;
  int     NumExtraCalcs=9; // Number of extra calculations possible, e.g., OUT_PE

  // other general variables
  gzFile *fin;
  FILE *flist;

  char ***GridFileName;
  char filename[512];
  char GridName[512];
  char StatType;
  char TmpType[20];
  char TmpStr[20];
  char OVERWRITE;
  char PrtAllGrids = FALSE;
  char ArcGrid = FALSE;
  char ProcessFile;

  int filenum;
  int Nfile;
  int row;
  int nrows;
  int col;
  int ncols;
  int cidx, didx;
  int DONE, LAST;
  int Nrec;
  int ErrNum;
  int NumSets, SumSets, Nvals;
  int sidx, sumidx;
  int Ncells;
  int ValCnt;

  double maxlat;
  double minlat;
  double maxlng;
  double minlng;
  double tmplat;
  double tmplng;
  double findlat;
  double findlng;
  double cellsize;
  double **gridval;
  int NODATA = -9999;
  double CheckSum;

  double *GridLat;
  double *GridLng;
  double ****GridValues, **OutGrid;

  DATE_STRUCT startdate;
  DATE_STRUCT enddate;
  DATE_STRUCT nextdate, lastdate;
  DATE_STRUCT tmpdate;
  DATE_STRUCT periodstart, periodend;

  char  **RawData;
  char  **ColNames, *ColTypes, *ColAggTypes;
  char  **VarList;
  double  *ColMults;
  int     TimeStep, NumLayers, NumNodes, NumBands, NumFrostFronts;
  int     NumLakeNodes, NumCols, NumOut, OutStart;
  int     NumBytes, NumRead=0, ReadBytes, RawPtr;
  int    *date;
  double  *data;
  int     BinaryFile, lidx;
  int     periodyear[2];
  
  date      = (int *)calloc(NUM_DATE_VALS,sizeof(int));

  /** Process Command Line Arguments **/
  if ( argc != 11 ) {
    fprintf (stderr, "NumVars = %i\n", argc );
    fprintf(stderr,"\nUsage: %s <file list file> <Output ArcGrid: TRUE/FALSE> <output prefix> <grid resolution> <column list file> <start date> <end date> <V|A|S|M|W|D> <OVERWRITE: TRUE/FALSE> <PrtAllPeriods: TRUE/FALSE>\n",argv[0]);
    fprintf(stderr,"\n\tThis program produces either an ARC/INFO ASCII grid file \n\t(output ArcGrid = T) or an XYZ file (Output ArcGrid = F) of the average\n\tof the selected data column for the given time period.\n\n");
    fprintf(stderr,"\tThe following additional variables can be calculated by this \n\tprogram, and used with the standard list of summary metrics, the \n\trequired columns are in the VIC files being processed.  \n\tAdd the variable names below to your statistics control file \n\tto include them in the analysis:\n");
    fprintf(stderr,"\t- OUT_PE Penman potential evaporation, requires \"OUT_R_NET\", \n\t  \"OUT_GRND_FLUX\", \"OUT_WIND\", \"OUT_SURF_TEMP\", \"OUT_REL_HUMID\", \n\t  \"OUT_AIR_TEMP\",\n");
    fprintf(stderr,"\t- OUT_TOTAL_RUNOFF is the sum of \"OUT_RUNOFF\" and \"OUT_BASEFLOW\",\n");
    fprintf(stderr,"\t- OUT_TOTAL_SOIL_MOIST is the sum of \"OUT_SOIL_MOIST\" for each layer,\n");
    fprintf(stderr,"\t- OUT_MGDD modified growing degree day (from mrcc.isws.illinois.edu), \n\t  requires \"IN_TMIN\" and \"IN_TMAX\",\n");
    fprintf(stderr,"\t- OUT_PLANT_DAY is a Boolean variable indicating whether or not soil \n\t  surface conditions support planting, requires \"OUT_SOIL_MOIST_0\" \n\t  and \"OUT_SOIL_TEMP_0\".  Also requires that a file with top soil \n\t  layer field capacity [mm] be provided following the statistic \n\t  definition in the soil control file.\n");
    fprintf(stderr,"\t- OUT_CHILL_HR is the number of chilling hours (fruit trees), \n\t  requires \"IN_TMIN\" and \"IN_TMAX\".\n");
    fprintf(stderr,"\t- OUT_DMWD is the number of working days defined by DRAINMOD, \n\t  requires \"OUT_PREC\" and \"OUT_SOIL_MOIST_0\".  Also requires that \n\t  a file with top soil layer saturation [mm] be provided following \n\t  the statistic definition in the soil control file.\n");
    fprintf(stderr,"\t- OUT_DSFW is the number of of days suitable for field work \n\t  (Gramig et al, 2017), requires \"IN_PREC\", \"IN_TMIN\", \"IN_TMAX\", \n\t  and \"Soil Drainge Class\".\n");
    fprintf(stderr,"\t<file list> is a file containing the full grid file name and location, \n\t\tlatitude and longitude of the grid cell, for each grid cell to be\n\t\tincluded.\n");
    fprintf(stderr,"\t<output prefix> is the prefix (path and start of file name) for the \n\t\toutput files that will be generated by this program.  A suffix \n\t\twill be added to all file names to separate individual output \n\t\tfor each variable, and for each tperiod (year, season, month, \n\t\tetc). Files containing multi-year average statistics for annual, \n\t\tseasonal and monthly periods will use \"9999\" for \n\t\tthe date in the file name.\n");
    fprintf(stderr,"\t<grid resolution> is the resolution in degrees of the desired output\n\t\tgrid.\n");
    fprintf(stderr,"\t<column list file> is a multi-column ASCII file that lists the column\n\t\tname for each column to be output.  Followed by the statistic \n\t\tto be computed and a thresehold if required by the statistic.  \n\n\t\tStatistic options include: \n\t\t- Mean value \t\t\t\t\'mean\', \n\t\t- Cumulative value \t\t\t\'sum\', \n\t\t- Standard Deviation \t\t\t\'stdev\', \n\t\t- Maximum value \t\t\t\'max\', \n\t\t- Minimum value \t\t\t\'min\', \n\t\t- First value \t\t\t\t\'first\', \n\t\t- Last value \t\t\t\t\'last\', \n\t\t- Last day over thres.\t\t\t\'ldayo\' <thres>, \n\t\t- Last day under thres.\t\t\t\'ldayu\' <thres>, \n\t\t- First day over thres.\t\t\t\'fdayo\' <thres>, \n\t\t- First day under thres.\t\t\'fdayu\' <thres>, \n\t\t- Days over threshold\t\t\t\'othres\' <thres>,\n\t\t- Days under threshold\t\t\t\'uthres\' <thres>,\n\t\t- Average days over threshold\t\t\'avgdaysothres\' <thres>,\n\t\t- Average days under threshold\t\t\'avgdaysuthres\' <thres>,\n\t\t- Average value over threshold\t\t\'avgvalothres\' <thres>,\n\t\t- Average value under threshold\t\t\'avgvaluthres\' <thres>,\n\t\t- Consecutive days over threshold\t\'daysothres\' <thres>,\n\t\t- Consecutive days under threshold\t\'daysuthres\' <thres>,\n\t\t- Last day over thres before middle\t\'ldaymido\' <thres>,\n\t\t- Last day under thres before middle\t\'ldaymidu\' <thres>,\n\t\t- First day over thres after middle\t\'fdaymido\' <thres>,\n\t\t- First day under thres after middle\t\'fdaymidu\' <thres>,\n\t\t- Number of times threshold is crossed\t\'crossthres\' <thres>,\n\t\t- RB Index (flashiness)\t\t\t\'RBI\',\n\t\t- TQ mean (days spent above mean)\t\'Tqmean\',\n\t\t- Seven day low value\t\t\t\'7daylow\',\n\t\t- Quantile value\t\t\t\'quan\' <quantile>.\n\n");
    fprintf(stderr,"\t\tNOTE: Any single value threshold can be replaced with \n\t\t\"FILE <shortname> <filename>\" to provide spatially \n\t\tdistributed thresholds.  <shortname> is used to replace the \n\t\tthreshold value in the output filename, and <filename> must \n\t\tbe the same format (Arc Grid/XYZ) being used \n\t\tfor overall processing.\n");
    fprintf(stderr,"\t<start date> and <end date> are the starting and ending dates of the\n\t\tperiod of interest in MMDDYYYY format (MM = month, DD = day,\n\t\tYYYY = year - date must be 8 characters).\n");
    fprintf(stderr,"\t<V|A|S|M|W|D> export Annual a(V)erage, (A)nnual, (S)easonal, (M)onthly, \n\t\t(W)eekly, (D)aily or annual repeating (P)eriod grids for all \n\t\tyears in the file\n\t\t- Annual uses given start date to start year; \n\t\t- Seasonal uses Winter = DJF, Spring = MAM, Summer = JJA, and \n\t\t  Autumn = SON; \n\t\t- Weekly parses data into 7 day weeks starting with the given \n\t\t  start date.  \n\t\t- Period requires a range in MM-DD_MM-DD format or key phrase \n\t\t  \"grow\" for dynamic growing season, for example \n\t\t  \"P08-15_11-01\" to define the fall working window of \n\t\t  August 15th to November 1st.\n");
    fprintf(stderr,"\t<OVERWRITE> if set to TRUE then the output grid files will be \n\t\toverwritten.  The default is to replace values in existing \n\t\tfiles with new data.\n");
    fprintf(stderr,"\t<PrtAllPeriods> if set to TRUE then the program will write all output \n\t\tfiles (all years, all seasons, etc), if set to FALSE only \n\t\tannual averages will be written.\n\n");
    exit(0);
  }

  /********************************
    Read grid cell file list and create template for output files
  ********************************/
  if((flist=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open file list %s.\n",argv[1]);
    exit(0);
  }
  Nfile=0;
  // find extent of spatial coverage
  fscanf(flist,"%*s %lf %lf",&tmplat,&tmplng);
  maxlat = minlat = tmplat;
  maxlng = minlng = tmplng;
  while(!feof(flist)) {
    if(tmplat>maxlat) maxlat = tmplat;
    if(tmplat<minlat) minlat = tmplat;
    if(tmplng>maxlng) maxlng = tmplng;
    if(tmplng<minlng) minlng = tmplng;
    Nfile++;
    fscanf(flist,"%*s %lf %lf",&tmplat,&tmplng);
  }
  Nfile++;
  rewind(flist);
  // check if output file is in arcgrid format
  if ( strcasecmp( argv[2], "TRUE" ) == 0 ) ArcGrid = TRUE;
  else ArcGrid = FALSE;
  // check output file cellsize
  if ( ( cellsize = atof(argv[4]) ) <= 0 ) {
    fprintf(stderr, "ERROR: Cell size (%f) must be a positive doubling value!\n\n", cellsize);
    return(-1);
  }
  if ( ArcGrid ) {
    // compute latitude and longitude limits
    maxlat += cellsize / 2.;
    minlat -= cellsize / 2.;
    maxlng += cellsize / 2.;
    minlng -= cellsize / 2.;
    // estimate number of rows and columns
    nrows = (int)((maxlat - minlat) / cellsize);
    ncols = (int)((maxlng - minlng) / cellsize);
    Ncells = nrows*ncols;
    // check that input file list has spatial extent
    if ( nrows == 1 || ncols == 1 ) {
      fprintf( stderr, "WARNING: Defined cells are at best one-dimensional with %i rows and %i columns.  If you are expecting a raster output, then double check your file list format.\n", nrows, ncols );
    }
  }
  else {
    // Using XYZ file, so adjust nrows and ncols to reflect number of cells 
    // being processed
    nrows = Nfile;
    ncols = 1;
    Ncells = Nfile;
  }

  /*****************************
    Setup input and output directories
  *****************************/
  strcpy(GridName,argv[3]);

  /** Build temporary flux file name and open to read header information **/
  fscanf(flist,"%s %lf %lf",filename,&tmplat,&tmplng);
  fprintf(stdout,"Getting header information from %s...\n",filename);
  rewind(flist);
  if ( ( fin  = gzopen(filename, "rb") ) == NULL ) {
    // Try opening a compressed file of the same name
    strcat(filename,".gz");
    puts(filename);
    if ( ( fin  = gzopen(filename, "rb") ) == NULL ) {
      printf("File %s does not exist\n", filename);
      exit(1);
    }
  }
  /** Read output file header **/
  BinaryFile = get_header( &fin, &ColNames, &ColTypes, &ColMults, 
			   &ColAggTypes, &TimeStep, &NumLayers, &NumNodes, 
			   &NumBands, &NumFrostFronts, &NumLakeNodes, 
			   &NumCols, &NumBytes );
  ColNames = (char **)realloc(ColNames,(NumCols+NumExtraCalcs)*sizeof(char *));

  /***** Add Additional Columns for Calculation *****/

  // set OUT_PE as a column name
  ColNames[NumCols] = (char *)calloc( 7, sizeof(char) );
  strcpy( ColNames[NumCols], "OUT_PE" );
  // set OUT_TOTAL_RUNOFF as a column name
  ColNames[NumCols+1] = (char *)calloc( 17, sizeof(char) );
  strcpy( ColNames[NumCols+1], "OUT_TOTAL_RUNOFF" );
  // set OUT_TOTAL_SOIL_MOIST as a column name
  ColNames[NumCols+2] = (char *)calloc( 17, sizeof(char) );
  strcpy( ColNames[NumCols+2], "OUT_TOTAL_SOIL_MOIST" );
  // set OUT_MGDD as a column name
  ColNames[NumCols+3] = (char *)calloc( 9, sizeof(char) );
  strcpy( ColNames[NumCols+3], "OUT_MGDD" );
  // set OUT_PLANT_DAY as a column name
  ColNames[NumCols+4] = (char *)calloc( 13, sizeof(char) );
  strcpy( ColNames[NumCols+4], "OUT_PLANT_DAY" );
  // set OUT_CHILL_HR as a column name
  ColNames[NumCols+5] = (char *)calloc( 13, sizeof(char) );
  strcpy( ColNames[NumCols+5], "OUT_CHILL_HR" );
  // set OUT_BUD_FREEZE as a column name
  ColNames[NumCols+6] = (char *)calloc( 14, sizeof(char) );
  strcpy( ColNames[NumCols+6], "OUT_BUD_FREEZE" );
  // set OUT_DMWD as a column name
  ColNames[NumCols+7] = (char *)calloc( 8, sizeof(char) );
  strcpy( ColNames[NumCols+7], "OUT_DMWD" );
  // set OUT_DSFW as a column name
  ColNames[NumCols+8] = (char *)calloc( 8, sizeof(char) );
  strcpy( ColNames[NumCols+8], "OUT_DSFW" );

  /***** Check for VIC model output file time step *****/

  // check for hourly column in dataset
  if ( strcmp( ColNames[3], "HOUR" ) == 0 ) { // hour is present offset by one column
    OutStart = 4;
    VarList = &ColNames[OutStart];
    NumOut  = NumCols - OutStart + NumExtraCalcs;
  }
  else { // hour is not present, do not offset daily data
    OutStart = 3;
    VarList = &ColNames[OutStart];
    NumOut = NumCols - OutStart + NumExtraCalcs;
  }
  data = (double *)calloc((NumOut+2),sizeof(double));
  gzclose(fin);

  /****************************
    Work on determining statistics to be computed 
  ****************************/
  /// get statistic type
  strcpy(TmpType,argv[8]);
  StatType = TmpType[0];

  /** Read statistics control file **/
  ErrNum = GetStatInfo( argv[5], &StatInfo, VarList, NumOut, &CalcPE, 
			&CalcTR, &CalcTSM, &CalcMGDD, &CalcPD, &CalcCH, 
			&CalcBF, &CalcDMWD, &CalcDSFW, ArcGrid, GridLat, 
			GridLng, nrows, ncols, minlat, minlng, cellsize, 
			NODATA );

  /** Setup for calculation of PE **/
  if ( CalcPE )
    CalcPE = GetPenmanInfo( &PenInfo, VarList, NumOut );

  /** Setup for calculation of total runoff (runoff + baseflow) **/
  if ( CalcTR )
    CalcTR = GetTotalRunoffInfo( &TRoffInfo, VarList, NumOut );

  /** Setup for calculation of total soil moisture (all soil moisture layers) **/
  if ( CalcTSM )
    CalcTSM = GetTotalSoilMoistureInfo( &TSMInfo, VarList, NumOut, NumLayers );

  /** Setup for calculation of modified growing degree day (Tmin and Tmax) **/
  if ( CalcMGDD )
    if ( StatType == 'A' || StatType == 'P' )
      CalcMGDD = GetModifiedGrowingDegreeDayInfo( &MGDDInfo, &StatInfo, VarList, NumOut );
    else {
      fprintf( stderr, "MGDD not computed because period must be annual.\n" );
      CalcMGDD = 0;
    }

  /** Setup for calculation of first planting day (Soil Surface Temp and Moisture) **/
  if ( CalcPD )
    if ( StatType == 'A' || StatType == 'P' )
      CalcPD = GetPlantingDayInfo( &PDayInfo, &StatInfo, VarList, NumOut );
    else {
      fprintf( stderr, "Planting Day not computed because period must be annual.\n" );
      CalcPD = 0;
    }

  /** Setup for calculation of chilling hours using Tmin and Tmax **/
  if ( CalcCH )
    if ( StatType == 'A' || StatType == 'P' ) {
      CalcCH = GetChillingHoursInfo( &ChillHrInfo, &StatInfo, VarList, NumOut );
    }
    else {
      fprintf( stderr, "Chilling Hours not computed because period must be annual.\n" );
      CalcCH = 0;
    }

  /** Setup for calculation of bud freeze potential **/
  if ( CalcBF )
    if ( StatType == 'A' || StatType == 'P' ) {
      CalcCH = GetBudFreezeInfo( &BudFreezeInfo, &StatInfo, VarList, NumOut );
    }
    else {
      fprintf( stderr, "Bud Freeze potential not computed because period must be annual.\n" );
      CalcCH = 0;
    }

  /** Setup for calculation of DRAINMOD working days using PREC and surface MOIST **/
  if ( CalcDMWD )
    CalcDMWD = GetDrainModWorkingDaysInfo( &DModWDaysInfo, &StatInfo, VarList, NumOut );
  else CalcDMWD = FALSE;

  /** Setup for calculation of DSFW using PREC, TMIN, TMAX and drainage class **/
  if ( CalcDSFW )
    CalcDSFW = GetDaysSuitable4FieldWorkInfo( &DSFWInfo, &StatInfo, VarList, NumOut );
  else CalcDSFW = FALSE;

  /** Check that requested columns are in the given output files **/
  for ( didx=NumOut-1; didx>=0; didx-- ) {
    for( cidx=0; cidx<NumCols+NumExtraCalcs; cidx++ ) {
      if ( strcasecmp( ColNames[cidx], VarList[didx] ) == 0 ) break;
    }
    if ( cidx >= NumCols+NumExtraCalcs ) {
      fprintf( stderr, "WARNING: Unable to find output column %s in the given output files, will remove from analysis, check the VIC model output file header to determine what variables are in the file %s.\n", VarList[didx], filename );
      for( cidx=didx; cidx<NumOut-1; cidx++ )
	VarList[cidx] = VarList[cidx+1]; // remove missing statistic
      NumOut--;
    }
  }

  /**********************************
     Initialize date structure
  **********************************/

  // get start date
  startdate.month   = (int)(atof(argv[6])/1000000);
  startdate.day     = (int)(atof(argv[6])/10000) - startdate.month*100;
  startdate.year    = (int)(atof(argv[6])) - startdate.month*1000000 - startdate.day*10000;
  startdate.hour    = 0;
  startdate.juldate = calc_juldate(startdate.year,startdate.month,startdate.day,0.);
  startdate         = get_season(startdate);
  /// get end date
  enddate.month     = (int)(atof(argv[7])/1000000);
  enddate.day       = (int)(atof(argv[7])/10000) - enddate.month*100;
  enddate.year      = (int)(atof(argv[7])) - enddate.month*1000000 - enddate.day*10000;
  enddate.hour      = 0;
  enddate.juldate   = calc_juldate(enddate.year,enddate.month,enddate.day,0);
  enddate           = get_season(enddate);
  
  // Process statistic type - adjust start and end date to encompass full periods
  if ( StatType == 'A' ) {
    // Annual statistics
    if ( enddate.month == startdate.month && enddate.day == startdate.day )
      NumSets = enddate.year - startdate.year;
    else
      NumSets = enddate.year - startdate.year + 1;
    SumSets = 1;
    strcat( GridName, "_\%i_\%s.asc" ); // build annual output file name
  }
  else if ( StatType == 'S' ) {
    // Seasonal statistics
    if ( startdate.month > 2 && startdate.month < 6 )
      // Start in Spring
      startdate.juldate = calc_juldate( startdate.year, 3, 1, 0 );
    else if ( startdate.month > 5 && startdate.month < 9 )
      // Start in Summer
      startdate.juldate = calc_juldate( startdate.year, 6, 1, 0 );
    else if ( startdate.month > 5 && startdate.month < 9 )
      // Start in Autumn
      startdate.juldate = calc_juldate( startdate.year, 9, 1, 0 );
    else 
      // Start in Winter
      startdate.juldate = calc_juldate( startdate.year, 12, 1, 0 );
    if ( enddate.month > 2 && enddate.month < 6 )
      // End in Winter
      enddate.juldate = calc_juldate( enddate.year, 3, 1, 0 );
    else if ( enddate.month > 5 && enddate.month < 9 )
      // End in Spring
      enddate.juldate = calc_juldate( enddate.year, 6, 1, 0 );
    else if ( enddate.month > 5 && enddate.month < 9 )
      // End in Summer
      enddate.juldate = calc_juldate( enddate.year, 9, 1, 0 );
    else if ( enddate.month >= 1 && enddate.month < 2 )
      // End in Winter
      enddate.juldate = calc_juldate( enddate.year-1, 12, 1, 0 );
    else 
      // End in Winter
      enddate.juldate = calc_juldate( enddate.year, 12, 1, 0 );
    // get number of seasons to process
    sidx = 0;
    nextdate = startdate;
    while ( nextdate.juldate < enddate.juldate ) {
      nextdate = get_next_season( nextdate );
      sidx++;
    }
    NumSets = sidx;
    SumSets = 4;
    strcat( GridName, "_\%i_\%s_\%s.asc" ); // build seasonal output file name
  }
  else if ( StatType == 'M' ) {
    // Monthly Statistics
    startdate.juldate = calc_juldate( startdate.year, startdate.month, 1, 0 );
    enddate.juldate = calc_juldate( enddate.year, enddate.month, 1, 0 );
    if ( enddate.month == startdate.month && enddate.day == startdate.day )
      NumSets = (enddate.year - startdate.year) * 12;
    else
      NumSets = (enddate.year - startdate.year + 1) * 12;
    SumSets = 12;
    strcat( GridName, "_\%i_\%s_\%s.asc" ); // build monthly output file name
  }
  else if ( StatType == 'W' ) {
    // Weekly Statistics
    NumSets = (int)((enddate.juldate - startdate.juldate)/7.);
    SumSets = 0;
    strcat( GridName, "_\%i_\%i_\%i_\%s.asc" ); // build weekly output file name
  }
  else if ( StatType == 'V' ) {
    // Overall Average Statistics
    NumSets = 1;
    SumSets = 0;
    strcat( GridName, "_\%s.asc" ); // build overall average output file name
  }
  else if ( StatType == 'D' ) {
    // Daily Statistics
    NumSets = (int)(enddate.juldate - startdate.juldate);
    SumSets = 0;
    strcat( GridName, "_\%i_\%u_\%s.asc" ); // build daily output file name
  }
  else if ( StatType == 'P' ) {
    // Annual period statistics
    //** Starting with fixed periods, then will try flexible **
    if ( enddate.month == startdate.month && enddate.day == startdate.day )
      NumSets = enddate.year - startdate.year;
    else
      NumSets = enddate.year - startdate.year + 1;
    SumSets = 1;
    // get start and end dates for annual period
    sscanf( TmpType+1, "%2s", &TmpStr ); 
    periodstart.month = atoi(TmpStr);
    sscanf( TmpType+4, "%2s", &TmpStr ); 
    periodstart.day = atoi(TmpStr);
    sscanf( TmpType+7, "%2s", &TmpStr ); 
    periodend.month = atoi(TmpStr);
    sscanf( TmpType+10, "%2s", &TmpStr ); 
    periodend.day = atoi(TmpStr);

    // set period start date
    tmpdate.year = startdate.year;
    tmpdate.month = startdate.month;
    tmpdate.day = periodstart.day;
    tmpdate.hour = 0;
    tmpdate = get_juldate(tmpdate);
    while ( periodstart.month != tmpdate.month ) 
      tmpdate = get_next_month( tmpdate );
    periodstart.year = tmpdate.year;
    periodstart.hour = 0;
    periodstart = get_juldate(periodstart);
    
    // set period end date
    tmpdate.day = periodend.day;
    tmpdate = get_juldate(tmpdate);
    while ( periodend.month != tmpdate.month ) 
      tmpdate = get_next_month( tmpdate );
    periodend.year = tmpdate.year;
    periodend.hour = 0;
    periodend = get_juldate(periodend);

    // check that length of period is less than one year
    if ( ( periodend.juldate - periodstart.juldate ) > 366 ) {
      fprintf( stderr, "ERROR: Period defined is longer than a year, that option is not currently supported.  Program will exit with an error.\n" );
      return(-1);
    }

    // store starting period years
    periodyear[0] = periodstart.year;
    periodyear[1] = periodend.year;

    // build output file name
    strcat( GridName, "_" ); // change P to _ for file name
    strcat( GridName, TmpType ); // add period descriptor
    strcat( GridName, "_\%i_\%s.asc" ); // build annual output file name

  }
  else {
    fprintf( stderr, "ERROR: Unknown statistics period variable (%s) program will exit.\n", StatType );
    return (-1);
  }

  fprintf( stderr, "GridName for output files = %s\n", GridName );

  // Get file overwrite flag
  if ( strcasecmp( argv[9], "TRUE" ) == 0 ) OVERWRITE = TRUE;
  else OVERWRITE = FALSE;

  // Get write all periods flag
  if ( strcasecmp( argv[10], "TRUE" ) == 0 ) PrtAllGrids = TRUE;
  else PrtAllGrids = FALSE;

  // Correct PrtAllGrids for situation where statistics are being calculated on less than a full year of data
  if ( SumSets == 0 ) PrtAllGrids = TRUE;

  /**********************************
    Initialize Data Array
  ***********************************/
  gridval = (double **)calloc(StatInfo.Ncols,sizeof(double *));  
  for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
    gridval[cidx] = (double *)calloc(366,sizeof(double));  
  }
  GridFileName = (char ***)calloc((NumSets+SumSets),sizeof(char **));
  for ( sidx = 0; sidx < NumSets+SumSets; sidx++ ) {
    GridFileName[sidx] = (char **)calloc(StatInfo.Ncols,sizeof(char *));
    for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) 
      GridFileName[sidx][cidx] = (char *)calloc(1024,sizeof(char));
  }

  /***************************************************
    Define variables required for processing
  ***************************************************/
  if ( !ArcGrid ) {
    // define latitude and longitude arrays for storing positions.
    //  This is only needed for XYZ files, in which case nrows is set 
    //  to the number of files, while ncols is set to 1.
    GridLat = (double *)calloc(nrows, sizeof(double));
    GridLng = (double *)calloc(nrows, sizeof(double));
    for ( row = 0; row < nrows; row++ ) 
      fscanf(flist,"%*s %lf %lf",&GridLat[row],&GridLng[row]);
    rewind(flist);
  }
 
  // define grid value array for storing all data
  GridValues = (double ****)calloc(NumSets,sizeof(double ***)); 
  for( sidx =0; sidx<NumSets; sidx++) {
    GridValues[sidx] = (double ***)calloc(StatInfo.Ncols,sizeof(double **));
    for( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
      GridValues[sidx][cidx] = (double **)calloc(nrows,sizeof(double *));
      OutGrid = (double **)calloc(nrows,sizeof(double *));
      for ( row = 0; row < nrows; row++ ) {
	GridValues[sidx][cidx][row] = (double *)calloc(ncols,sizeof(double));
	OutGrid[row] = (double *)calloc(ncols,sizeof(double));
	for ( col = 0; col < ncols; col++ )
	  GridValues[sidx][cidx][row][col] = NODATA;
      }
    }
  }

  // establish next date as equal to start date for now
  nextdate.day     = startdate.day;
  nextdate.month   = startdate.month;
  nextdate.year    = startdate.year;
  nextdate.hour    = startdate.hour;
  nextdate.juldate = startdate.juldate;
  nextdate         = get_season(nextdate);

  /******************************************
    Build blank files for overall statistics
  ******************************************/
  for ( sidx = 0; sidx < SumSets; sidx++ ) {
    for( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
      if ( strcasecmp( StatInfo.ColStatList[cidx], "none" ) != 0 ) {
  	// Going to output statistics as an ArcInfo grid
  	ErrNum = CheckExistingOutputFiles(ArcGrid, StatType, GridName,
					  StatInfo, GridFileName,
					  nextdate, NumSets+sidx, cidx,
					  NumSets, nrows, ncols,
					  minlat, minlng,
					  cellsize, NODATA, Ncells, 
					  OVERWRITE);
      }
    }
    if ( StatType == 'A' ) nextdate = get_next_year( nextdate );
    else if ( StatType == 'S' ) nextdate = get_next_season( nextdate );
    else if ( StatType == 'M' ) nextdate = get_next_month( nextdate );
    else if ( StatType == 'W' ) nextdate = get_next_week( nextdate );
    else if ( StatType == 'V' ) nextdate = enddate;
    else if ( StatType == 'P' ) nextdate = get_next_year( nextdate );
    else nextdate = get_next_day( nextdate );
    nextdate = get_juldate( nextdate );
  }

  /*******************/
  /** Process Files **/
  /*******************/

  for ( filenum = 0; filenum < Nfile; filenum++ ) {

    ProcessFile = TRUE;

    /* Build flux file name and open */
    fscanf(flist,"%s %lf %lf",filename,&tmplat,&tmplng);
    fprintf(stdout,"%i: processing %s...\n",filenum,filename);
    if ( ( fin  = gzopen(filename, "rb") ) == NULL ) {
      // Try opening a compressed file of the same name
      strcat(filename,".gz");
      puts(filename);
      if ( ( fin  = gzopen(filename, "rb") ) == NULL ) {
	printf("File %s does not exist\n", filename);
	exit(1);
      }
    }
    // read file header
    BinaryFile = get_header( &fin, &ColNames, &ColTypes, &ColMults, 
			     &ColAggTypes, &TimeStep, &NumLayers, &NumNodes, 
			     &NumBands, &NumFrostFronts, &NumLakeNodes, 
			     &NumCols, &NumBytes );
    DONE = FALSE;
      
    /**** Need to determine row and column of current file in grid file ****/
    if ( ArcGrid ) {
      // identify row and column of current input file
      row = nrows - (int)( ( tmplat - minlat ) / cellsize ) - 1;
      col = (int)( ( tmplng - minlng ) / cellsize );

      // Check that current row and column are in the Arc/Info file boundaries
      if ( row >= nrows || col >= ncols ) {
	fprintf(stderr,"WARNING: %s located at (%i,%i) is outside the range of the grid file (%i,%i). Will not be processed.\n",
		filename,row,col,nrows,ncols);
	ProcessFile = FALSE;
      }
      fprintf(stdout,"\trow = %i,\tcol = %i\n",row,col);
    }
    else {
      // for xyz data output files, col will alwyas be 0, row equal to the filenum
      col = 0;
      row = filenum;
      fprintf(stdout,"\tline = %i\n",row);
    }

    if ( ProcessFile ) {

      // Initialize accounting variables
      DONE = FALSE;
      LAST = FALSE;
      sidx = NODATA;
      if ( StatType == 'P' ) {
 	periodstart.year = periodyear[0];
	periodstart = get_juldate( periodstart );
	periodend.year = periodyear[0];
	periodend = get_juldate( periodend );
      }
      if ( filenum == 0 ) {
	// Use first file to setup all additional processing
	// determine number of bytes in full flux file, in case user has requested shorter time period
	if ( BinaryFile ) {
	  // determine number of bytes in full binary flux file, in case user has requested shorter time period
	  RawData = (char **)calloc(1,sizeof(char*));
	  RawData[0] = (char *)calloc(NumBytes,sizeof(char));
	  while ( ( ReadBytes = gzread(fin,RawData[0],NumBytes*sizeof(char)) ) > 0 ) {
	    NumRead += ReadBytes;
	  }
	}
	else {
	  // determine number of lines in full ASCII flux file, in case shorter time period requested
	  RawData = (char **)calloc(1,sizeof(char *));
	  RawData[0] = (char *)calloc(MaxCharData,sizeof(char));
	  while ( ( gzgets( fin, RawData[0], MaxCharData) ) != NULL ) {
	    NumRead ++;  // in this case a line counter
	  }
	}
	// reset start of flux file
	gzrewind( fin );
	BinaryFile = get_header( &fin, &ColNames, &ColTypes, &ColMults, 
				 &ColAggTypes, &TimeStep, &NumLayers, 
				 &NumNodes, &NumBands, &NumFrostFronts, 
				 &NumLakeNodes, &NumCols, &NumBytes );
	// reallocate Raw data for full size of flux file
	if ( BinaryFile ) {
	  // allocation for binary file
	  if ( ( RawData[0] = (char *)realloc( RawData[0], NumRead ) ) == NULL ) {
	    fprintf( stderr, "ERROR: Unable to allocate sufficient memory to read in full flux file.\n" );
	    return (-1);
	  }
	}
	else {
	  // free previous allocation for ASCII file
	  free ((char*)RawData[0]);
	  free ((char*)RawData);
	  // allocation for ASCII file
	  if ( ( RawData = (char **)calloc( NumRead,sizeof(char * ) ) ) == NULL ) {
	    fprintf( stderr, "ERROR: Unable to allocate sufficient memory to read in full flux file.\n" );
	    return (-1);
	  }
	  for ( lidx = 0; lidx < NumRead; lidx ++ ) {
	    if ( ( RawData[lidx] = (char *)calloc( MaxCharData,sizeof(char) ) ) == NULL ) {
	      fprintf( stderr, "ERROR: Unable to allocate sufficient memory to read in full flux file.\n" );
	      return (-1);
	    }
	  }
	}
      }

      // Read in full flux file
      if ( BinaryFile ) {
	// read in entire binary data file
	if ( ( ReadBytes = gzread(fin,RawData[0],NumRead) ) != NumRead ) {
	  fprintf( stderr, "WARNING: Unable to read in complete flux file %s, got %d of %ld bytes.\n", filename, ReadBytes, NumRead*sizeof(char) );
	  //return (-1);
	}
      }
      else {
	// read in entire ASCII data file
	for ( lidx = 0; lidx < NumRead; lidx ++ ) {
	  if ( gzgets(fin,RawData[lidx],MaxCharData) == NULL ) {
	    fprintf( stderr, "WARNING: Unable to read in complete flux file %s, got only %d of %d lines.\n", filename, lidx, NumRead );
	    //return (-1);
	  }
	}
      }

      // Process file until complete
      RawPtr = 0;
      
      // read date from current line
      ErrNum = get_record_PEN( BinaryFile, RawData, &RawPtr, ColNames, ColTypes, 
			       ColMults, date, data, NumCols, NumOut, &OutStart, 
			       CalcPE, PenInfo, CalcTR, TRoffInfo, CalcTSM, TSMInfo );
      
      // reset to correct ASCII position before processing
      if ( !BinaryFile ) RawPtr = 0; 
      
      // compute julian day of current record
      tmpdate.year = date[0];
      tmpdate.month = date[1];
      tmpdate.day = date[2];
      tmpdate.hour = 0;
      tmpdate = get_juldate(tmpdate);
      
      while( !DONE && ErrNum >= 0 ) { 
	// process all lines in the current file
	
	if ( tmpdate.juldate >= startdate.juldate && sidx == NODATA ) {
	  // Initialize date accounting for new grid cell
	  sidx = 0;
	  lastdate = copy_date( startdate );
	  if ( StatType == 'A' ) nextdate = get_next_year( startdate );
	  else if ( StatType == 'S' ) nextdate = get_next_season( startdate );
	  else if ( StatType == 'M' ) nextdate = get_next_month( startdate );
	  else if ( StatType == 'W' ) nextdate = get_next_week( startdate );
	  else if ( StatType == 'V' ) nextdate = enddate;
	  else if ( StatType == 'P' ) nextdate = get_next_year( startdate );
	  else nextdate = get_next_day( startdate );
	  nextdate = get_juldate( nextdate );
	  // Create output file names and initialize grid files
	  for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
	    if ( strcasecmp( StatInfo.ColStatList[cidx], "none" ) != 0 ) {
	      if ( filenum == 0 && PrtAllGrids ) {
	  	ErrNum = CheckExistingOutputFiles(ArcGrid, StatType, GridName,
						  StatInfo, GridFileName,
						  startdate, sidx, cidx,
						  NumSets, nrows, ncols,
						  minlat, minlng,
						  cellsize, NODATA, Ncells,
						  OVERWRITE);
	      }
	    }
	    // Initialize daily statistic storage variables for each new grid cell
	    for ( didx = 0; didx < 366; didx++ ) {
	      gridval[cidx][didx] = (double)NODATA;  
	    }
	  }
	  Nrec = 0;
	}
	else if ( tmpdate.juldate >= nextdate.juldate && sidx >= 0 ) {
	  // process next set (annual, seasonal, monthly, etc) of data
	  Nvals = (int)(nextdate.juldate+0.5)-(int)(lastdate.juldate+0.5);
	  // calculate modified growing degree days, if possible
	  if ( CalcMGDD ) {
	    CalcModifiedGrowingDegreeDays( gridval[MGDDInfo.ColNumList[2]], 
					   gridval[MGDDInfo.ColNumList[0]], 
					   gridval[MGDDInfo.ColNumList[1]], 
					   lastdate, Nvals );
	    for ( cidx = 1; cidx < MGDDInfo.Noutput; cidx++ )
	      // extra values, so copy statistic
	      for ( didx = 0; didx < Nvals; didx++ )
		gridval[MGDDInfo.OutputCols[cidx]][didx] 
		  = gridval[MGDDInfo.ColNumList[2]][didx];
	  }
	  // calculate first day of planting, if possible
	  if ( CalcPD ) {
	    CalcPlantingDay( gridval[PDayInfo.ColNumList[2]], 
			     gridval[PDayInfo.ColNumList[0]], 
			     gridval[PDayInfo.ColNumList[1]], 
			     StatInfo.Extra[PDayInfo.OutputCols[0]][row][col], 
			     lastdate, Nvals );
	    for ( cidx = 1; cidx < PDayInfo.Noutput; cidx++ )
	      // extra values, so copy statistic
	      for ( didx = 0; didx < Nvals; didx++ )
		gridval[PDayInfo.OutputCols[PDayInfo.OutputCols[0]]][didx] 
		  = gridval[PDayInfo.ColNumList[2]][didx];
	  }
	  // calculate chilling hours, if possible
	  if ( CalcCH ) {
	    CalcChillingHours( gridval[ChillHrInfo.ColNumList[2]], 
			       gridval[ChillHrInfo.ColNumList[0]], 
			       gridval[ChillHrInfo.ColNumList[1]], 
			       lastdate, Nvals );
	    for ( cidx = 1; cidx < ChillHrInfo.Noutput; cidx++ )
	      // extra values, so copy statistic
	      for ( didx = 0; didx < Nvals; didx++ )
		gridval[ChillHrInfo.OutputCols[cidx]][didx] 
		  = gridval[ChillHrInfo.ColNumList[2]][didx];
	  }
	  // calculate bud freeze potential, if possible
	  if ( CalcBF ) {
	    CalcBudFreeze( gridval[BudFreezeInfo.ColNumList[2]], 
			   gridval[BudFreezeInfo.ColNumList[0]], 
			   gridval[BudFreezeInfo.ColNumList[1]], 
			   lastdate, Nvals );
	    for ( cidx = 1; cidx < BudFreezeInfo.Noutput; cidx++ )
	      // extra values, so copy statistic
	      for ( didx = 0; didx < Nvals; didx++ )
		gridval[BudFreezeInfo.OutputCols[cidx]][didx] 
		  = gridval[BudFreezeInfo.ColNumList[2]][didx];
	  }
	  // calculate DRAINMOD working days, if possible
	  if ( CalcDMWD ) {
	    CalcDrainModWorkingDays( gridval[DModWDaysInfo.ColNumList[2]], 
				     gridval[DModWDaysInfo.ColNumList[0]], 
				     gridval[DModWDaysInfo.ColNumList[1]], 
				     StatInfo.Extra[DModWDaysInfo.OutputCols[0]][row][col], 
				     lastdate, Nvals );
	    for ( cidx = 1; cidx < DModWDaysInfo.Noutput; cidx++ )
	      // extra values, so copy statistic
	      for ( didx = 0; didx < Nvals; didx++ )
		gridval[DModWDaysInfo.OutputCols[DModWDaysInfo.OutputCols[0]]][didx] 
		  = gridval[DModWDaysInfo.ColNumList[2]][didx];
	  }
	  // calculate Days Suitable for Field Work (DSFW), if possible
	  if ( CalcDSFW ) {
	    CalcDaysSuitable4FieldWork( gridval[DSFWInfo.ColNumList[3]], 
					gridval[DSFWInfo.ColNumList[0]], 
					gridval[DSFWInfo.ColNumList[1]], 
					gridval[DSFWInfo.ColNumList[2]], 
					lastdate, Nvals );
	    for ( cidx = 1; cidx < DSFWInfo.Noutput; cidx++ )
	      // extra values, so copy statistic
	      for ( didx = 0; didx < Nvals; didx++ )
		gridval[DSFWInfo.OutputCols[cidx]][didx] 
		  = gridval[DSFWInfo.ColNumList[2]][didx];
	  }
	  // now process all requested output statistics
	  for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
	    // Check is processing annual periods
	    if ( StatType == 'P' ) {
	      // Set values outside the period to NoData
	      for ( didx = 0; didx < 366; didx++ ) {
		if ( lastdate.juldate + didx < periodstart.juldate ||
		     lastdate.juldate + didx > periodend.juldate )
		  gridval[cidx][didx] = (double)NODATA;
	      }
	    }
	    // compute base statistics
	    if ( strcasecmp(StatInfo.ColStatList[cidx], "max") == 0 )
	      // Find the maximum value
	      GridValues[sidx][cidx][row][col] = get_max(gridval[cidx], Nvals, (double)NODATA );
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "min") == 0 ) {
	      // Find the minimum value
	      GridValues[sidx][cidx][row][col] = get_min(gridval[cidx], Nvals, (double)NODATA );
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "othres", 6 ) == 0 ) {
	      // Count all days a value is under the given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_count_over_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "uthres", 6 ) == 0 ) {
	      // Count all days a value is under the given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_count_under_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "avgdaysothres", 13 ) == 0 ) {
	      // Count all days a value is under the given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_average_days_over_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "avgdaysuthres", 13 ) == 0 ) {
	      // Count all days a value is under the given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_average_days_under_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "avgvalothres", 12 ) == 0 ) {
	      // average value of events over the given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_average_value_over_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "avgvaluthres", 12 ) == 0 ) {
	      // average value of events under the given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_average_value_under_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "daysothres", 10 ) == 0 ) {
	      // Count all days a value is under the given threshold between first and last occurance centerd on the middle of the record
	      GridValues[sidx][cidx][row][col] = (double)get_consecutive_days_over_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "daysuthres", 10 ) == 0 ) {
	      // Count consecutive days a value is under the given threshold between first and last occurance centerd on the middle of the record
	      GridValues[sidx][cidx][row][col] = (double)get_consecutive_days_under_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "ldayo", 5 ) == 0 ) {
	      // find the last day a value is over a given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_last_day_over_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA, lastdate);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "ldayu", 5 ) == 0 ) {
	      // find the last day a value is under a given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_last_day_under_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA, lastdate);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "fdayo", 5 ) == 0 ) {
	      // find the first day a value is over a given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_first_day_over_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA, lastdate);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "fdayu", 5 ) == 0 ) {
	      // find the first day a value is under a given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_first_day_under_thres(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA, lastdate);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "ldaymido", 8 ) == 0 ) {
	      // find the last day a value is over a given threshold from middle record
	      GridValues[sidx][cidx][row][col] = (double)get_last_day_over_thres_from_middle(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA, lastdate);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "ldaymidu", 8 ) == 0 ) {
	      // find the last day a value is under a given threshold from middle record
	      GridValues[sidx][cidx][row][col] = (double)get_last_day_under_thres_from_middle(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA, lastdate);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "fdaymido", 8 ) == 0 ) {
	      // find the first day a value is over a given threshold from middle record
	      GridValues[sidx][cidx][row][col] = (double)get_first_day_over_thres_from_middle(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA, lastdate);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "fdaymidu", 8 ) == 0 ) {
	      // find the first day a value is under a given threshold from middle record
	      GridValues[sidx][cidx][row][col] = (double)get_first_day_under_thres_from_middle(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA, lastdate);
	    }
	    else if ( strncasecmp( StatInfo.ColStatList[cidx], "crossthres", 10 ) == 0 ) {
	      // count the number of times a threshold is crossed
	      GridValues[sidx][cidx][row][col] = (double)get_count_times_thres_crossed(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA);
	    }
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "stdev") == 0 ) {
	      // Find the standard deviation
	      GridValues[sidx][cidx][row][col] = get_stdev(gridval[cidx], get_mean(gridval[cidx], Nvals, (double)NODATA), Nvals, (double)NODATA );
	    }
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "sum") == 0 ) {
	      // Find the cumulative sum
	      GridValues[sidx][cidx][row][col] = get_sum(gridval[cidx], Nvals, (double)NODATA );
	    }
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "first") == 0 ) {
	      // Find the first value of the period
	      GridValues[sidx][cidx][row][col] = gridval[cidx][0];
	    }
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "last") == 0 ) {
	      // Find the last value of the period
	      GridValues[sidx][cidx][row][col] = gridval[cidx][Nvals-1];
	    }
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "RBI") == 0 ) {
	      // Find the cumulative sum
	      GridValues[sidx][cidx][row][col] = get_RB_index(gridval[cidx], Nvals, (double)NODATA );
	    }
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "Tqmean") == 0 ) {
	      // Find the cumulative sum
	      GridValues[sidx][cidx][row][col] = get_TQ_mean(gridval[cidx], Nvals, (double)NODATA );
	    }
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "7daylow") == 0 ) {
	      // Find the cumulative sum
	      GridValues[sidx][cidx][row][col] = get_seven_day_low_value(gridval[cidx], Nvals, (double)NODATA );
	    }
	    else if ( strncasecmp(StatInfo.ColStatList[cidx], "quan", 4) == 0 ) {
	      // Find the cumulative sum
	      GridValues[sidx][cidx][row][col] = get_quantile(gridval[cidx], StatInfo.Thres[cidx][row][col], Nvals, (double)NODATA );
	    }
	    else
	    // Find the mean
	      GridValues[sidx][cidx][row][col] = get_mean(gridval[cidx], Nvals, (double)NODATA );
	    // Create and initialize output file name for next increment
	    if ( filenum == 0 && tmpdate.juldate < enddate.juldate && PrtAllGrids ) {
	      if ( strcasecmp( StatInfo.ColStatList[cidx], "none" ) != 0 ) {
	    	ErrNum = CheckExistingOutputFiles(ArcGrid, StatType, GridName,
						  StatInfo, GridFileName,
						  nextdate, sidx+1, cidx,
						  NumSets, nrows, ncols,
						  minlat, minlng,
						  cellsize, NODATA, Ncells, 
						  OVERWRITE);
	      }
	    }
	    // Initialize daily statistic storage variables for each new grid cell
	    for ( didx = 0; didx < 366; didx++ ) {
	      gridval[cidx][didx] = (double)NODATA;  
	    }
	  }
	  Nrec = 0;
	  // Update Date Information
	  sidx++;
	  lastdate = copy_date( nextdate );
	  if ( StatType == 'A' ) nextdate = get_next_year( nextdate );
	  else if ( StatType == 'S' ) nextdate = get_next_season( nextdate );
	  else if ( StatType == 'M' ) nextdate = get_next_month( nextdate );
	  else if ( StatType == 'W' ) nextdate = get_next_week( nextdate );
	  else if ( StatType == 'P' ) {
	    nextdate = get_next_year( nextdate );
	    periodstart = get_next_year( periodstart );
	    periodend = get_next_year( periodend );
	  }
	  else nextdate = get_next_day( nextdate );
	  nextdate = get_juldate( nextdate );
	}
	
	// process line if date is between start and next dates
	if ( tmpdate.juldate >= startdate.juldate 
	     && tmpdate.juldate < nextdate.juldate ) {
	  for ( cidx = 0; cidx < StatInfo.Ncols; cidx ++ ) {
	    // store current value
	    gridval[cidx][(int)(tmpdate.juldate+0.5)-(int)(lastdate.juldate+0.5)] = data[StatInfo.ColNumList[cidx]];
	  }
	  Nrec++;
	}
	
	// read data from current line if there is still data left to read
	if ( LAST ) DONE = TRUE; // don't read potentially empty line
	else if ( !BinaryFile && RawPtr == NumRead ) {
	  // reached end of file, stop reading data but one more pass to process
	  LAST = TRUE;
	  tmpdate = get_next_day( tmpdate );
	}
	else if ( tmpdate.juldate < enddate.juldate ) {
	  // Still data to read so process the next line of data
	  ErrNum = get_record_PEN( BinaryFile, RawData, &RawPtr, ColNames, 
				   ColTypes, ColMults, date, data, NumCols, 
				   NumOut, &OutStart, CalcPE, PenInfo, CalcTR, 
				   TRoffInfo, CalcTSM, TSMInfo );
	  
	  if ( ErrNum >= 0 ) {
	    // compute julian day of current record
	    tmpdate.year = date[0];
	    tmpdate.month = date[1];
	    tmpdate.day = date[2];
	    tmpdate.hour = 0;
	    tmpdate = get_juldate(tmpdate);
	  }
	  else if ( tmpdate.juldate+1 == enddate.juldate ) {
	    tmpdate.year = enddate.year;
	    tmpdate.month = enddate.month;
	    tmpdate.day = enddate.day;
	    tmpdate.hour = 0;
	    tmpdate = get_juldate(tmpdate);
	    ErrNum = 0;
	  }	

	  // Check that current cell had valid data if processed
	  if ( sidx >= 0 ) {
	    CheckSum = 0;
	    for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ )
	      if ( gridval[cidx][0] == (double)NODATA ) CheckSum ++;
	    if ( CheckSum > 0 ) 
	      fprintf(stderr,"WARNING: Possible problem with cell at %f, %f, since %i summary statistics are NODATA.\n", tmplat, tmplng, CheckSum);
	  }

	}
	else { 
	  DONE = TRUE; 
	}
	
      }

    }
    
    /* Close grid file */
    gzclose(fin);
    
  }

  //
  // Output final statistics
  //

  for ( cidx = 0; cidx < StatInfo.Ncols; cidx ++ ) {

    if ( strcasecmp( StatInfo.ColStatList[cidx], "none" ) != 0 ) {

      // Plot individual time period grids
      if ( PrtAllGrids ) 
	for ( sidx = 0; sidx < NumSets; sidx++ ) // GRIDFILENAME should be formatted for the next date, so that that is not required here.
	  ErrNum = WriteOutputFiles( ArcGrid, GridFileName[sidx][cidx], 
				     GridValues[sidx][cidx], GridLat,
				     GridLng, nrows, ncols, minlat, 
				     minlng, cellsize, NODATA, Ncells,
				     OVERWRITE);

      // Compute average for all time periods
      for ( sumidx = 0; sumidx < SumSets; sumidx++ ) {
	for ( row = 0; row < nrows; row++ ) {
	  for ( col = 0; col < ncols; col++ ) {
	    // compute average value
	    OutGrid[row][col] = ValCnt = 0;
	    for ( sidx = sumidx; sidx < NumSets; sidx+=SumSets ) {
	      if ( fabs( GridValues[sidx][cidx][row][col] - NODATA ) > SMALL_CHK_VAL ) {
		OutGrid[row][col] += GridValues[sidx][cidx][row][col];
		ValCnt++;
	      }
	    }
	    if ( ValCnt == 0 )
	      // no values present, set to no data value
	      OutGrid[row][col] = NODATA;
	    else
	      // otherwise store the average value
	      OutGrid[row][col] /= (double)ValCnt;
	  }
	}

	// Update Grid file for current value
	ErrNum = WriteOutputFiles( ArcGrid, GridFileName[NumSets+sumidx][cidx], 
				   OutGrid, GridLat, GridLng, nrows, ncols, 
				   minlat, minlng, cellsize, NODATA, Ncells,
				   OVERWRITE);
      }
    }
  }

  return(0);

}

/**********************************************************************
**********************************************************************
**********************************************************************

  Subroutines and Functions Defined Below

**********************************************************************
**********************************************************************
**********************************************************************/

int GetStatInfo( char            *StatInfoParam, 
		 StatInfoStruct  *StatInfo, 
		 char           **ColNames, 
		 int              NumCols, 
		 int             *CalcPE, 
		 int             *CalcTR, 
		 int             *CalcTSM, 
		 int             *CalcMGDD, 
		 int             *CalcPD, 
		 int             *CalcCH,
		 int             *CalcBF,
		 int             *CalcDMWD,
		 int             *CalcDSFW,
		 char             ArcGrid, 
		 double          *GridLat,
		 double          *GridLng,
		 int              nrows, 
		 int              ncols,
		 double           minlat,
		 double           minlng,
		 double           cellsize,
		 int              NODATA) {
/**********************************************************************
  GetStatInfo                Keith Cherkauer           

  Read statistics definition file, and populated the StatInfoStruct
  with all required information.

  Modifications
  2008-May-13 to remove column number from input file.
  2017-May-31 to add flags for additional calculations.
  2017-Jun-12 to add two dimensional storage array for 
       threshold values, so files can be read in to set thresholds
       specific to each VIC model output cell. KAC

**********************************************************************/

  FILE *fin;
  int cidx, idx, tidx, Ncols, row, col, tmpCells;
  char tmpstr[MaxCharData], chkFile[MaxCharData];
  char shortName[MaxCharData], thresName[MaxCharData];
  char *last, *token;
  double tmpThres;

  // Open column statistics definition file
  if((fin=fopen(StatInfoParam,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open column statistics definition file %s for reading.\n",StatInfoParam);
    exit(0);
  }

  // Determine columns to be processed
  StatInfo->Ncols = 0;
  fgets(tmpstr, MaxCharData, fin);
  while ( !feof(fin) ) {
    StatInfo->Ncols++;
    fgets(tmpstr, MaxCharData, fin);
  }
  rewind(fin);
  StatInfo->ColNumList  = (int *)calloc(StatInfo->Ncols,sizeof(int));
  StatInfo->ColNameList = (char **)calloc(StatInfo->Ncols,sizeof(char *));
  StatInfo->ColStatList = (char **)calloc(StatInfo->Ncols,sizeof(char *));
  StatInfo->Thres       = (double ***)calloc(StatInfo->Ncols,sizeof(double**));
  StatInfo->Extra       = (double ***)calloc(StatInfo->Ncols,sizeof(double**));
  for ( cidx=0; cidx < StatInfo->Ncols; cidx++ ) {
    // allocate space for gridded threshold values
    StatInfo->Thres[cidx] = (double **)calloc(nrows,sizeof(double*));
    StatInfo->Extra[cidx] = (double **)calloc(nrows,sizeof(double*));
    for ( row = 0; row < nrows; row++ ) {
      StatInfo->Thres[cidx][row] = (double *)calloc(ncols,sizeof(double));
      StatInfo->Extra[cidx][row] = (double *)calloc(ncols,sizeof(double));
    }
  }
  StatInfo->MaxColNum   = -9;
  Ncols = 0;
  for ( cidx = StatInfo->Ncols-1; cidx >= 0; cidx-- ) {
    StatInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    StatInfo->ColStatList[cidx] = (char *)calloc(25,sizeof(char));
    fgets(tmpstr, MaxCharData, fin);
    sscanf( tmpstr, "%s %s", StatInfo->ColNameList[cidx], 
	    StatInfo->ColStatList[cidx] );
    if ( strncasecmp( StatInfo->ColStatList[cidx], "othres", 6 ) == 0 
	 || strncasecmp( StatInfo->ColStatList[cidx], "uthres", 6 ) == 0 
	 || strncasecmp( StatInfo->ColStatList[cidx], "lday", 4 ) == 0 
	 || strncasecmp( StatInfo->ColStatList[cidx], "fday", 4 ) == 0 
	 || strncasecmp( StatInfo->ColStatList[cidx], "days", 4 ) == 0 
	 || strncasecmp( StatInfo->ColStatList[cidx], "crossthres", 10 ) == 0
	 || strncasecmp( StatInfo->ColStatList[cidx], "avgdays", 7 ) == 0
	 || strncasecmp( StatInfo->ColStatList[cidx], "avgval", 6 ) == 0
	 || strncasecmp( StatInfo->ColStatList[cidx], "quan", 4 ) == 0 ) {
      // check if threshold is a file
      sscanf( tmpstr, "%*s %*s %s", chkFile );
      if ( strcasecmp( chkFile, "FILE" ) == 0 ) {
	// threshold is a valid file name, so read contents into threshold array
	sscanf( tmpstr, "%*s %*s %*s %s %s", shortName, thresName );
	// read spatially distributed thresholds from existing file
	tmpCells = ReadSpatialDataFiles( ArcGrid, thresName, 
					 StatInfo->Thres[cidx], GridLat, 
					 GridLng, nrows, ncols, minlat,
					 minlng, cellsize, NODATA);
	if ( tmpCells != nrows*ncols ) {
	  fprintf( stderr, "ERROR: Number of cells from threshold file does not match the number defined for the current domain.\n" );
	  exit (FAIL);
	}
	// build metric name with shortened name from stats control file
	sprintf( StatInfo->ColStatList[cidx], "%s_%s", StatInfo->ColStatList[cidx], shortName );
      }
      else {
	// threshold is not a valid filename, check for a constant value
	tmpThres = atof( chkFile );
	for ( row = 0; row < nrows; row ++ )
	  for ( col = 0; col < ncols; col++ )
	    StatInfo->Thres[cidx][row][col] = tmpThres;
	// add thres value to the statistic name so that output files can be differentiated
	sprintf( StatInfo->ColStatList[cidx], "%s_%g", StatInfo->ColStatList[cidx], tmpThres );

      }
    }
    // check if requested variable requires additional input to be computed
    if ( ( strncasecmp( StatInfo->ColNameList[cidx], "OUT_DMWD", 8 ) == 0 ) ||
	 ( strncasecmp( StatInfo->ColNameList[cidx], "OUT_PLANT_DAY", 13 ) == 0 ) ) {
      last = token = strtok(tmpstr, " ");
      for (;(token = strtok(NULL, " ")) != NULL; last = token);
      last = strip_copy(last);
      // read spatially distributed thresholds from existing file
      tmpCells = ReadSpatialDataFiles( ArcGrid, last, StatInfo->Extra[cidx],
				       GridLat, GridLng, nrows, ncols, minlat,
				       minlng, cellsize, NODATA);
      if ( tmpCells != nrows*ncols ) {
	fprintf( stderr, "ERROR: Number of cells from threshold file does not match the number defined for the current domain.\n" );
	exit (FAIL);
      }
    }
    // Make sure that request variable is in the input file
    for ( idx = 0; idx < NumCols; idx++ ) {
      if ( strcasecmp( ColNames[idx], StatInfo->ColNameList[cidx] ) == 0 ) {
	StatInfo->ColNumList[cidx] = idx;
	break;
      }
    }
    if ( idx == NumCols ) {
      if ( StatInfo->ColNameList[cidx] != NULL && StatInfo->ColNameList[cidx][0] != '#' )
	fprintf( stderr, "WARNING: Unable to find %s in the current output file, removing from analysis.\n", StatInfo->ColNameList[cidx] );
      for( tidx = cidx; tidx < StatInfo->Ncols-1; tidx++ ) {
	// shift all stats entries down one to remove the problem
	strcpy( StatInfo->ColNameList[tidx], StatInfo->ColNameList[tidx+1] );
	strcpy( StatInfo->ColStatList[tidx], StatInfo->ColStatList[tidx+1] );
	StatInfo->ColNumList[tidx] = StatInfo->ColNumList[tidx+1];
	StatInfo->Thres[tidx][row][col] = StatInfo->Thres[tidx+1][row][col];
      }
      StatInfo->Ncols--;
    }
    else if ( StatInfo->ColNumList[cidx] > StatInfo->MaxColNum ) 
      StatInfo->MaxColNum = StatInfo->ColNumList[cidx]+1;
  }
  fclose(fin);
  StatInfo->UseColFile = TRUE;
  StatInfo->Ncols -= Ncols;

  // check if any special calculations have been requested
  (*CalcPE) = (*CalcTR) = (*CalcTSM) = (*CalcMGDD) = (*CalcPD) = FALSE;
  (*CalcCH) = (*CalcBF) = (*CalcDMWD) = (*CalcDSFW) = FALSE;
  for ( cidx = StatInfo->Ncols-1; cidx >= 0; cidx-- ) {
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_PE" ) == 0 ) 
      (*CalcPE) = TRUE;
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_TOTAL_RUNOFF" ) == 0 ) 
      (*CalcTR) = TRUE;
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_TOTAL_SOIL_MOIST" ) == 0 ) 
      (*CalcTSM) = TRUE;
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_MGDD" ) == 0 ) 
      (*CalcMGDD) = TRUE;
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_PLANT_DAY" ) == 0 ) 
      (*CalcPD) = TRUE;
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_CHILL_HR" ) == 0 ) 
      (*CalcCH) = TRUE;
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_BUD_FREEZE" ) == 0 ) 
      (*CalcBF) = TRUE;
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_DMWD" ) == 0 ) 
      (*CalcDMWD) = TRUE;
    if ( strcasecmp( StatInfo->ColNameList[cidx], "OUT_DSFW" ) == 0 ) 
      (*CalcDSFW) = TRUE;
  }

  return(0);
}

int GetPenmanInfo( PenInfoStruct  *PenInfo, 
		   char          **ColNames, 
		   int             NumCols ) {
/**********************************************************************
  GetPenmanInfo     Keith Cherkauer           May 14, 2008

  Determine if the required columns are in the set of model output
  files to complete an off-line calculation of Potential Evaporation 
  using the Penman method.  If so, then store the required information
  for use later in the processing script.

  Modifications
  2008-May-14 copied from GetStatInfo.
  2008-May-28 to remove need for separate file to define penman variables.
  2016-Jan-12 to add OUT_PE as a variable column name

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation for PE!
  static char *PenmanList[] = { "OUT_R_NET", "OUT_GRND_FLUX", "OUT_WIND", "OUT_SURF_TEMP", "OUT_REL_HUMID", "OUT_AIR_TEMP", "OUT_PE" };

  int cidx, idx;

  // Determine columns to be processed
  PenInfo->Ncols = 7;
  PenInfo->ColNumList  = (int *)calloc(PenInfo->Ncols,sizeof(int));
  PenInfo->ColNameList = (char **)calloc(PenInfo->Ncols,sizeof(char *));
  PenInfo->MaxColNum   = -9;
 
  for ( cidx = 0; cidx < PenInfo->Ncols; cidx++ ) {
    PenInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    strcpy( PenInfo->ColNameList[cidx], PenmanList[cidx] );
    for ( idx = 0; idx < NumCols; idx++ ) {
      if ( strcasecmp( ColNames[idx], PenInfo->ColNameList[cidx] ) == 0 ) {
	PenInfo->ColNumList[cidx] = idx;
	break;
      }
    }
    if ( idx == NumCols ) {
      // handle missing column by turning off calculation
      fprintf( stderr, "Columns required for PE calculation not found, so will not be computed.\n" );
      PenInfo->Ncols = 0;
      return (FALSE);
    }
    else if ( PenInfo->ColNumList[cidx] > PenInfo->MaxColNum ) 
      PenInfo->MaxColNum = PenInfo->ColNumList[cidx]+1;
  }
  
  return (TRUE);

}

int GetTotalRunoffInfo( PenInfoStruct  *TRoffInfo, 
			char          **ColNames, 
			int             NumCols ) {
/**********************************************************************
  GetTotalRunoffInfo     Keith Cherkauer           November 12, 2008

  Determine if the required columns are in the set of model output
  files to complete an off-line calculation of Total Runoff (surface
  runoff + baseflow).  If so, then store the required information
  for use later in the processing script.

  Modifications
  2008-Nov-12 copied from GetPenInfo.

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation for total runoff!
  static char *TotalRunoffList[] = { "OUT_RUNOFF", "OUT_BASEFLOW", "OUT_TOTAL_RUNOFF" };

  int cidx, idx;

  // Determine columns to be processed
  TRoffInfo->Ncols = 3;
  TRoffInfo->ColNumList  = (int *)calloc(TRoffInfo->Ncols,sizeof(int));
  TRoffInfo->ColNameList = (char **)calloc(TRoffInfo->Ncols,sizeof(char *));
  TRoffInfo->MaxColNum   = -9;
 
 for ( cidx = 0; cidx < TRoffInfo->Ncols; cidx++ ) {
   TRoffInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
   strcpy( TRoffInfo->ColNameList[cidx], TotalRunoffList[cidx] );
   for ( idx = 0; idx < NumCols; idx++ ) {
     if ( strcasecmp( ColNames[idx], TRoffInfo->ColNameList[cidx] ) == 0 ) {
       TRoffInfo->ColNumList[cidx] = idx;
       break;
     }
   }
   if ( idx == NumCols ) {
      // handle missing column by turning off calculation
      fprintf( stderr, "Columns required for TOTAL_RUNOFF calculation not found, so will not be computed.\n" );
      TRoffInfo->Ncols = 0;
      return (FALSE);
   }
   if ( TRoffInfo->ColNumList[cidx] > TRoffInfo->MaxColNum ) 
      TRoffInfo->MaxColNum = TRoffInfo->ColNumList[cidx]+1;
  }

  return (TRUE);
}

int GetTotalSoilMoistureInfo( PenInfoStruct  *TSMInfo, 
			      char          **ColNames, 
			      int             NumCols,
			      int             NumLayers ) {
/**********************************************************************
  GetTotalSoilMoistureInfo     Keith Cherkauer           June 16, 2017

  Determine if the required columns are in the set of model output
  files to compute total soil moisture (sum of soil moisture in all
  soil layers - NOTE only works if soil moisutre reported in mm, not as
  a fraction).  If so, then store the required information for use 
  later in the processing script.

  Modifications
  2017-Jun-16 copied from GetTotalRunoffInfo.

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation for total soil moisture!
  static char **TotalSoilMoistureList;
  char tmpstr[30];

  int cidx, idx;

  // Determine columns to be processed
  TSMInfo->Ncols = NumLayers+1;
  TSMInfo->ColNumList  = (int *)calloc(TSMInfo->Ncols,sizeof(int));
  TSMInfo->ColNameList = (char **)calloc(TSMInfo->Ncols,sizeof(char *));
  TSMInfo->MaxColNum   = -9;
 
  // define column names - DO NOT MODIFY THIS LIST - unless you are also changing the calculation
  TotalSoilMoistureList = (char **)calloc( TSMInfo->Ncols, sizeof(char *) );
  for ( idx = 0; idx < TSMInfo->Ncols; idx++ ) 
    TotalSoilMoistureList[idx] = (char *)calloc( 30, sizeof(char) );
  strcpy( TotalSoilMoistureList[NumLayers], "OUT_TOTAL_SOIL_MOIST" );
  for ( idx = 0; idx < NumLayers; idx++ )
    sprintf( TotalSoilMoistureList[idx], "OUT_SOIL_MOIST_%i", idx );

 for ( cidx = 0; cidx < TSMInfo->Ncols; cidx++ ) {
   TSMInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
   strcpy( TSMInfo->ColNameList[cidx], TotalSoilMoistureList[cidx] );
   for ( idx = 0; idx < NumCols; idx++ ) {
     if ( strcasecmp( ColNames[idx], TSMInfo->ColNameList[cidx] ) == 0 ) {
       TSMInfo->ColNumList[cidx] = idx;
       break;
     }
   }
   if ( idx == NumCols ) {
      // handle missing column by turning off calculation
      fprintf( stderr, "Columns required for TOTAL_SOIL_MOIST calculation not found, so will not be computed.\n" );
      TSMInfo->Ncols = 0;
      return (FALSE);
   }
   if ( TSMInfo->ColNumList[cidx] > TSMInfo->MaxColNum ) 
      TSMInfo->MaxColNum = TSMInfo->ColNumList[cidx]+1;
  }

  return (TRUE);
}

int GetModifiedGrowingDegreeDayInfo( PenInfoStruct   *MGDDInfo, 
				     StatInfoStruct  *StatInfo, 
				     char           **ColNames, 
				     int              NumCols ) {
/**********************************************************************
  GetModifiedGrowingDegreeDayInfo     Keith Cherkauer           May 30, 2017

  Determine if the required columns are in the set of model output
  files to complete an off-line calculation of modified Growing Degree 
  Days using the method defined by the Midwestern Regional Climate Center 
  (MRCC).  If so, then store the required information for use later in 
  the processing script.

  Modifications
  2017-May-30 Copied from GetTotalRunoffInfo   KAC

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation for total runoff!
  static char *MGDDList[] = { "IN_TMIN", "IN_TMAX", "OUT_MGDD" };

  int cidx, sidx, idx;

  // Determine columns to be processed
  MGDDInfo->Ncols = 3;
  MGDDInfo->ColNumList  = (int *)calloc(MGDDInfo->Ncols,sizeof(int));
  MGDDInfo->ColNameList = (char **)calloc(MGDDInfo->Ncols,sizeof(char *));
  MGDDInfo->MaxColNum   = -9;

  // Identify required columns in the input file
  for ( cidx = 0; cidx < MGDDInfo->Ncols; cidx++ ) {
    MGDDInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    strcpy( MGDDInfo->ColNameList[cidx], MGDDList[cidx] );
    for ( sidx=0; sidx<StatInfo->Ncols; sidx++ )
      // check for variable in StatInfo table
      if ( strcasecmp( StatInfo->ColNameList[sidx], MGDDInfo->ColNameList[cidx] ) == 0 ) 
	break;
    if ( sidx < StatInfo->Ncols ) 
      // variable found in current table
      MGDDInfo->ColNumList[cidx] = sidx;
    else {
      // variable not found in current table, so check if it is in input file
      for ( idx = 0; idx < NumCols; idx++ )
	if ( strcasecmp( ColNames[idx], MGDDInfo->ColNameList[cidx] ) == 0 )
	  break;
      if ( idx == NumCols ) {
	// handle missing column by turning off calculation
	fprintf( stderr, "Columns required for Modified Growing Degree Day calculation not found, so will not be computed.\n" );
	MGDDInfo->Ncols = 0;
	return (FALSE);
      }
      // add missing variable to the StatInfo table (must reallocate memory)
      StatInfo->Ncols++;
      if ( idx > StatInfo->MaxColNum ) StatInfo->MaxColNum = idx;
      StatInfo->ColNameList = (char **)realloc( StatInfo->ColNameList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColNameList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColNameList[StatInfo->Ncols-1], MGDDInfo->ColNameList[cidx] );
      StatInfo->ColStatList = (char **)realloc( StatInfo->ColStatList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColStatList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColStatList[StatInfo->Ncols-1], "none" );
      StatInfo->ColNumList = (int *)realloc( StatInfo->ColNumList, StatInfo->Ncols*sizeof(int) );
      StatInfo->ColNumList[StatInfo->Ncols-1] = idx;
      MGDDInfo->ColNumList[cidx] = StatInfo->Ncols-1;
    }
    // reset max column setting for current statistic
    if ( MGDDInfo->ColNumList[cidx] > MGDDInfo->MaxColNum ) 
      MGDDInfo->MaxColNum = MGDDInfo->ColNumList[cidx]+1;
  }

  // check is the calculated statistics is in Statinfo more than once
  MGDDInfo->OutputCols = (int *)calloc(StatInfo->Ncols,sizeof(int));
  MGDDInfo->Noutput = 0;
  for ( sidx=0; sidx<StatInfo->Ncols; sidx++ ) {
    if ( strcasecmp( StatInfo->ColNameList[sidx], MGDDInfo->ColNameList[MGDDInfo->Ncols-1] ) == 0 ) {
      // for each occurance of stat name add column to list
      MGDDInfo->OutputCols[MGDDInfo->Noutput] = sidx;
      MGDDInfo->Noutput++;
    }
  }

  return (TRUE);

}

int CalcModifiedGrowingDegreeDays( double      *MGDD, 
				   double      *TMIN, 
				   double      *TMAX, 
				   DATE_STRUCT  startdate, 
				   int          Nvals ) {
/**********************************************************************
  CalcModifiedGrowingDegreeDays     Keith Cherkauer           May 30, 2017

  This subroutine calculates Modified Growing Degree Days using the 
  method defined by the Midwestern Regional Climate Center (MRCC) as 
  of May 2017
  Source: mrcc.isws.illinois.edu/cliwatch/mgdd/
  Based on critical temperature of 50F for corn and soybean

  Modifications

**********************************************************************/

  int ErrNum=0;
  int didx;
  double tmpTMIN, tmpTMAX;
  double tmpStartJulDate, tmpEndJulDate;
  double BaseT = 50.;

  tmpStartJulDate = calc_juldate( startdate.year, 4, 1, 0 );
  tmpEndJulDate = calc_juldate( startdate.year, 12, 1, 0 );
  for ( didx = 0; didx < Nvals; didx++ ) {
    // initialize MGDD
    if ( didx == 0 ) MGDD[didx] = 0;
    else MGDD[didx] = MGDD[didx-1];
    // process MGDD within the growing season
    if ( startdate.juldate + didx >= tmpStartJulDate 
	 && startdate.juldate + didx < tmpEndJulDate ) {
      // within the valid season, so convert temps to degree F
      tmpTMIN = ( TMIN[didx] * 9. ) / 5. + 32.;
      tmpTMAX = ( TMAX[didx] * 9. ) / 5. + 32.;
      // check for temperature extremes that exceed calculation range
      if ( tmpTMAX < BaseT ) tmpTMAX = BaseT;
      if ( tmpTMIN < BaseT ) tmpTMIN = BaseT;
      if ( tmpTMIN > 86. ) tmpTMIN = 86.;
      if ( tmpTMAX > 86. ) tmpTMAX = 86.;
      // calculate current MGDD
      MGDD[didx] += (tmpTMAX + tmpTMIN) / 2. - BaseT;
    }
  }

  return 0;

}

int GetPlantingDayInfo( PenInfoStruct   *PDayInfo, 
			StatInfoStruct  *StatInfo, 
			char           **ColNames, 
			int              NumCols ) {
/**********************************************************************
  GetPlantingDayInfo     Keith Cherkauer           June 19, 2017


  Modifications
  2017-May-30 Copied from GetTotalRunoffInfo   KAC

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation!
  static char *PDayList[] = { "OUT_SOIL_TEMP_0", "OUT_SOIL_MOIST_0", "OUT_PLANT_DAY" };

  int cidx, sidx, idx;

  // Determine columns to be processed
  PDayInfo->Ncols = 3;
  PDayInfo->ColNumList  = (int *)calloc(PDayInfo->Ncols,sizeof(int));
  PDayInfo->ColNameList = (char **)calloc(PDayInfo->Ncols,sizeof(char *));
  PDayInfo->MaxColNum   = -9;

  // Identify required columns in the input file
  for ( cidx = 0; cidx < PDayInfo->Ncols; cidx++ ) {
    PDayInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    strcpy( PDayInfo->ColNameList[cidx], PDayList[cidx] );
    for ( sidx=0; sidx<StatInfo->Ncols; sidx++ )
      // check for variable in StatInfo table
      if ( strcasecmp( StatInfo->ColNameList[sidx], PDayInfo->ColNameList[cidx] ) == 0 ) 
	break;
    if ( sidx < StatInfo->Ncols ) 
      // variable found in current table
      PDayInfo->ColNumList[cidx] = sidx;
    else {
      // variable not found in current table, so check if it is in input file
      for ( idx = 0; idx < NumCols; idx++ )
	if ( strcasecmp( ColNames[idx], PDayInfo->ColNameList[cidx] ) == 0 )
	  break;
      if ( idx == NumCols ) {
	// handle missing column by turning off calculation
	fprintf( stderr, "Columns required for Modified Growing Degree Day calculation not found, so will not be computed.\n" );
	PDayInfo->Ncols = 0;
	return (FALSE);
      }
      // add missing variable to the StatInfo table (must reallocate memory)
      StatInfo->Ncols++;
      if ( idx > StatInfo->MaxColNum ) StatInfo->MaxColNum = idx;
      StatInfo->ColNameList = (char **)realloc( StatInfo->ColNameList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColNameList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColNameList[StatInfo->Ncols-1], PDayInfo->ColNameList[cidx] );
      StatInfo->ColStatList = (char **)realloc( StatInfo->ColStatList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColStatList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColStatList[StatInfo->Ncols-1], "none" );
      StatInfo->ColNumList = (int *)realloc( StatInfo->ColNumList, StatInfo->Ncols*sizeof(int) );
      StatInfo->ColNumList[StatInfo->Ncols-1] = idx;
      PDayInfo->ColNumList[cidx] = StatInfo->Ncols-1;
    }
    // reset max column setting for current statistic
    if ( PDayInfo->ColNumList[cidx] > PDayInfo->MaxColNum ) 
      PDayInfo->MaxColNum = PDayInfo->ColNumList[cidx]+1;
  }

  // check is the calculated statistics is in Statinfo more than once
  PDayInfo->OutputCols = (int *)calloc(StatInfo->Ncols,sizeof(int));
  PDayInfo->Noutput = 0;
  for ( sidx=0; sidx<StatInfo->Ncols; sidx++ ) {
    if ( strcasecmp( StatInfo->ColNameList[sidx], PDayInfo->ColNameList[PDayInfo->Ncols-1] ) == 0 ) {
      // for each occurance of stat name add column to list
      PDayInfo->OutputCols[PDayInfo->Noutput] = sidx;
      PDayInfo->Noutput++;
    }
  }

  return (TRUE);

}

int CalcPlantingDay( double      *PDay, 
		     double      *TEMP, 
		     double      *MOIST, 
		     double       FieldCapacity,
		     DATE_STRUCT  startdate, 
		     int          Nvals ) {
/**********************************************************************
  CalcPlantingDay     Keith Cherkauer           June 19, 2017


  Modifications

**********************************************************************/

  int ErrNum=0;
  int didx;

  for ( didx = 0; didx < Nvals; didx++ ) {
    if ( TEMP[didx] >= 10. && MOIST[didx] < FieldCapacity ) PDay[didx] = 1;
  }

  return 0;

}

int GetChillingHoursInfo( PenInfoStruct   *ChillHrInfo, 
			  StatInfoStruct  *StatInfo, 
			  char           **ColNames, 
			  int              NumCols ) {
/**********************************************************************
  GetChillingHoursInfo     Keith Cherkauer           May 30, 2017

  Determine if the required columns are in the set of model output
  files to complete an off-line calculation of Chilling Hours as defined
  by Janna Beckerman (Purdeu, HORT) for the INCCIA 2017.   Chilling hours
  refer to required low temperatures to set fruit tree buds.  If so, then 
  store the required information for use later in the processing script.

  Modifications
  2017-May-30 copied from GetTotalRunoffInfo   KAC

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation for total runoff!
  static char *ChillHrList[] = { "IN_TMIN", "IN_TMAX", "OUT_CHILL_HR" };

  int cidx, sidx, idx;

  // Determine columns to be processed
  ChillHrInfo->Ncols = 3;
  ChillHrInfo->ColNumList  = (int *)calloc(ChillHrInfo->Ncols,sizeof(int));
  ChillHrInfo->ColNameList = (char **)calloc(ChillHrInfo->Ncols,sizeof(char *));
  ChillHrInfo->MaxColNum   = -9;
 
  // Identify required columns in the input file
  for ( cidx = 0; cidx < ChillHrInfo->Ncols; cidx++ ) {
    ChillHrInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    strcpy( ChillHrInfo->ColNameList[cidx], ChillHrList[cidx] );
    for ( sidx=0; sidx<StatInfo->Ncols; sidx++ )
      // check for variable in StatInfo table
      if ( strcasecmp( StatInfo->ColNameList[sidx], ChillHrInfo->ColNameList[cidx] ) == 0 ) 
	break;
    if ( sidx < StatInfo->Ncols ) 
      // variable found in current table
      ChillHrInfo->ColNumList[cidx] = sidx;
    else {
      // variable not found in current table, so check if it is in input file
      for ( idx = 0; idx < NumCols; idx++ )
	if ( strcasecmp( ColNames[idx], ChillHrInfo->ColNameList[cidx] ) == 0 )
	  break;
      if ( idx == NumCols ) {
	// handle missing column by turning off calculation
	fprintf( stderr, "Columns required for Chilling Hour calculation not found, so will not be computed.\n" );
	ChillHrInfo->Ncols = 0;
	return (FALSE);
      }
      // add missing variable to the StatInfo table (must reallocate memory)
      StatInfo->Ncols++;
      if ( idx > StatInfo->MaxColNum ) StatInfo->MaxColNum = idx;
      StatInfo->ColNameList = (char **)realloc( StatInfo->ColNameList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColNameList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColNameList[StatInfo->Ncols-1], ChillHrInfo->ColNameList[cidx] );
      StatInfo->ColStatList = (char **)realloc( StatInfo->ColStatList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColStatList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColStatList[StatInfo->Ncols-1], "none" );
      StatInfo->ColNumList = (int *)realloc( StatInfo->ColNumList, StatInfo->Ncols*sizeof(int) );
      StatInfo->ColNumList[StatInfo->Ncols-1] = idx;
      ChillHrInfo->ColNumList[cidx] = StatInfo->Ncols-1;
    }
    // reset max column setting for current statistic
    if ( ChillHrInfo->ColNumList[cidx] > ChillHrInfo->MaxColNum ) 
      ChillHrInfo->MaxColNum = ChillHrInfo->ColNumList[cidx]+1;
 }
 
  // check is the calculated statistics is in Statinfo more than once
  ChillHrInfo->OutputCols = (int *)calloc(StatInfo->Ncols,sizeof(int));
  ChillHrInfo->Noutput = 0;
  for ( sidx=0; sidx<StatInfo->Ncols; sidx++ ) {
    if ( strcasecmp( StatInfo->ColNameList[sidx], ChillHrInfo->ColNameList[ChillHrInfo->Ncols-1] ) == 0 ) {
      // for each occurance of stat name add column to list
      ChillHrInfo->OutputCols[ChillHrInfo->Noutput] = sidx;
      ChillHrInfo->Noutput++;
    }
  }

 return (TRUE);

}

int CalcChillingHours( double *ChillHr, 
		       double *TMIN, 
		       double *TMAX, 
		       DATE_STRUCT startdate, 
		       int Nvals ) {
/**********************************************************************
  CalcChillingHours     Keith Cherkauer           May 30, 2017

  This subroutine calculates Chilling Hours using the definition 
  developed as part of the Agriculture Working Group of the 2017 INCCIA.
  Relies on a simple triangular (Sanders, 1975) approximation to get 
  hourly temperatures from daily maximum and minimum values.  

  Chilling hours evalaute the accumulation of cold required to set the 
  buds of fruit trees and bushes.

  Source: Janna Beckerman, Purdue

  Modifications

**********************************************************************/

  int ErrNum=0;
  int didx, hidx;
  double tmpTMIN, tmpTMAX;
  double tmpStartJulDate, tmpEndJulDate;
  double *tmpTair;

  // estimate hourly temperatures based on sawtooth method (Sanders, 1975)
  // assumes TMIN occurs at 0500 and TMAX at 1500 and temperature change 
  // is linear between those extremes.
  tmpTair = (double *)calloc(Nvals*24,sizeof(double));
  for ( didx=0; didx<Nvals; didx++ ) {
    for( hidx=0; hidx<5; hidx++ ) {
      // fill in evening to morning minimum
      if ( didx == 0 ) tmpTair[didx*24+hidx] = TMIN[didx];
      else tmpTair[didx*24+hidx] = ((double)hidx-29.)/(15.-29.)*(TMAX[didx-1]-TMIN[didx])+TMIN[didx];
    }
    tmpTair[didx*24+5]=TMIN[didx]; // set minimum temperature
    for ( hidx=6; hidx<15; hidx++ ) {
      // fill in daily temperature increase
      tmpTair[didx*24+hidx] = ((double)hidx-15.)/(5.-15.)*(TMIN[didx]-TMAX[didx])+TMAX[didx];
    }
    tmpTair[didx*24+15]=TMAX[didx]; // set maximum temperature
    for( hidx=16; hidx<24; hidx++ ) {
      // fill in afternoon to evening temperatures
      if ( didx == Nvals-1 ) tmpTair[didx*24+hidx] = TMAX[didx];
      else tmpTair[didx*24+hidx] = ((double)hidx-29.)/(15.-29.)*(TMAX[didx]-TMIN[didx+1])+TMIN[didx+1];
    }    
  }

  // Compute chilling hours per day
  for ( didx=0; didx<Nvals; didx++ ) {
    ChillHr[didx] = 0;
    for ( hidx=0; hidx<24; hidx++ ) {
      // convert temps to degree F
      tmpTair[didx*24+hidx] = ( tmpTair[didx*24+hidx] * 9. ) / 5. + 32.;
      // compute chilling hours
      if ( tmpTair[didx*24+hidx] >= 35. && tmpTair[didx*24+hidx] <= 45 )
	ChillHr[didx]++;
    }
    if ( didx>0 ) ChillHr[didx] += ChillHr[didx-1];
  }

  free(tmpTair);

  return ErrNum;

}

int GetBudFreezeInfo( PenInfoStruct   *BudFreezeInfo, 
		      StatInfoStruct  *StatInfo, 
		      char           **ColNames, 
		      int              NumCols ) {
/**********************************************************************
  GetBudFreezeInfo     Keith Cherkauer           June 20, 2017

  Determine if the required columns are in the set of model output
  files to complete an off-line calculation of Bud Freeze, which 
  requires calculation of Chilling Hours as defined by Janna Beckerman 
  (Purdue, HORT) for the INCCIA 2017.   If so, then store the required 
  information for use later in the processing script.

  Modifications
  2017-May-30 copied from GetTotalRunoffInfo   KAC

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation for total runoff!
  static char *BudFreezeList[] = { "IN_TMIN", "IN_TMAX", "OUT_BUD_FREEZE" };

  int cidx, sidx, idx;

  // Determine columns to be processed
  BudFreezeInfo->Ncols = 3;
  BudFreezeInfo->ColNumList  = (int *)calloc(BudFreezeInfo->Ncols,sizeof(int));
  BudFreezeInfo->ColNameList = (char **)calloc(BudFreezeInfo->Ncols,sizeof(char *));
  BudFreezeInfo->MaxColNum   = -9;
 
  // Identify required columns in the input file
  for ( cidx = 0; cidx < BudFreezeInfo->Ncols; cidx++ ) {
    BudFreezeInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    strcpy( BudFreezeInfo->ColNameList[cidx], BudFreezeList[cidx] );
    for ( sidx=0; sidx<StatInfo->Ncols; sidx++ )
      // check for variable in StatInfo table
      if ( strcasecmp( StatInfo->ColNameList[sidx], BudFreezeInfo->ColNameList[cidx] ) == 0 ) 
	break;
    if ( sidx < StatInfo->Ncols ) 
      // variable found in current table
      BudFreezeInfo->ColNumList[cidx] = sidx;
    else {
      // variable not found in current table, so check if it is in input file
      for ( idx = 0; idx < NumCols; idx++ )
	if ( strcasecmp( ColNames[idx], BudFreezeInfo->ColNameList[cidx] ) == 0 )
	  break;
      if ( idx == NumCols ) {
	// handle missing column by turning off calculation
	fprintf( stderr, "Columns required for Chilling Hour calculation not found, so will not be computed.\n" );
	BudFreezeInfo->Ncols = 0;
	return (FALSE);
      }
      // add missing variable to the StatInfo table (must reallocate memory)
      StatInfo->Ncols++;
      if ( idx > StatInfo->MaxColNum ) StatInfo->MaxColNum = idx;
      StatInfo->ColNameList = (char **)realloc( StatInfo->ColNameList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColNameList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColNameList[StatInfo->Ncols-1], BudFreezeInfo->ColNameList[cidx] );
      StatInfo->ColStatList = (char **)realloc( StatInfo->ColStatList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColStatList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColStatList[StatInfo->Ncols-1], "none" );
      StatInfo->ColNumList = (int *)realloc( StatInfo->ColNumList, StatInfo->Ncols*sizeof(int) );
      StatInfo->ColNumList[StatInfo->Ncols-1] = idx;
      BudFreezeInfo->ColNumList[cidx] = StatInfo->Ncols-1;
    }
    // reset max column setting for current statistic
    if ( BudFreezeInfo->ColNumList[cidx] > BudFreezeInfo->MaxColNum ) 
      BudFreezeInfo->MaxColNum = BudFreezeInfo->ColNumList[cidx]+1;
 }
 
  // check is the calculated statistics is in Statinfo more than once
  BudFreezeInfo->OutputCols = (int *)calloc(StatInfo->Ncols,sizeof(int));
  BudFreezeInfo->Noutput = 0;
  for ( sidx=0; sidx<StatInfo->Ncols; sidx++ ) {
    if ( strcasecmp( StatInfo->ColNameList[sidx], BudFreezeInfo->ColNameList[BudFreezeInfo->Ncols-1] ) == 0 ) {
      // for each occurance of stat name add column to list
      BudFreezeInfo->OutputCols[BudFreezeInfo->Noutput] = sidx;
      BudFreezeInfo->Noutput++;
    }
  }

 return (TRUE);

}

int CalcBudFreeze( double *BudFreeze, 
		   double *TMIN, 
		   double *TMAX, 
		   DATE_STRUCT startdate, 
		   int Nvals ) {
/**********************************************************************
  CalcBudFreeze     Keith Cherkauer           June 20, 2017

  This subroutine calculates Chilling Hours then evaluates the 
  potential of bud freeze, by looking at 5 day average temperatures
  cold enough to freeze (kill) the bud that has been set once chilling
  hours have accumulated to 1200.  Developed as part of the Agriculture 
  Working Group of the 2017 INCCIA.  

  Chilling hours evalaute the accumulation of cold required to set the 
  buds of fruit trees and bushes.

  Source: Janna Beckerman, Purdue

  Modifications

**********************************************************************/

  char ColdEnough = FALSE;
  char WarmEnough = FALSE;
  char FrozenBud = FALSE;
  int ErrNum=0;
  int didx, idx, WarmDays;
  double *tmpChillHr;

  tmpChillHr = (double *)calloc( Nvals, sizeof(double) );
  ErrNum = CalcChillingHours( tmpChillHr, TMIN, TMAX, startdate, Nvals );

  for ( idx = 0; idx < Nvals; idx++ ) BudFreeze[didx] = 0;

  // Look for periods with potential killing frost after accumultion of 1200 chilling hours
  didx = 0;
  // check that at least 1200 chilling hours have accumulated
  while ( didx < Nvals && tmpChillHr[didx] < 1200 ) didx++;
  WarmDays = 0;
  if ( didx < Nvals ) {
    ColdEnough = TRUE;
    while ( didx < Nvals && !WarmEnough ) {
      // check that at least 14 days of warm weather have occurred
      if ( TMIN[didx] > 7.2 ) WarmDays++;
      if ( WarmDays >= 14 ) WarmEnough = TRUE;
      didx++;
    }
    if ( WarmEnough ) {
      // check for freezing conditions with sensitive buds
      while ( didx < Nvals && !FrozenBud ) {
	if ( TMIN[didx] < -2.78 ) {
	  FrozenBud = TRUE;
	  BudFreeze[didx] = 1;
	}
	didx++;
      }
      if ( FrozenBud )
	for ( idx = didx; idx < Nvals; idx++ ) BudFreeze[didx] = 1;
    }
  }

  free(tmpChillHr);

  return ErrNum;

}

int GetDrainModWorkingDaysInfo( PenInfoStruct   *DModWDaysInfo, 
				StatInfoStruct  *StatInfo, 
				char           **ColNames, 
				int              NumCols ) {
  /**********************************************************************
  GetDrainModWorkingDaysInfo     Keith Cherkauer           June 19, 2017

  Number of days per week when soil moisture in the top layer < 
  max moisture - 20 mm, Pi < 6 mm & Pi-1 < 6cm (Ale et al. 2009)

  Modifications
  2017-May-30 copied from GetTotalRunoffInfo   KAC

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation for total runoff!
  static char *DModWDaysList[] = { "OUT_PREC", "OUT_SOIL_MOIST_0", "OUT_DMWD" };

  int cidx, sidx, idx;

  // Determine columns to be processed
  DModWDaysInfo->Ncols = 3;
  DModWDaysInfo->ColNumList  = (int *)calloc(DModWDaysInfo->Ncols,sizeof(int));
  DModWDaysInfo->ColNameList = (char **)calloc(DModWDaysInfo->Ncols,sizeof(char *));
  DModWDaysInfo->MaxColNum   = -9;
 
  // Identify required columns in the input file
  for ( cidx = 0; cidx < DModWDaysInfo->Ncols; cidx++ ) {
    DModWDaysInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    strcpy( DModWDaysInfo->ColNameList[cidx], DModWDaysList[cidx] );
    for ( sidx=0; sidx<StatInfo->Ncols; sidx++ )
      // check for variable in StatInfo table
      if ( strcasecmp( StatInfo->ColNameList[sidx], DModWDaysInfo->ColNameList[cidx] ) == 0 ) 
	break;
    if ( sidx < StatInfo->Ncols ) 
      // variable found in current table
      DModWDaysInfo->ColNumList[cidx] = sidx;
    else {
      // variable not found in current table, so check if it is in input file
      for ( idx = 0; idx < NumCols; idx++ )
	if ( strcasecmp( ColNames[idx], DModWDaysInfo->ColNameList[cidx] ) == 0 )
	  break;
      if ( idx == NumCols ) {
	// handle missing column by turning off calculation
	fprintf( stderr, "Columns required for Chilling Hour calculation not found, so will not be computed.\n" );
	DModWDaysInfo->Ncols = 0;
	return (FALSE);
      }
      // add missing variable to the StatInfo table (must reallocate memory)
      StatInfo->Ncols++;
      if ( idx > StatInfo->MaxColNum ) StatInfo->MaxColNum = idx;
      StatInfo->ColNameList = (char **)realloc( StatInfo->ColNameList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColNameList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColNameList[StatInfo->Ncols-1], DModWDaysInfo->ColNameList[cidx] );
      StatInfo->ColStatList = (char **)realloc( StatInfo->ColStatList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColStatList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColStatList[StatInfo->Ncols-1], "none" );
      StatInfo->ColNumList = (int *)realloc( StatInfo->ColNumList, StatInfo->Ncols*sizeof(int) );
      StatInfo->ColNumList[StatInfo->Ncols-1] = idx;
      DModWDaysInfo->ColNumList[cidx] = StatInfo->Ncols-1;
    }
    // reset max column setting for current statistic
    if ( DModWDaysInfo->ColNumList[cidx] > DModWDaysInfo->MaxColNum ) 
      DModWDaysInfo->MaxColNum = DModWDaysInfo->ColNumList[cidx]+1;
 }
 
  // check is the calculated statistics is in Statinfo more than once
  DModWDaysInfo->OutputCols = (int *)calloc(StatInfo->Ncols,sizeof(int));
  DModWDaysInfo->Noutput = 0;
  for ( sidx=0; sidx<StatInfo->Ncols; sidx++ ) {
    if ( strcasecmp( StatInfo->ColNameList[sidx], DModWDaysInfo->ColNameList[DModWDaysInfo->Ncols-1] ) == 0 ) {
      // for each occurance of stat name add column to list
      DModWDaysInfo->OutputCols[DModWDaysInfo->Noutput] = sidx;
      DModWDaysInfo->Noutput++;
    }
  }

 return (TRUE);

}

int CalcDrainModWorkingDays( double *DModWDays, 
			     double *PREC, 
			     double *MOIST, 
			     double  MaxMoist,
			     DATE_STRUCT startdate, 
			     int Nvals ) {
/**********************************************************************
  CalcDrainModWorkingDays     Keith Cherkauer           June 19, 2017


  Source: Ale et al., 2009

  Modifications

**********************************************************************/

  int ErrNum=0;
  int didx, hidx;
  double tmpTMIN, tmpTMAX;
  double tmpStartJulDate, tmpEndJulDate;
  double *tmpTair;

  // Compute DRAINMOD working days
  DModWDays[0] = 0;
  for ( didx=1; didx<Nvals; didx++ ) {
    DModWDays[didx] = 0;
    if ( MOIST[didx] < MaxMoist - 20. ) {
      if ( PREC[didx] < 6. && PREC[didx-1] < 60. ) {
	// check precipitation condition - soil is considered workable
	DModWDays[didx] = 1;
      }
    }
  }

  return ErrNum;

}

int GetDaysSuitable4FieldWorkInfo( PenInfoStruct   *DSFWInfo, 
				   StatInfoStruct  *StatInfo, 
				   char           **ColNames, 
				   int              NumCols ) {
  /**********************************************************************
  GetDaysSuitable4FieldWorkInfo     Keith Cherkauer           June 19, 2017

  Ben Gramig's regression model to determine the workability of 
  agricultural fields.

  Source: Gramig, Massey and Yun (2017). Nitrogen application 
    decision-making under climate risk in the US Corn Belt, 
    Climate Risk Management, 15, 82-89.

  Modifications
  2017-May-30 copied from GetTotalRunoffInfo   KAC

**********************************************************************/

  // DO NOT MODIFY THIS LIST - unless you are also changing the calculation for total runoff!
  static char *DSFWList[] = { "IN_PREC", "IN_TMAX", "IN_TMIN", "OUT_DSFW" };

  int cidx, sidx, idx;

  // Determine columns to be processed
  DSFWInfo->Ncols = 4;
  DSFWInfo->ColNumList  = (int *)calloc(DSFWInfo->Ncols,sizeof(int));
  DSFWInfo->ColNameList = (char **)calloc(DSFWInfo->Ncols,sizeof(char *));
  DSFWInfo->MaxColNum   = -9;
 
  // Identify required columns in the input file
  for ( cidx = 0; cidx < DSFWInfo->Ncols; cidx++ ) {
    DSFWInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    strcpy( DSFWInfo->ColNameList[cidx], DSFWList[cidx] );
    for ( sidx=0; sidx<StatInfo->Ncols; sidx++ )
      // check for variable in StatInfo table
      if ( strcasecmp( StatInfo->ColNameList[sidx], DSFWInfo->ColNameList[cidx] ) == 0 ) 
	break;
    if ( sidx < StatInfo->Ncols ) 
      // variable found in current table
      DSFWInfo->ColNumList[cidx] = sidx;
    else {
      // variable not found in current table, so check if it is in input file
      for ( idx = 0; idx < NumCols; idx++ )
	if ( strcasecmp( ColNames[idx], DSFWInfo->ColNameList[cidx] ) == 0 )
	  break;
      if ( idx == NumCols ) {
	// handle missing column by turning off calculation
	fprintf( stderr, "Columns required for Chilling Hour calculation not found, so will not be computed.\n" );
	DSFWInfo->Ncols = 0;
	return (FALSE);
      }
      // add missing variable to the StatInfo table (must reallocate memory)
      StatInfo->Ncols++;
      if ( idx > StatInfo->MaxColNum ) StatInfo->MaxColNum = idx;
      StatInfo->ColNameList = (char **)realloc( StatInfo->ColNameList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColNameList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColNameList[StatInfo->Ncols-1], DSFWInfo->ColNameList[cidx] );
      StatInfo->ColStatList = (char **)realloc( StatInfo->ColStatList, StatInfo->Ncols*sizeof(char *) );
      StatInfo->ColStatList[StatInfo->Ncols-1] = (char *)calloc(25,sizeof(char));
      strcpy( StatInfo->ColStatList[StatInfo->Ncols-1], "none" );
      StatInfo->ColNumList = (int *)realloc( StatInfo->ColNumList, StatInfo->Ncols*sizeof(int) );
      StatInfo->ColNumList[StatInfo->Ncols-1] = idx;
      DSFWInfo->ColNumList[cidx] = StatInfo->Ncols-1;
    }
    // reset max column setting for current statistic
    if ( DSFWInfo->ColNumList[cidx] > DSFWInfo->MaxColNum ) 
      DSFWInfo->MaxColNum = DSFWInfo->ColNumList[cidx]+1;
 }
 
  // check is the calculated statistics is in Statinfo more than once
  DSFWInfo->OutputCols = (int *)calloc(StatInfo->Ncols,sizeof(int));
  DSFWInfo->Noutput = 0;
  for ( sidx=0; sidx<StatInfo->Ncols; sidx++ ) {
    if ( strcasecmp( StatInfo->ColNameList[sidx], DSFWInfo->ColNameList[DSFWInfo->Ncols-1] ) == 0 ) {
      // for each occurance of stat name add column to list
      DSFWInfo->OutputCols[DSFWInfo->Noutput] = sidx;
      DSFWInfo->Noutput++;
    }
  }

 return (TRUE);

}

int CalcDaysSuitable4FieldWork( double *DSFW, double *PREC, double *TMAX, double *TMIN, DATE_STRUCT startdate, int Nvals ) {
/**********************************************************************
  CalcDaysSuitable4FieldWork     Keith Cherkauer           June 19, 2017

  Ben Gramig's regression model to determine the workability of 
  agricultural fields.

  Source: Gramig, Massey and Yun (2017). Nitrogen application 
    decision-making under climate risk in the US Corn Belt, 
    Climate Risk Management, 15, 82-89.

  Modifications

**********************************************************************/

  int ErrNum=0;
  int didx, widx;
  double avgTMIN[52], avgTMAX[52], sumPREC[52];
  double tmpStartJulDate, tmpEndJulDate;
  double *tmpTair;

  // Compute days suitable for field work
  DSFW[0] = 0;
  for ( didx=1; didx<Nvals; didx++ ) {
  }
  for ( widx = 1; widx < 52; widx++ ) {
    /* DSFW[didx] = -5.56285 + 0.210292*(avgTMAX) - 0.00237*pow(avgTMAX,2.)  */
    /*   -0.23367*(avgTMIN) - 0.00251*pow(avgTMIN,2)  */
    /*   + 0.005212*((avgTMAX)*(avgTMIN)) - 1.12549*(sumPREC)  */
    /*   + 0.113922*(sumPREC,2) + 0.055012*(avgTMAX[widx-1])  */
    /*   - 0.03594*(avgTMIN[widx-1]) - 1.01346*(sumPREC[widx-1])  */
    /*   + 0.00868*((sumPREC[widx-1])*(avgTMIN[widx]))  */
    /*   + 0.021476*(annual time trend) - 0.00025*pow(annual time trend,2)  */
    /*   - 0.45365*(spring indicator) + 0.374839*(fall indicator)  */
    /*   - 0.42513*(winter indicator) + 3.065181*(soil drainage class index) */
    /*   - 0.43322*pow(soil drainage class index,2); */
  }

  return ErrNum;

}

char CheckExistingOutputFiles ( char              ArcGrid,
				char              StatType,
				char             *GridName,
				StatInfoStruct    StatInfo, 
				char           ***GridFileName,
				DATE_STRUCT       date,
				int               sidx,
				int               cidx,
				int               NumSets,
				int               nrows,
				int               ncols,
				double            ll_lat,
				double            ll_lng,
				double            cellsize,
				int               NODATA,
				int               Ncells, // same as Nfiles for XYZ case
				char              OVERWRITE) {
/**********************************************************************
  CheckExistingOutputFiles     Keith Cherkauer           June 14, 2017

  This subroutine opens an existing data file in ArcInfo Grid or XYZ
  format, extracts it contents, and then checks that the data in the
  file corresponds to that in the previously defined analysis domain.

  Modifications

**********************************************************************/
  int ErrNum;
  char idstr[MaxCharData];

  // build output file name
  if ( sidx < NumSets ) {
    // this is an annual output file
    if ( StatInfo.UseColFile )
      sprintf(idstr, "%s_%s", StatInfo.ColNameList[cidx], 
	      StatInfo.ColStatList[cidx]);
    else 
      sprintf(idstr, "%i", StatInfo.ColNumList[cidx] );	
    if ( StatType == 'A' || StatType == 'P' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, idstr);
    else if ( StatType == 'S' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, SeasonNames[date.season], idstr);
    else if ( StatType == 'M' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, MonthNames[date.month-1], idstr);
    else if ( StatType == 'W' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, date.doy, date.doy+6, idstr);
    else 
      sprintf(GridFileName[sidx][cidx], GridName, date.year, date.doy, idstr);
  }
  else {
    // this is an overall annual average statistic file
    if ( StatInfo.UseColFile )
      sprintf(idstr, "%s_%s", StatInfo.ColNameList[cidx], 
	      StatInfo.ColStatList[cidx]);
    else 
      sprintf(idstr, "%i", StatInfo.ColNumList[cidx] );	
    if ( StatType == 'A' || StatType == 'P' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, idstr);
    else if ( StatType == 'S' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, SeasonNames[date.season], idstr);
    else if ( StatType == 'M' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, MonthNames[date.month-1], idstr);
    else if ( StatType == 'W' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, date.doy, date.doy+6, idstr);
    else 
      sprintf(GridFileName[sidx][cidx], GridName, 9999, date.doy, idstr);
  }

  /* ErrNum = OVERWRITE; */
  /* if ( OVERWRITE ) */
  /*   ErrNum = CheckExistingArcInfoGridFiles ( StatType, GridName,  */
  /* 					     StatInfo, GridFileName,  */
  /* 					     startdate, 0, 0, &nrows,  */
  /* 					     &ncols, &minlat, &minlng,  */
  /* 					     &cellsize, &NODATA, OVERWRITE ); */
  /* else */
  /*   // make sure that row equals number of files, while col is always 1 */
  /*   ErrorNum = CheckExistingXyzFiles ( StatType, GridName, StatInfo,  */
  /* 				       GridFileName, startdate, 0, */
  /* 				       0, &nrows, &ncols, &minlat, */
  /* 				       &minlng, &cellsize, &NODATA,  */
  /* 				       OVERWRITE ); */

  return (OVERWRITE);

}

int read_arcinfo_grid(char     *filename,
		      double  **lat,
		      double  **lng,
		      int      *ncols,
		      int      *nrows,
		      double   *ll_lat,
		      double   *ll_lng,
		      double   *cellsize,
		      int      *NODATA,
		      double ***values) {
/**********************************************************************
  read_arcinfo_info           Keith Cherkauer           May 5, 1998

  This subroutine reads data from an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

  Modifications
  14-Jun-2017 Reversed sense of latitude storage so that north it up
    (first row) rather than down as previously stored.        KAC

**********************************************************************/

  FILE   *farc;
  int     i, j;
  double  tmp_lat;
  double  tmp_lng;
  int     cell;
  int     Ncells;
  double  tmpvalue;

  if((farc=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open ARCINFO grid %s for reading.\n",filename);
    exit(0);
  }

  /***** Read ARC/INFO Header *****/
  fscanf(farc,"%*s %i",  ncols );
  fscanf(farc,"%*s %i",  nrows);
  fscanf(farc,"%*s %lf", ll_lng );
  fscanf(farc,"%*s %lf", ll_lat);
  fscanf(farc,"%*s %lf", cellsize );
  fscanf(farc,"%*s %i",  NODATA );
  Ncells = ncols[0]*nrows[0];

  /***** Allocate memory for data storage *****/
  (*lat) = (double *)calloc((*nrows),sizeof(double));
  (*lng) = (double *)calloc((*ncols),sizeof(double));
  (*values) = (double **)calloc((*nrows),sizeof(double *));
  for ( j = 0; j < (*nrows); j++ )
    (*values)[j] = (double *)calloc((*ncols),sizeof(double));

  /***** Check for Valid Location *****/
  cell = 0;
  for ( j = 0; j < (*nrows); j++ ) {
    tmp_lat = (*ll_lat) + (double)( (*nrows) - j - 0.5 ) * (*cellsize);
    for ( i = 0; i < (*ncols); i++ ) {
      tmp_lng = (*ll_lng) + (double)( i + 0.5 ) * (*cellsize);
      fscanf(farc,"%lf",&tmpvalue);
      //(*lat)[(*nrows)-j-1] = tmp_lat;
      (*lat)[j] = tmp_lat;
      (*lng)[i] = tmp_lng;
      //(*values)[(*nrows)-j-1][i] = tmpvalue;
      (*values)[j][i] = tmpvalue;
      cell++;
    }
  }
  fclose(farc);
  Ncells = cell;

  return Ncells;

}

int write_arcinfo_grid(char   *filename,
		       double **values,
		       double  *lat,
		       double  *lng,
		       int     ncols,
		       int     nrows,
		       double   ll_lat,
		       double   ll_lng,
		       double   cellsize,
		       int     NODATA,
		       int     Ncells) {
/**********************************************************************
  write_arcinfo_info       Keith Cherkauer           April 14, 1999

  This subroutine writes data to an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

  Modifications:
  09-22-04 Modified to make write cell values even if they are not in
           the correct order.                                   KAC
  14-Jun-2017 Reversed sense of latitude storage so that north it up
    (first row) rather than down as previously stored.        KAC

**********************************************************************/

  FILE   *farc;
  int    i, j;
  int    cell;
  int    ErrNum;

  if((farc=fopen(filename,"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open ARCINFO grid %s for writing.\n",filename);
    ErrNum = -1;
    return( ErrNum );
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
  for ( j = 0; j < nrows; j++ ) {
    for ( i = 0; i < ncols; i++ ) {
      //if ( values[nrows-j-1][i] != NODATA )
      //fprintf(farc,"%lf",values[nrows-j-1][i]);
      if ( values[j][i] != NODATA )
	fprintf(farc,"%lf",values[j][i]);
      else
	fprintf(farc,"%i",(int)NODATA);
      cell++;
      if ( i < ncols - 1 ) fprintf(farc,"\t"); // within line
      else fprintf(farc,"\n"); // end of line
    }
  }
  fclose(farc);
  if(Ncells != cell) {
    fprintf(stderr,"WARNING: number of cells written (%i) does not equal number of cells defined (%i).\n", cell, Ncells);
  }

  return cell;

}

int write_blank_arcinfo_grid(char   *filename,
			     int     ncols,
			     int     nrows,
			     double   ll_lat,
			     double   ll_lng,
			     double   cellsize,
			     int     NODATA) {
/**********************************************************************
  write_blank_arcinfo_info       Keith Cherkauer           September 21, 2004

  This subroutine writes an ArcInfo grid filled with NODATA values.

  Modifications:
  14-Jun-2017 Reversed sense of latitude storage so that north it up
    (first row) rather than down as previously stored.        KAC

**********************************************************************/

  FILE   *farc;
  int    i, j;
  int    cell;

  if((farc=fopen(filename,"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open ARCINFO grid %s for writing.\n",filename);
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
  for ( j = 0; j < nrows; j++ ) {
    for ( i = 0; i < ncols; i++ ) {
      fprintf(farc,"%i",NODATA);
      if(i<ncols-1) fprintf(farc,"\t");
      else fprintf(farc,"\n");
    }
  }
  fclose(farc);

  return cell;

}

int WriteOutputFiles ( char              ArcGrid,
		       char             *GridFileName,
		       double          **GridValues,
		       double           *GridLat,
		       double           *GridLng,
		       int               nrows,
		       int               ncols,
		       double            minlat,
		       double            minlng,
		       double            cellsize,
		       int               NODATA,
		       int               Ncells,
		       char              OVERWRITE ) {
/**********************************************************************
  WriteOutputFiles       Keith Cherkauer           June 13, 2017

  This subroutine writes summary statistics to either ARC/INFO ascii 
  grid files, or ASCII XYZ files.  When OVERWRITE flag is set, the 
  program reads an existing file of he same format and uses values from 
  the file read to fill missing values from the current job.

**********************************************************************/
  int ErrNum;

  // Read in existing data and use to complete current data array
  if ( !OVERWRITE ) {
    // here the original file should be read in, and any missing values
    // in the new dataset shoudl be replaced using the values from 
    // the existing file.  In that way, original values are preserved, if
    // nothing new has been estimated.
  }

  // Write existing data array to output file
  if ( ArcGrid )
    // Going to output statistics as an ArcInfo grid
    ErrNum = write_arcinfo_grid(GridFileName, GridValues, 
				GridLat, GridLng, ncols, nrows, 
				minlat, minlng, cellsize, NODATA, Ncells);
  else
    // Output to an ASCII XYZ file
    
    ErrNum = WriteXyzFiles(GridFileName, GridValues, GridLat, GridLng, 
			   nrows, ncols, minlat, minlng, 
			   cellsize, NODATA, Ncells);
  
  return (ErrNum);
}

int WriteXyzFiles ( char    *filename,
		   double **values,
		   double  *lat,
		   double  *lng,
		   int      ncols,
		   int      nrows,
		   double   ll_lat,
		   double   ll_lng,
		   double   cellsize,
		   int      NODATA,
		   int      Ncells ) {
/**********************************************************************
  WriteXyzFiles       Keith Cherkauer           June 13, 2017

  This subroutine writes summary statistics to an ASCII XYZ file.  

  NOTE: XYZ files are stored with nrows = number of XY pairs in the 
    file, and ncols = 1.

**********************************************************************/
  FILE *fxyz;
  int ErrNum=0;
  int row;

  // Open output file
  if ( ( fxyz = fopen( filename, "w" ) ) == NULL ) {
    fprintf(stderr,"ERROR: Unable to open XYZ file %s for writing.\n",filename);
    ErrNum = -1; // cannot open file for some reason, report error
    return(ErrNum);
  }

  // Write data to output file
  for ( row = 0; row < nrows; row++ )
    fprintf( fxyz, "%f %f %f\n", lng[row], lat[row], values[row][0] );

  // close output file
  ErrNum = fclose(fxyz);

  return(ErrNum);

}

int ReadSpatialDataFiles ( char              ArcGrid,
			   char             *GridFileName,
			   double          **GridValues,
			   double           *GridLat,
			   double           *GridLng,
			   int               nrows,
			   int               ncols,
			   double            minlat,
			   double            minlng,
			   double            cellsize,
			   int               NODATA ) {
/**********************************************************************
  ReadSpatialDataFiles       Keith Cherkauer           June 14, 2017

  This subroutine opens an existing data file in ArcInfo Grid or XYZ
  format, extracts it contents, and then checks that the data in the
  file corresponds to that in the previously defined analysis domain.

**********************************************************************/
  int ErrNum, Ncells;
  double *tmp_lat, *tmp_lng;
  int tmp_nrows, tmp_ncols;
  double tmp_ll_lat, tmp_ll_lng, tmp_cellsize;
  int tmp_NODATA;
  double **tmp_values;
  int tmp_row, row, tmp_col, col, ridx, cidx, cell;
  double lat, lng;

  // Read existing data array from output file
  if ( ArcGrid ) {
    /**** Reading from an ArcInfo grid file *****/
    Ncells = read_arcinfo_grid(GridFileName, &tmp_lat, &tmp_lng, &tmp_ncols, 
			       &tmp_nrows, &tmp_ll_lat, &tmp_ll_lng, 
			       &tmp_cellsize, &tmp_NODATA, &tmp_values );
    if ( cellsize != NODATA ) {
      // Match coordinates from new input file with those already defined
      if ( fabs( tmp_cellsize - cellsize ) > SMALL_CHK_VAL ) {
	fprintf( stderr, "ERROR: Cellsize of %s (%f), does not match previously defined cellsize (%f).\n", GridFileName, tmp_cellsize, cellsize );
	exit(FAIL);
      }
      // find input file upper left coordinate in original grid domain
      tmp_col=0;
      col = (int)((tmp_ll_lng + (double)tmp_col*tmp_cellsize - minlng) / cellsize);
      if ( col < 0 ) {
	tmp_col = col * -1;
	col = 0;
	if ( tmp_ncols - tmp_col > ncols ) {
	  fprintf( stderr, "ERROR: Not enough columns in the file being read, %s, to match existing domain.\n", GridFileName );
	  exit(FAIL);
	}
      } 
      tmp_row = 0;
      row = (tmp_ll_lat + (double)tmp_row*tmp_cellsize - minlat) / cellsize;
      if ( row < 0 ) {
	tmp_row = row * -1;
	row = 0;
	if ( tmp_nrows - tmp_row > nrows ) {
	  fprintf( stderr, "ERROR: Not enough rows in the file being read, %s, to match existing domain.\n", GridFileName );
	  exit(FAIL);
	}
      }
    }
    // copy new data into provided array
    cell = 0;
    for ( ridx = 0; ridx < nrows; ridx++ )
      for ( cidx = 0; cidx < ncols; cidx++ ) {
	if ( fabs( tmp_values[ridx+tmp_row][cidx+tmp_col] - tmp_NODATA ) < SMALL_CHK_VAL )
	  // convert to standard no data value
	  GridValues[ridx][cidx] = NODATA;
	else
	  // store current value
	  GridValues[ridx][cidx] = tmp_values[ridx+tmp_row][cidx+tmp_col];
	cell++;
      }
    Ncells = cell;
  }
  else {
    /***** Read from an ASCII XYZ file *****/
    Ncells = ReadXyzFiles(GridFileName, &tmp_lat, &tmp_lng, &tmp_ncols, 
			  &tmp_nrows, &tmp_ll_lat, &tmp_ll_lng, 
			  &tmp_cellsize, &tmp_NODATA, &tmp_values );
    // Match coordinates from new input file with those already defined
    cell = 0;
    for ( tmp_row = 0; tmp_row < tmp_nrows; tmp_row++ ) {
      for ( row = 0; row < nrows; row++ ) {
	if ( ( fabs( GridLat[row] - tmp_lat[tmp_row] ) < SMALL_CHK_VAL ) &&
	     ( fabs( GridLng[row] - tmp_lng[tmp_row] ) < SMALL_CHK_VAL ) ) {
	  // found matching location, copy into provided data array
	  GridValues[row][0] = tmp_values[tmp_row][0];
	  cell++;
	  break;
	}
      }
    }
    // Check that the correct number of cells were found
    if ( cell < nrows ) {
      fprintf( stderr, "ERROR: Too few lines (%i < %i) in XYZ file %s, cannot match new values to existing cells, so program is exiting.\n", cell, nrows, GridFileName );
      exit(-1);
    }
    Ncells = cell;
  }

  /***** free memory allocated for this subroutine *****/
  free (tmp_lat);
  free (tmp_lng);
  for ( row = 0; row < tmp_nrows; row++ ) free (tmp_values[row]);
  free(tmp_values);

  return (Ncells);

}

int ReadXyzFiles ( char     *filename,
		   double  **lat,
		   double  **lng,
		   int      *ncols,
		   int      *nrows,
		   double   *ll_lat,
		   double   *ll_lng,
		   double   *cellsize,
		   int      *NODATA,
		   double ***values ) {
/**********************************************************************
  ReadXyzFiles           Keith Cherkauer           June 14, 2017

  This subroutine reads data from an XYZ ascii file and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cells.  Data are stored to a two-dimensional array where first
  index [row] is equal to the number of grid cells, and the second index
  [col] is set to 1.  This maintains compatability with the output from
  read_arcinfo_grid.

**********************************************************************/

  FILE   *fxyz;
  int     row, col, i, j;
  double  tmp_lat;
  double  tmp_lng;
  int     cell;
  int     Ncells;
  double  tmpvalue;
  char    tmpstr[MaxCharData];

  // Open input file
  if ( ( fxyz = fopen( filename, "r" ) ) == NULL ) {
    fprintf(stderr,"ERROR: Unable to open XYZ file %s for reading.\n",filename);
    return(FAIL);
  }

  // Count number of lines in the input file
  Ncells=0;
  while ( fgets(tmpstr, MaxCharData, fxyz) != NULL ) Ncells++;
  rewind(fxyz);
  (*nrows) = Ncells;
  (*ncols) = 1;

  // Allocate memory for data storage
  (*lat) = (double *)calloc(Ncells,sizeof(double));
  (*lng) = (double *)calloc(Ncells,sizeof(double));
  (*values) = (double **)calloc(Ncells,sizeof(double *));
  for ( j = 0; j < (*nrows); j++ )
    (*values)[j] = (double *)calloc(1,sizeof(double));

  // Write data to output file
  for ( row = 0; row < Ncells; row++ )
    fscanf( fxyz, "%f %f %f\n", &lng[row], &lat[row], &values[row][0] );

  // close output file
  fclose(fxyz);

  return (Ncells);

}


