#include <VicStuff.h>

int GetStatInfo( char *, StatInfoStruct *, char **, int );
int GetPenmanInfo( PenInfoStruct *, char **, int );
int GetTotalRunoffInfo( PenInfoStruct *, char **, int );
DATE_STRUCT GetJulianDate(int, int, int, double);
int CreateNewBlankArcInfoGridFiles ( char, char *, StatInfoStruct, char ***, 
				     DATE_STRUCT, int, int, int, int, int, double, 
				     double, double, int, char );
int CreateNewBlankXyzFiles ( char, char *, StatInfoStruct, char ***, 
			     DATE_STRUCT, int, int, int, int, int, double, 
			     double, double, int, char );
char CheckExistingArcInfoGridFiles ( char, char *, StatInfoStruct, char ***, 
				     DATE_STRUCT, int, int, int *, int *, 
				     double *, double *, double *, int * );
char CheckExistingXyzFiles ( char, char *, StatInfoStruct, char ***, DATE_STRUCT,
			     int, int, int *, int *, double *, double *, double *,
			     int * );
int read_arcinfo_grid(char *, double *, double *, int *, int *, double *, 
		      double *, double *, int *, double **);
int write_arcinfo_grid(char *, double **, double *, double *, int, int, double, 
		       double, double, int, int);
int write_blank_arcinfo_grid(char *, int, int, double, double, double, int);

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

***********************************************************************/

  StatInfoStruct StatInfo;
  PenInfoStruct PenInfo, TRoffInfo;

  gzFile *fin;
  FILE *flist;

  char ***GridFileName;
  char filename[512];
  char GridName[512];
  char StatType;
  char OVERWRITE;
  char PrtAllGrids = FALSE;
  char ArcGrid = FALSE;
  char ProcessFile;

  int filenum;
  int Nfile;
  int row;
  int nrow;
  int col;
  int ncol;
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
  int     CalcPE, CalcTR;

  date      = (int *)malloc(NUM_DATE_VALS*sizeof(int));

  /** Process Command Line Arguments **/
  if ( argc != 11 ) {
    fprintf (stderr, "NumVars = %i\n", argc );
    fprintf(stderr,"\nUsage: %s <file list file> <Output ArcGrid: TRUE/FALSE> <output prefix> <grid resolution> <column list file> <start date> <end date> <V|A|S|M|W|D> <OVERWRITE: TRUE/FALSE> <PrtAllPeriods: TRUE/FALSE>\n",argv[0]);
    fprintf(stderr,"\n\tThis program produces either an ARC/INFO ASCII grid file \n\t(output ArcGrid = T) or an XYZ file (Output ArcGrid = F) of the average\n\tof the selected data column for the given time period.\n");
    fprintf(stderr,"\n\t<file list> is a file containing the full grid file name and location, \n\t\tlatitude and longitude of the grid cell, for each grid cell to be\n\t\tincluded.\n");
    fprintf(stderr,"\t<output prefix> is the prefix (path and start of file name) for the \n\t\toutput files that will be generated by this program.  A suffix \n\t\twill be added to all file names to separate individual output \n\t\tfor each variable, and for each tperiod (year, season, month, \n\t\tetc). Files containing multi-year average statistics for annual, \n\t\tseasonal and monthly periods will use \"9999\" for \n\t\tthe date in the file name.\n");
    fprintf(stderr,"\t<grid resolution> is the resolution in degrees of the desired output\n\t\tgrid.\n");
    fprintf(stderr,"\t<column list file> is a multi-column ASCII file that lists the column\n\t\tname for each column to be output.  Followed by the statistic \n\t\tto be computed and a thresehold if required by the statistic.  \n\n\t\tStatistic options include: \n\t\t\tMean value \t\t\'mean\', \n\t\t\tCumulative value \t\'sum\', \n\t\t\tStandard Deviation \t\'stdev\', \n\t\t\tMaximum value \t\t\'max\', \n\t\t\tMinimum value \t\t\'min\', \n\t\t\tFirst value \t\t\'first\', \n\t\t\tLast value \t\t\'last\', \n\t\t\tLast day over thres. \t\'ldayo\' \t<thres>, \n\t\t\tLast day under thres. \t\'ldayu\' \t<thres>, \n\t\t\tFirst day over thres. \t\'fdayo\' \t<thres>, \n\t\t\tFirst day under thres. \t\'fdayu\' \t<thres>, \n\t\t\tDays over threshold \t\'othres\' \t<thres>, and \n\t\t\tDays under threshold \t\'uthres\' \t<thres>.\n\n\t\tNOTE: To output potential evapotranspiration, include \"OUT_PE\" \n\t\tas a variable.  To output total runoff (runoff + baseflow), \n\t\tinclude \"OUT_TOTAL_RUNOFF\" as a variable.\n\n");
    fprintf(stderr,"\t<start date> and <end date> are the starting and ending dates of the\n\t\tperiod of interest in MMDDYYYY format (MM = month, DD = day,\n\t\tYYYY = year - date must be 8 characters).\n");
    fprintf(stderr,"\t<V|A|S|M|W|D> export Annual a(V)erage, (A)nnual, (S)easonal, (M)onthly, \n\t\t(W)eekly or (D)aily grids for all years in the file\n\t\t- Annual uses given start date to start year; \n\t\t- Seasonal uses Winter = DJF, Spring = MAM, Summer = JJA, and \n\t\t  Autumn = SON; \n\t\t- Weekly parses data into 7 day weeks starting with the given \n\t\t  start date.\n");
    fprintf(stderr,"\t<OVERWRITE> if set to TRUE then the output grid files will be \n\t\toverwritten.  The default is to replace values in existing \n\t\tfiles with new data.\n");
    fprintf(stderr,"\t<PrtAllPeriods> if set to TRUE then the program will write all output \n\t\tfiles (all years, all seasons, etc), if set to FALSE only \n\t\tannual averages will be written.\n\n");
    exit(0);
  }

  // Read grid cell file list and create template for ArcInfo grid file
  if((flist=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open file list %s.\n",argv[1]);
    exit(0);
  }
  Nfile=0;
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

  if ( strcasecmp( argv[2], "TRUE" ) == 0 ) ArcGrid = TRUE;
  else ArcGrid = FALSE;

  if ( ( cellsize = atof(argv[4]) ) <= 0 ) {
    fprintf(stderr, "ERROR: Cell size (%f) must be a positive doubling value!\n\n", cellsize);
    return(-1);
  }

  maxlat += cellsize / 2.;
  minlat -= cellsize / 2.;
  maxlng += cellsize / 2.;
  minlng -= cellsize / 2.;

  nrow = (int)((maxlat - minlat) / cellsize);
  ncol = (int)((maxlng - minlng) / cellsize);
  Ncells = nrow*ncol;
  if ( nrow == 1 || ncol == 1 ) {
    fprintf( stderr, "WARNING: Defined cells are at best one-dimensions with %i rows and %i columns.  If you are expecting a raster output, then double check your file list format.\n", nrow, ncol );
  }

  // Setup input and output directories
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
  ColNames = (char **)realloc(ColNames,(NumCols+2)*sizeof(char *));
  // set OUT_PE as a column name
  ColNames[NumCols] = (char *)calloc( 7, sizeof(char) );
  strcpy( ColNames[NumCols], "OUT_PE" );
  // set OUT_TOTAL_RUNOFF as a column name
  ColNames[NumCols+1] = (char *)calloc( 17, sizeof(char) );
  strcpy( ColNames[NumCols+1], "OUT_TOTAL_RUNOFF" );
  // check for hourly column in dataset
  if ( strcmp( ColNames[3], "HOUR" ) == 0 ) { // hour is present offset by one column
    OutStart = 4;
    VarList = &ColNames[OutStart];
    NumOut  = NumCols - OutStart + 2;
  }
  else { // hour is not present, do not offset daily data
    OutStart = 3;
    VarList = &ColNames[OutStart];
    NumOut = NumCols - OutStart + 2;
  }
  data = (double *)malloc((NumOut+2)*sizeof(double));
  gzclose(fin);

  /** Read statistics control file **/
  ErrNum = GetStatInfo( argv[5], &StatInfo, VarList, NumOut );

  /** Setup for calculation of PE **/
  CalcPE = GetPenmanInfo( &PenInfo, VarList, NumOut );

  /** Setup for calculation of total runoff (runoff + baseflow) **/
  CalcTR = GetTotalRunoffInfo( &TRoffInfo, VarList, NumOut );

  /** Check that requests columns are in the given output files **/
  for ( didx=NumOut-1; didx>=0; didx-- ) {
    for( cidx=0; cidx<NumCols+2; cidx++ ) {
      if ( strcasecmp( ColNames[cidx], VarList[didx] ) == 0 ) break;
    }
    if ( cidx >= NumCols+2 ) {
      fprintf( stderr, "WARNING: Unable to find output column %s in the given output files, will remove from analysis, check the VIC model output file header to determine what variables are in the file %s.\n", VarList[didx], filename );
      for( cidx=didx; cidx<NumOut-1; cidx++ )
	VarList[cidx] = VarList[cidx+1]; // remove missing statistic
      NumOut--;
    }
  }

  // Initialize dates
  startdate.month   = (int)(atof(argv[6])/1000000);
  startdate.day     = (int)(atof(argv[6])/10000) - startdate.month*100;
  startdate.year    = (int)(atof(argv[6])) - startdate.month*1000000 - startdate.day*10000;
  startdate.hour    = 0;
  startdate.juldate = calc_juldate(startdate.year,startdate.month,startdate.day,0.);
  startdate         = get_season(startdate);

  enddate.month     = (int)(atof(argv[7])/1000000);
  enddate.day       = (int)(atof(argv[7])/10000) - enddate.month*100;
  enddate.year      = (int)(atof(argv[7])) - enddate.month*1000000 - enddate.day*10000;
  enddate.hour      = 0;
  enddate.juldate   = calc_juldate(enddate.year,enddate.month,enddate.day,0);
  enddate           = get_season(enddate);
  
  // Process statistic type - adjust start and end date to encompass full periods
  StatType    = argv[8][0];
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
  else {
    // Daily Statistics
    NumSets = (int)(enddate.juldate - startdate.juldate);
    SumSets = 0;
    strcat( GridName, "_\%i_\%u_\%s.asc" ); // build daily output file name
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

  // Initialize Data Array
  gridval = (double **)malloc(StatInfo.Ncols*sizeof(double *));  
  for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
    gridval[cidx] = (double *)malloc(366*sizeof(double));  
  }
  GridFileName = (char ***)malloc((NumSets+SumSets)*sizeof(char **));
  for ( sidx = 0; sidx < NumSets+SumSets; sidx++ ) {
    GridFileName[sidx] = (char **)malloc(StatInfo.Ncols*sizeof(char *));
    for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) 
      GridFileName[sidx][cidx] = (char *)malloc(1024*sizeof(char));
  }

  /***************************************************
    If not overwriting file validate previous files
  ***************************************************/
  if ( !OVERWRITE ) 
    if ( ArcGrid )
      OVERWRITE = CheckExistingArcInfoGridFiles ( StatType, GridName, StatInfo, 
						  GridFileName, startdate, 0,
						  0, &nrow, &ncol, &minlat,
						  &minlng, &cellsize, &NODATA );
    else
      // make sure that row equals number of files, while col is always 1
      OVERWRITE = CheckExistingXyzFiles ( StatType, GridName, StatInfo, 
					  GridFileName, startdate, 0,
					  0, &nrow, &ncol, &minlat,
					  &minlng, &cellsize, &NODATA );

  GridLat = (double *)calloc(nrow, sizeof(double));
  GridLng = (double *)calloc(ncol, sizeof(double));
 
  GridValues = (double ****)malloc(NumSets*sizeof(double ***)); 
  for( sidx =0; sidx<NumSets; sidx++) {
    GridValues[sidx] = (double ***)malloc(StatInfo.Ncols*sizeof(double **));
    for( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
      GridValues[sidx][cidx] = (double **)malloc(nrow*sizeof(double *));
      OutGrid = (double **)malloc(nrow*sizeof(double *));
      for ( row = 0; row < nrow; row++ ) {
	GridValues[sidx][cidx][row] = (double *)malloc(ncol*sizeof(double));
	OutGrid[row] = (double *)malloc(ncol*sizeof(double));
	for ( col = 0; col < ncol; col++ )
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

  // Build blank files for overall statistics
  for ( sidx = 0; sidx < SumSets; sidx++ ) {
    for( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
      if ( ArcGrid )
	// Going to output statistics as an ArcInfo grid
	ErrNum = CreateNewBlankArcInfoGridFiles(StatType, GridName,
						StatInfo, GridFileName,
						nextdate, NumSets+sidx, cidx, 
						NumSets, nrow, ncol, 
						minlat, minlng,
						cellsize, NODATA, 
						OVERWRITE);
      else
	// Output to an ASCII XYZ file
	ErrNum = CreateNewBlankXyzFiles(StatType, GridName,
					StatInfo, GridFileName,
					nextdate, NumSets+sidx, cidx, 
					NumSets, nrow, ncol, 
					minlat, minlng,
					cellsize, NODATA, 
					OVERWRITE);
    }
    if ( StatType == 'A' ) nextdate = get_next_year( nextdate );
    else if ( StatType == 'S' ) nextdate = get_next_season( nextdate );
    else if ( StatType == 'M' ) nextdate = get_next_month( nextdate );
    else if ( StatType == 'W' ) nextdate = get_next_week( nextdate );
    else if ( StatType == 'V' ) nextdate = enddate;
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
    
    /**** Need to determine row and column of current file in grid file ****/
    if ( ArcGrid ) {
      DONE = 0;
      row = 0;
      col = 0;
      findlat = minlat + (row + 0.5) * cellsize;
      findlng = minlng + (col + 0.5) * cellsize;
      
      // Determine the row position of the current cell
      while(findlat + cellsize*0.25 < tmplat) {
	row++;
	findlat = minlat + (row + 0.5) * cellsize;
      }
      
      // Determine column position of current cell
      while ( findlng + cellsize*0.25 < tmplng ) {
	col++;
	findlng = minlng + (col + 0.5) * cellsize;
      }

      // Check that current row and column are in the Arc/Info file boundaries
      if ( row >= nrow || col >= ncol ) {
	fprintf(stderr,"WARNING: %s located at (%i,%i) is outside the range of the grid file (%i,%i). Will not be processed.\n",
		filename,row,col,nrow,ncol);
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

      if ( filenum == 0 ) {
	// Use first file to setup all additional processing
	// determine number of bytes in full flux file, in case user has requested shorter time period
	if ( BinaryFile ) {
	  // determine number of bytes in full binary flux file, in case user has requested shorter time period
	  RawData = (char **)malloc(sizeof(char*));
	  RawData[0] = (char *)malloc(NumBytes*sizeof(char));
	  while ( ( ReadBytes = gzread(fin,RawData[0],NumBytes*sizeof(char)) ) > 0 ) {
	    NumRead += ReadBytes;
	  }
	}
	else {
	  // determine number of lines in full ASCII flux file, in case shorter time period requested
	  RawData = (char **)malloc(sizeof(char *));
	  RawData[0] = (char *)malloc(MaxCharData*sizeof(char));
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
	  if ( ( RawData = (char **)malloc( NumRead*sizeof(char * ) ) ) == NULL ) {
	    fprintf( stderr, "ERROR: Unable to allocate sufficient memory to read in full flux file.\n" );
	    return (-1);
	  }
	  for ( lidx = 0; lidx < NumRead; lidx ++ ) {
	    if ( ( RawData[lidx] = (char *)malloc( MaxCharData*sizeof(char) ) ) == NULL ) {
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
			       CalcPE, PenInfo, CalcTR, TRoffInfo );
      
      // reset to correct ASCII position before processing
      if ( !BinaryFile ) RawPtr = 0; 
      
      // compute julian day of current record
      tmpdate.year = date[0];
      tmpdate.month = date[1];
      tmpdate.day = date[2];
      tmpdate.juldate = calc_juldate(tmpdate.year,tmpdate.month,tmpdate.day,0);
      
      while( !DONE && ErrNum >= 0 ) { 
	// process all lines in the current file
	
	if ( tmpdate.juldate >= startdate.juldate && sidx == NODATA ) {
	  // Initialize date accounting for new grid cell
	  sidx = 0;
	  lastdate.juldate = startdate.juldate;
	  if ( StatType == 'A' ) nextdate = get_next_year( startdate );
	  else if ( StatType == 'S' ) nextdate = get_next_season( startdate );
	  else if ( StatType == 'M' ) nextdate = get_next_month( startdate );
	  else if ( StatType == 'W' ) nextdate = get_next_week( startdate );
	  else if ( StatType == 'V' ) nextdate = enddate;
	  else nextdate = get_next_day( startdate );
	  nextdate = get_juldate( nextdate );
	  // Create output file names and initialize grid files
	  for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
	    if ( filenum == 0 && PrtAllGrids ) {
	      ErrNum = CreateNewBlankArcInfoGridFiles(StatType, GridName,
						      StatInfo, GridFileName,
						      startdate, sidx, cidx, 
						      NumSets, nrow, ncol, 
						      minlat, minlng,
						      cellsize, NODATA, 
						      OVERWRITE);
	    }
	    // Initialize daily statistic storage variables for each new grid cell
	    for ( didx = 0; didx < 366; didx++ ) {
	      gridval[cidx][didx] = (double)NODATA;  
	    }
	  }
	  Nrec = 0;
	}
	else if ( tmpdate.juldate >= nextdate.juldate && sidx >= 0 ) {
	  Nvals = (int)(nextdate.juldate+0.5)-(int)(lastdate.juldate+0.5);
	  for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ ) {
	    // compute base statistics
	    if ( strcasecmp(StatInfo.ColStatList[cidx], "max") == 0 )
	      // Find the maximum value
	      GridValues[sidx][cidx][row][col] = get_max(gridval[cidx], Nvals, (double)NODATA );
	    else if ( strcasecmp(StatInfo.ColStatList[cidx], "min") == 0 ) {
	      // Find the minimum value
	      GridValues[sidx][cidx][row][col] = get_min(gridval[cidx], Nvals, (double)NODATA );
	    }
	    else if ( strcasecmp( StatInfo.ColStatList[cidx], "othres" ) == 0 ) {
	      // Count all days a value is under the given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_count_over_thres(gridval[cidx], StatInfo.Thres[cidx], Nvals, (double)NODATA);
	    }
	    else if ( strcasecmp( StatInfo.ColStatList[cidx], "uthres" ) == 0 ) {
	      // Count all days a value is under the given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_count_under_thres(gridval[cidx], StatInfo.Thres[cidx], Nvals, (double)NODATA);
	    }
	    else if ( strcasecmp( StatInfo.ColStatList[cidx], "ldayo" ) == 0 ) {
	      // find the last day a value is over a given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_last_day_over_thres(gridval[cidx], StatInfo.Thres[cidx], Nvals, (double)NODATA);
	      //
	      // SHOULD I ADD DATES FROM January 1 to get DOY?  Do NOT SEEM TO HAVE SUCH A FUNCTION
	      //
	    }
	    else if ( strcasecmp( StatInfo.ColStatList[cidx], "ldayu" ) == 0 ) {
	      // find the last day a value is under a given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_last_day_under_thres(gridval[cidx], StatInfo.Thres[cidx], Nvals, (double)NODATA);
	    }
	    else if ( strcasecmp( StatInfo.ColStatList[cidx], "fdayo" ) == 0 ) {
	      // find the first day a value is over a given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_first_day_over_thres(gridval[cidx], StatInfo.Thres[cidx], Nvals, (double)NODATA);
	    }
	    else if ( strcasecmp( StatInfo.ColStatList[cidx], "fdayu" ) == 0 ) {
	      // find the first day a value is under a given threshold
	      GridValues[sidx][cidx][row][col] = (double)get_first_day_under_thres(gridval[cidx], StatInfo.Thres[cidx], Nvals, (double)NODATA);
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
	    else
	    // Find the mean
	      GridValues[sidx][cidx][row][col] = get_mean(gridval[cidx], Nvals, (double)NODATA );
	    // Create and initialize output file name for next increment
	    if ( filenum == 0 && tmpdate.juldate < enddate.juldate && PrtAllGrids ) {
	      ErrNum = CreateNewBlankArcInfoGridFiles(StatType, GridName,
						      StatInfo, GridFileName,
						      nextdate, sidx+1, cidx, 
						      NumSets, nrow, ncol, minlat, minlng,
						      cellsize, NODATA, 
						      OVERWRITE);
	    }
	    // Initialize daily statistic storage variables for each new grid cell
	    for ( didx = 0; didx < 366; didx++ ) {
	      gridval[cidx][didx] = (double)NODATA;  
	    }
	  }
	  Nrec = 0;
	  // Update Date Information
	  sidx++;
	  lastdate.juldate = nextdate.juldate;
	  if ( StatType == 'A' ) nextdate = get_next_year( nextdate );
	  else if ( StatType == 'S' ) nextdate = get_next_season( nextdate );
	  else if ( StatType == 'M' ) nextdate = get_next_month( nextdate );
	  else if ( StatType == 'W' ) nextdate = get_next_week( nextdate );
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
				   TRoffInfo );
	  
	  if ( ErrNum >= 0 ) {
	    // compute julian day of current record
	    tmpdate.year = date[0];
	    tmpdate.month = date[1];
	    tmpdate.day = date[2];
	    tmpdate.juldate = calc_juldate(tmpdate.year,tmpdate.month,tmpdate.day,0);
	  }
	  else if ( tmpdate.juldate+1 == enddate.juldate ) {
	    tmpdate.year = enddate.year;
	    tmpdate.month = enddate.month;
	    tmpdate.day = enddate.day;
	    tmpdate.juldate = calc_juldate(tmpdate.year,tmpdate.month,tmpdate.day,0);
	    ErrNum = 0;
	  }	

	  // Check that current cell had valid data
	  CheckSum = 0;
	  for ( cidx = 0; cidx < StatInfo.Ncols; cidx++ )
	    if ( gridval[cidx][0] != (double)NODATA )
	      CheckSum += gridval[cidx][0];
	  if ( CheckSum == 0 ) 
	    fprintf(stderr,"WARNING: Possible problem with cell at %f, %f, since summary statistics are all 0 or NODATA.\n", tmplat, tmplng);

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

    // Plot individual time period grids
    if ( PrtAllGrids ) {
      for ( sidx = 0; sidx < NumSets; sidx++ )
	ErrNum = write_arcinfo_grid(GridFileName[sidx][cidx], 
				    GridValues[sidx][cidx], GridLat, GridLng, 
				    ncol, nrow, minlat, minlng, cellsize, 
				    NODATA, Ncells);
    }

    // Compute average for all time periods
    for ( sumidx = 0; sumidx < SumSets; sumidx++ ) {
      for ( row = 0; row < nrow; row++ ) {
	for ( col = 0; col < ncol; col++ ) {
	  if ( GridValues[0][cidx][row][col] == NODATA ) {
	    // no values present, set to no data value
	    OutGrid[row][col] = NODATA;
	  }
	  else {
	    // compute average value
	    OutGrid[row][col] = ValCnt = 0;
	    for ( sidx = sumidx; sidx < NumSets; sidx+=SumSets ) {
	      OutGrid[row][col] += GridValues[sidx][cidx][row][col];
	      ValCnt++;
	    }
	    OutGrid[row][col] /= (double)ValCnt;
	  }
	}
      }

      // Update Grid file for current value
      ErrNum = write_arcinfo_grid(GridFileName[NumSets+sumidx][cidx], OutGrid, 
				  GridLat, GridLng, ncol, nrow, minlat, 
				  minlng, cellsize, NODATA, Ncells);
    }
  }

  return(0);

}

int GetStatInfo( char *StatInfoParam, StatInfoStruct *StatInfo, char **ColNames, int NumCols ) {

  // Modified 2008-May-13 to remove column number from input file.

  FILE *fin;
  int cidx, idx, tidx, Ncols;
  char tmpstr[1024];

  // Open column statistics definition file
  if((fin=fopen(StatInfoParam,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open column statistics definition file %s for reading.\n",StatInfoParam);
    exit(0);
  }

  // Determine columns to be processed
  StatInfo->Ncols = 0;
  fgets(tmpstr, 1024, fin);
  while ( !feof(fin) ) {
    StatInfo->Ncols++;
    fgets(tmpstr, 1024, fin);
  }
  rewind(fin);
  StatInfo->ColNumList  = (int *)calloc(StatInfo->Ncols,sizeof(int));
  StatInfo->ColNameList = (char **)calloc(StatInfo->Ncols,sizeof(char *));
  StatInfo->ColStatList = (char **)calloc(StatInfo->Ncols,sizeof(char *));
  StatInfo->Thres       = (float *)calloc(StatInfo->Ncols,sizeof(float));
  StatInfo->MaxColNum   = -9;
  Ncols = 0;
  for ( cidx = StatInfo->Ncols-1; cidx >= 0; cidx-- ) {
    StatInfo->ColNameList[cidx] = (char *)calloc(25,sizeof(char));
    StatInfo->ColStatList[cidx] = (char *)calloc(7,sizeof(char));
    fgets(tmpstr, 1024, fin);
    sscanf( tmpstr, "%s %s", StatInfo->ColNameList[cidx], 
	    StatInfo->ColStatList[cidx] );
    if ( strcasecmp( StatInfo->ColStatList[cidx], "othres" ) == 0 
	 || strcasecmp( StatInfo->ColStatList[cidx], "uthres" ) == 0 
	 || strcasecmp( StatInfo->ColStatList[cidx], "lday" ) == 0 
	 || strcasecmp( StatInfo->ColStatList[cidx], "fday" ) == 0 ) {
      sscanf( tmpstr, "%*i %*s %*s %f", &StatInfo->Thres[cidx] ); 
    }
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
	StatInfo->Thres[tidx] = StatInfo->Thres[tidx+1];
      }
      StatInfo->Ncols--;
    }
    else if ( StatInfo->ColNumList[cidx] > StatInfo->MaxColNum ) 
      StatInfo->MaxColNum = StatInfo->ColNumList[cidx]+1;
  }
  fclose(fin);
  StatInfo->UseColFile = TRUE;
  StatInfo->Ncols -= Ncols;

  return(0);
}

int GetPenmanInfo( PenInfoStruct *PenInfo, char **ColNames, int NumCols ) {

  // Modified 2008-May-14 copied from GetStatInfo.
  // Modified 2008-May-28 to remove need for separate file to define penman variables.
  // Modified 2016-Jan-12 to add OUT_PE as a variable column name

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

int GetTotalRunoffInfo( PenInfoStruct *TRoffInfo, char **ColNames, int NumCols ) {

  // Modified 2008-Nov-12 copied from GetPenInfo.

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

int CreateNewBlankArcInfoGridFiles ( char StatType,
				     char *GridName,
				     StatInfoStruct StatInfo, 
				     char ***GridFileName,
				     DATE_STRUCT date,
				     int sidx,
				     int cidx,
				     int NumSets,
				     int nrow,
				     int ncol,
				     double minlat,
				     double minlng,
				     double cellsize,
				     int NODATA,
				     char OVERWRITE ) {

  FILE *fin;
  char idstr[250];
  int ErrNum;

  if ( sidx < NumSets ) {
    if ( StatInfo.UseColFile )
      sprintf(idstr, "%s_%s", StatInfo.ColNameList[cidx], 
	      StatInfo.ColStatList[cidx]);
    else 
      sprintf(idstr, "%i", StatInfo.ColNumList[cidx] );	
    if ( StatType == 'A' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      idstr);
    else if ( StatType == 'S' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      SeasonNames[date.season], idstr);
    else if ( StatType == 'M' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      MonthNames[date.month-1], idstr);
    else if ( StatType == 'W' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      date.doy, date.doy+6, idstr);
    else 
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      date.doy, idstr);
  }
  else {
    if ( StatInfo.UseColFile )
      sprintf(idstr, "%s_%s", StatInfo.ColNameList[cidx], 
	      StatInfo.ColStatList[cidx]);
    else 
      sprintf(idstr, "%i", StatInfo.ColNumList[cidx] );	
    if ( StatType == 'A' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      idstr);
    else if ( StatType == 'S' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      SeasonNames[date.season], idstr);
    else if ( StatType == 'M' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      MonthNames[date.month-1], idstr);
    else if ( StatType == 'W' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      date.doy, date.doy+6, idstr);
    else 
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      date.doy, idstr);
  }

  // Create blank file
  if ( OVERWRITE ) 
    ErrNum = write_blank_arcinfo_grid(GridFileName[sidx][cidx], ncol, nrow, 
				      minlat, minlng, cellsize, NODATA);
  else if ( ( fin = fopen(GridFileName[sidx][cidx], "r") ) == NULL )
    ErrNum = write_blank_arcinfo_grid(GridFileName[sidx][cidx], ncol, nrow, 
				      minlat, minlng, cellsize, NODATA);
  else {
    fclose ( fin ); 
    ErrNum = 0;
  }

  return(ErrNum);

}

char CheckExistingArcInfoGridFiles ( char              StatType,
				     char             *GridName,
				     StatInfoStruct    StatInfo, 
				     char           ***GridFileName,
				     DATE_STRUCT        date,
				     int               sidx,
				     int               cidx,
				     int              *nrows,
				     int              *ncols,
				     double            *ll_lat,
				     double            *ll_lng,
				     double            *cellsize,
				     int              *NODATA ) {

  FILE *fin;
  char idstr[250];

  if ( StatInfo.UseColFile )
    sprintf(idstr, "%s_%s", StatInfo.ColNameList[cidx], 
	    StatInfo.ColStatList[cidx]);
  else 
    sprintf(idstr, "%i", StatInfo.ColNumList[cidx] );	
  if ( StatType == 'A' )
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    idstr);
  else if ( StatType == 'S' )
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    SeasonNames[date.season], idstr);
  else if ( StatType == 'M' )
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    MonthNames[date.month-1], idstr);
  else if ( StatType == 'W' )
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    date.doy, date.doy+6, idstr);
  else 
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    date.doy, idstr);

  // Create blank file
  if ( ( fin = fopen(GridFileName[sidx][cidx], "r") ) == NULL )
    return (TRUE);
  else {
    fscanf(fin,"%*s %i",&ncols[0]);
    fscanf(fin,"%*s %i",&nrows[0]);
    fscanf(fin,"%*s %lf",&ll_lng[0]);
    fscanf(fin,"%*s %lf",&ll_lat[0]);
    fscanf(fin,"%*s %lf",&cellsize[0]);
    fscanf(fin,"%*s %i",&NODATA[0]);
    fclose ( fin ); 
    return (FALSE);
  }

  return(-1);

}

int read_arcinfo_grid(char    *filename,
		      double  *lat,
		      double  *lng,
		      int     *ncols,
		      int     *nrows,
		      double  *ll_lat,
		      double  *ll_lng,
		      double  *cellsize,
		      int    *NODATA,
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
  double  tmp_lat;
  double  tmp_lng;
  int    cell;
  int    Ncells;
  double  tmpvalue;

  if((farc=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open ARCINFO grid %s for reading.\n",filename);
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

  /***** Check for Valid Location *****/
  cell = 0;
  for ( j = 0; j < (*nrows); j++ ) {
    tmp_lat = (*ll_lat) + (double)( (*nrows) - j - 0.5 ) * (*cellsize);
    for ( i = 0; i < (*ncols); i++ ) {
      tmp_lng = (*ll_lng) + (double)( i + 0.5 ) * (*cellsize);
      fscanf(farc,"%lf",&tmpvalue);
      lat[(*nrows)-j-1]    = tmp_lat;
      lng[i]    = tmp_lng;
      values[(*nrows)-j-1][i] = tmpvalue;
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

  This subroutine reads data from an ARC/INFO ascii grid and returns 
  the central latitude and longitude of all grid cells that are not 
  assigned the value of NODATA.  Routine returns the number of defined 
  grid cell.

  09-22-04 Modified to make write cell values even if they are not in
           the correct order.                                   KAC

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
      if ( values[nrows-j-1][i] != NODATA )
	fprintf(farc,"%lf",values[nrows-j-1][i]);
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
  write_arcinfo_info       Keith Cherkauer           September 21, 2004

  This subroutine writes an ArcInfo grid filled with NODATA values.

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

int CreateNewBlankXyzFiles ( char StatType,
			     char *GridName,
			     StatInfoStruct StatInfo, 
			     char ***GridFileName,
			     DATE_STRUCT date,
			     int sidx,
			     int cidx,
			     int NumSets,
			     int nrow,
			     int ncol,
			     double minlat,
			     double minlng,
			     double cellsize,
			     int NODATA,
			     char OVERWRITE ) {

  FILE *fin;
  char idstr[250];
  int ErrNum;

  if ( sidx < NumSets ) {
    if ( StatInfo.UseColFile )
      sprintf(idstr, "%s_%s", StatInfo.ColNameList[cidx], 
	      StatInfo.ColStatList[cidx]);
    else 
      sprintf(idstr, "%i", StatInfo.ColNumList[cidx] );	
    if ( StatType == 'A' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      idstr);
    else if ( StatType == 'S' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      SeasonNames[date.season], idstr);
    else if ( StatType == 'M' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      MonthNames[date.month-1], idstr);
    else if ( StatType == 'W' )
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      date.doy, date.doy+6, idstr);
    else 
      sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	      date.doy, idstr);
  }
  else {
    if ( StatInfo.UseColFile )
      sprintf(idstr, "%s_%s", StatInfo.ColNameList[cidx], 
	      StatInfo.ColStatList[cidx]);
    else 
      sprintf(idstr, "%i", StatInfo.ColNumList[cidx] );	
    if ( StatType == 'A' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      idstr);
    else if ( StatType == 'S' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      SeasonNames[date.season], idstr);
    else if ( StatType == 'M' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      MonthNames[date.month-1], idstr);
    else if ( StatType == 'W' )
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      date.doy, date.doy+6, idstr);
    else 
      sprintf(GridFileName[sidx][cidx], GridName, 9999, 
	      date.doy, idstr);
  }

  // set default no error
  ErrNum = 0;

  // Create blank file
  if ( OVERWRITE ) // obliterate any existing file
    if ( ( fin = fopen( GridFileName[sidx][cidx], "w" ) ) == NULL )
      ErrNum = -1; // cannot open file for some reason, report error
  else if ( ( fin = fopen(GridFileName[sidx][cidx], "r") ) == NULL ) // no previous file
    if ( ( fin = fopen( GridFileName[sidx][cidx], "w" ) ) == NULL ) // create new file
      ErrNum = -1; // Failed to create new file, report error
  fclose ( fin ); 

  return(ErrNum);

}

char CheckExistingXyzFiles ( char              StatType,
			     char             *GridName,
			     StatInfoStruct    StatInfo, 
			     char           ***GridFileName,
			     DATE_STRUCT       date,
			     int               sidx,
			     int               cidx,
			     int              *nrows,
			     int              *ncols,
			     double           *ll_lat,
			     double           *ll_lng,
			     double           *cellsize,
			     int              *NODATA ) {

  FILE *fin;
  char idstr[250];
  float MaxLat, MinLat, MaxLng, MinLng, tmplat, tmplng;

  if ( StatInfo.UseColFile )
    sprintf(idstr, "%s_%s", StatInfo.ColNameList[cidx], 
	    StatInfo.ColStatList[cidx]);
  else 
    sprintf(idstr, "%i", StatInfo.ColNumList[cidx] );	
  if ( StatType == 'A' )
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    idstr);
  else if ( StatType == 'S' )
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    SeasonNames[date.season], idstr);
  else if ( StatType == 'M' )
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    MonthNames[date.month-1], idstr);
  else if ( StatType == 'W' )
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    date.doy, date.doy+6, idstr);
  else 
    sprintf(GridFileName[sidx][cidx], GridName, date.year, 
	    date.doy, idstr);

  // Create blank file
  if ( ( fin = fopen(GridFileName[sidx][cidx], "r") ) == NULL )
    return (TRUE);
  else {
    MaxLat = -LARGE_CHK_VAL;
    MinLat = +LARGE_CHK_VAL;
    MaxLng = -LARGE_CHK_VAL;
    MinLng = +LARGE_CHK_VAL;
    while ( fscanf( fin,"%*s %f %f", &tmplat, &tmplng ) != EOF ) {
      if (tmplat > MaxLat ) MaxLat = tmplat;
      if (tmplat < MinLat ) MinLat = tmplat;
      if (tmplng > MaxLng ) MaxLng = tmplng;
      if (tmplng < MinLng ) MinLng = tmplng; 
    }

    /***** THIS IS WHERE I WAS IN EDITING, WANT THIS TO ESTIMATE LAT AND LONG RANGE *****/
    fclose ( fin ); 
    return (FALSE);
  }

  return(-1);

}

