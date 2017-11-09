#include <VicStuff.h>

int main(int argc, char *argv[])
/**********************************************************************
  VicBinaryDump2ASCII.c     Keith Cherkauer        February 11, 2008

  This program reads VIC model output files with headers (post 4.1.0 r4)
  and extracts actual values for all model output variables or for a 
  subset of those variables if provided.  Note that columns will be in 
  the same order as they appear in the VIC model output file no matter
  what order they were listed on the command line.  But only those 
  variables listed on the command line will be in the output file.

  Modified:
  11-Fab-2008 Modified from VicLdasDump2ASCII.c             KAC  
  11-Feb-2008 Removed calls to reset_spcaes to leave "_"s in 
              variable names.   KAC
  16-Jan-2017 Updated to make use of lastest version of get_header
              that will work with both binary and ASCII VIC model
              output formats.    KAC
  18-Jan-2017 Updated code to work correctly and without errors
              when given an ASCII data file.             KAC
  28-Sep-2017 Updated get_record_PEN to refelct addition of one
              more set of possible calculations due to recent 
              updates to VicOutputASMStats.c.  Also added SKIP
              flag, which is set to TRUE with gz function tries
              but fails to read from the input file.  The SKIP 
              flag then allows follow-up processing to be skipped 
              for the incomplete line.  Without it, was getting an 
              error when the last line in an ASCII file included a 
              \n, so loop passed the check for EOF, but there was
              no additional information to be read.           KAC

  ***** Seems to only work if variables in list appear in the same order 
        as in the original file. Why????? - note from VicLdasDump2ASCII.c
        *****

**********************************************************************/
{
  gzFile *InFile;
  char  **RawData;
  int    *date, UseVarList;
  double *data;
  char    name[600];
  char  **ColNames, *ColTypes, *ColAggTypes;
  double *ColMults;
  char    VarList[MaxColumns][256];
  char   *VarName, *ptr;
  int    *OutputCol;
  int     DONE, DataStart, SKIP;
  int     TimeStep, NumLayers, NumNodes, NumBands, NumFrostFronts;
  int     NumLakeNodes, NumCols, NumOut;
  int     ErrNum;
  int     idx, nidx, maxdate;
  int     NumBytes, NumRead, day, ReadBytes, RawPtr;
  int     BinaryFile;
  PenInfoStruct TempInfo; // used as a placeholder for the get_record command

  date      = (int *)malloc(NUM_DATE_VALS*sizeof(int));

  // check if usage message should be printed
  if ( argc < 2 || argc > 3 ) {
    fprintf(stderr,"\nUsage: %s <filename> [<var1>[,<var2>[,<varN>]]]\n",argv[0]);
    fprintf(stderr,"\n\tThis program reads an LDAS binary file and its header and dumps the \n\tfile contents to stdout.  If a variable list is included then only \n\tthose variables, plus date information is exported, otherwise all \n\tvariables are dumped.\n\n\tThe output file is in ASCII column format: \n");
    fprintf(stderr,"\t\t<year> <month> <day> <prec> <evap> ...\n\n");
    exit(0);
  }

  // handle command line arguments
  strcpy(name,argv[1]);
  NumOut = 0;
  UseVarList = FALSE;
  if ( argc == 3 ) {
    UseVarList = TRUE;
    // Parse variable list
    VarName = strtok_r( argv[2], ",", &ptr );
    strcpy( VarList[NumOut], VarName );
    NumOut++;
    do { 
      VarName = strtok_r( '\0', ",", &ptr );
      if ( !VarName ) break;
      strcpy( VarList[NumOut], VarName );
      NumOut++;
    } while ( VarName );
  }

  if ( ( InFile  = gzopen(name, "rb") ) == NULL ) {
    // Try opening a compressed file of the same name
    strcat(name,".gz");
    puts(name);
    if ( ( InFile  = gzopen(name, "rb") ) == NULL ) {
      printf("File %s does not exist\n", name);
      exit(1);
    }
  }

  /** Process grid flux file **/
  DONE = FALSE;
  BinaryFile = get_header( &InFile, &ColNames, &ColTypes, &ColMults, 
			   &ColAggTypes, &TimeStep, &NumLayers, &NumNodes, 
			   &NumBands, &NumFrostFronts, &NumLakeNodes, 
			   &NumCols, &NumBytes );
  if ( ErrNum < 0 ) {
    fprintf( stderr, "\t... exiting program.\n" );
    exit(-1);
  }
  data = (double *)malloc(NumCols*sizeof(double));
  fprintf( stderr, "Number of bytes per line: %i\n", NumBytes );

  // allocate arrays for reading each line of data
  NumRead = NumBytes*(int)(24/TimeStep);
  if ( BinaryFile ) {
    // Binary file format, so NumRead is actaully the number of bytes
    RawData = (char **)malloc(sizeof(char*));
    RawData[0] = (char *)malloc(NumRead*sizeof(char));
  }
  else {
    // ASCII file format, so NumRead is te number of strings to read
    RawData = (char **)malloc(sizeof(char *));
    RawData[0] = (char *)malloc(MaxCharData*sizeof(char));
  }
  OutputCol = (int *)calloc( NumCols, sizeof( int ) );
  day = 0;
  if ( TimeStep < 24 ) maxdate = 4;
  else maxdate = 3;

  // identify which columns are to be output, if specified
  if ( UseVarList ) {
    for ( idx = 0; idx < NumCols; idx++ )
      for ( nidx = 0; nidx < NumOut; nidx++ )
	if ( strcasecmp( ColNames[idx], VarList[nidx] ) == 0 )
	  OutputCol[idx] = TRUE;
  }

  // Read all data records
  while ( !gzeof( InFile ) ) {
    SKIP=FALSE;
    if ( BinaryFile ) {
      //ReadBytes = gzread(InFile,RawData,NumBytes*sizeof(char));
      if ( ( ReadBytes = gzread(InFile,RawData[0],NumRead) ) != NumRead ) {
	fprintf( stderr, "WARNING: Unable to read in complete flux file %s, got %d of %ld bytes.\n", name, ReadBytes, NumRead*sizeof(char) );
	SKIP=TRUE;
      }
    }
    else {
      if ( gzgets(InFile,RawData[0],MaxCharData) == NULL ) {
	SKIP=TRUE;
      }
    }
    if ( !SKIP ) {
      RawPtr = 0;
      ErrNum = get_record_PEN( BinaryFile, RawData, &RawPtr, ColNames, ColTypes, 
			       ColMults, date, data, NumCols, NumOut, &DataStart, 
			       FALSE, TempInfo, FALSE, TempInfo, FALSE, TempInfo );
      for ( idx = 0; idx < maxdate-1; idx++ ) fprintf( stdout, "%i\t", date[idx] );
      fprintf( stdout, "%i", date[idx] );
      for ( idx = 0; idx < NumCols-maxdate; idx++ )
	if ( !UseVarList || OutputCol[idx+maxdate] ) 
	  fprintf( stdout, "\t%f", data[idx] );
      fprintf( stdout, "\n" );
    }
  }

  // close input file
  gzclose(InFile);

  // free memory
  free((char *)date);
  free((char *)data);
  free((char *)OutputCol);

  return (0);

}
  
