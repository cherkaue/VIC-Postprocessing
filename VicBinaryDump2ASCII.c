#include <VicStuff.h>

int main(int argc, char *argv[])
/**********************************************************************
  VicBinaryDump2ASCII.c     Keith Cherkauer        February 11, 2008

  This program reads VIC model binary files with headers (post 4.1.0 r4)
  and extracts actual values for all model output variables for the 
  specified period of interest.  

  Modified:
  11-Fab-2008 Modified from VicLdasDump2ASCII.c             KAC  
  11-Feb-2008 Removed calls to reset_spcaes to leave "_"s in 
              variable names.   KAC

  ***** Seems to only work if variables in list appear in the same order 
        as in the original file. Why????? - note from VicLdasDump2ASCII.c
        *****

**********************************************************************/
{
  gzFile *InFile;
  char   *RawData;
  int    *date, UseVarList;
  float  *data;
  char    name[600];
  char  **ColNames, *ColTypes, *ColAggTypes;
  float  *ColMults;
  char    VarList[MaxColumns][256];
  char   *VarName, *ptr;
  int    *OutputCol;
  int     DONE, DataStart;
  int     TimeStep, NumLayers, NumNodes, NumBands, NumFrostFronts;
  int     NumLakeNodes, NumCols, NumOut;
  int     ErrNum;
  int     idx, nidx, maxdate;
  int     NumBytes, NumRead, day, ReadBytes, RawPtr;

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
  ErrNum = get_header_NEW( &InFile, &ColNames, &ColTypes, &ColMults, 
			     &ColAggTypes, &TimeStep, &NumLayers, &NumNodes, 
			     &NumBands, &NumFrostFronts, &NumLakeNodes, 
			     &NumCols, &NumBytes );
  if ( ErrNum < 0 ) {
    fprintf( stderr, "\t... exiting program.\n" );
    exit(-1);
  }
  data = (float *)malloc(NumCols*sizeof(float));
  fprintf( stderr, "Number of bytes per line: %i\n", NumBytes );

  // allocate arrays for reading each line of data
  NumRead = NumBytes*(int)(24/TimeStep);
  RawData = (char *)malloc(NumRead*sizeof(char));
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
    ReadBytes = gzread(InFile,RawData,NumBytes*sizeof(char));
    RawPtr = 0;
    ErrNum = get_record_NEW( RawData, &RawPtr, ColNames, ColTypes, 
			     ColMults, date, data, NumCols, &DataStart );
    for ( idx = 0; idx < maxdate-1; idx++ ) fprintf( stdout, "%i\t", date[idx] );
    fprintf( stdout, "%i", date[idx] );
    for ( idx = 0; idx < NumCols-maxdate; idx++ )
      if ( !UseVarList || OutputCol[idx+maxdate] ) 
	fprintf( stdout, "\t%f", data[idx] );
    fprintf( stdout, "\n" );
  }

  // close input file
  gzclose(InFile);

  // free memory
  free((char *)date);
  free((char *)data);
  free((char *)OutputCol);

  return (0);

}
  
