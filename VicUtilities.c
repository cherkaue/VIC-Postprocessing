#include <VicStuff.h>

#define MaxColumns 1000

int get_header( gzFile **fin, 
		char  ***ColNames, 
		char   **ColTypes, 
		double  **ColMults, 
		char   **ColAggTypes, 
		int     *TimeStep, 
		int     *NumLayers, 
		int     *NumNodes, 
		int     *NumBands, 
		int     *NumFrostFronts, 
		int     *NumLakeNodes,
		int     *NumCols,
		int     *NumBytes ) {
  /************************************************************************
    This routine checks the file header to assess whether the file is binary
    or ASCII, it sets a flag to that effect and then calls the appropriate 
    subroutine to read the header information.
  ************************************************************************/

  int vidx, UsedBytes;
  int BinaryFile;
  unsigned short int Identifier[4];  // 0xFFFF, repeated 4 times

  // Check for correct header form
  gzread(fin[0],Identifier,4*sizeof(unsigned short int));
  if ( Identifier[0] == 0xFFFF && Identifier[1] == 0xFFFF && Identifier[2] == 0xFFFF && Identifier[3] == 0xFFFF ) {
    // file is binary
    BinaryFile = TRUE;
    get_header_BINARY( fin, ColNames, ColTypes, ColMults, ColAggTypes, 
		       TimeStep, NumLayers, NumNodes, NumBands, 
		       NumFrostFronts, NumLakeNodes, NumCols, NumBytes );
  }
  else { 
    BinaryFile = FALSE;
    gzrewind( fin[0] );
    get_header_ASCII( fin, ColNames, ColTypes, ColMults, ColAggTypes, 
		      TimeStep, NumLayers, NumNodes, NumBands, 
		      NumFrostFronts, NumLakeNodes, NumCols, NumBytes );
  }

  return ( BinaryFile );

}

int get_header_BINARY( gzFile **fin, 
		       char  ***ColNames, 
		       char   **ColTypes, 
		       double  **ColMults, 
		       char   **ColAggTypes, 
		       int     *TimeStep, 
		       int     *NumLayers, 
		       int     *NumNodes, 
		       int     *NumBands, 
		       int     *NumFrostFronts, 
		       int     *NumLakeNodes,
		       int     *NumCols,
		       int     *NumBytes ) {
  /************************************************************************
    This routine reads the binary VIC model file header (VIC 4.1.0 r4 and later)
    that includes information on the model simulation and format of the output 
    file.
  ************************************************************************/

  int vidx, UsedBytes;

  unsigned short int Identifier[4];  // 0xFFFF, repeated 4 times
  unsigned short int Nbytes;         // Number of bytes in the header, including Identifier
  // Part 1: Global Attributes
  unsigned short int Nbytes1;        // Number of bytes in part 1
  int             nrecs;          // Number of records in the file
  //int             dt;             // Time step length in hours
  int             startyear;      // Year of first record
  int             startmonth;     // Month of first record
  int             startday;       // Day of first record
  int             starthour;      // Hour of first record
  char            ALMA_OUTPUT;    // 0 = standard VIC units; 1 = ALMA units
  char            Nvars;          // Number of variables in the file, including date fields
  // Part 2: Variables
  unsigned short    Nbytes2;        // Number of bytes in part 2
  // For each variable, the following fields: { len varname type mult }
  char            len;            // Number of characters in varname
  char           *varname;        // Variable name [len]
  char            type;           // Code identifying variable type
  double           mult;           // Multiplier for variable

  // Initialize model parameters
  (*NumLayers)      = 0;
  (*NumNodes)       = 0;
  (*NumBands)       = 0;
  (*NumFrostFronts) = 0;
  (*NumLakeNodes)   = 0;

  // Get header length, useful for skipping header
  gzread(fin[0],&Nbytes,sizeof(unsigned short int));

  // Process global attributes
  gzread(fin[0],&Nbytes1,sizeof(unsigned short int));
  gzread(fin[0],&nrecs,sizeof(int));
  gzread(fin[0],TimeStep,sizeof(int));
  gzread(fin[0],&startyear,sizeof(int));
  gzread(fin[0],&startmonth,sizeof(int));
  gzread(fin[0],&startday,sizeof(int));
  gzread(fin[0],&starthour,sizeof(int));
  gzread(fin[0],&ALMA_OUTPUT,sizeof(char));
  gzread(fin[0],&Nvars,sizeof(char));
  (*NumCols) = (int)Nvars;
  (*ColNames) = (char **)malloc(MaxColumns*sizeof(char *));
  (*ColTypes) = (char *) malloc(MaxColumns*sizeof(char));
  (*ColMults) = (double *)malloc(MaxColumns*sizeof(double));
  (*ColAggTypes) = (char *) malloc(MaxColumns*sizeof(char));

  // Process Variables
  gzread(fin[0],&Nbytes2,sizeof(unsigned short int));
  UsedBytes = sizeof(unsigned short int);
  (*NumBytes) = 0;
  vidx = 0;
  while( UsedBytes < Nbytes2 ) {
    if ( vidx >= MaxColumns ) {
      fprintf( stderr, "ERROR: Too many columns in the current file (%i), need to increase MaxColumns from %i.\n", vidx, MaxColumns );
      return (-1);
    }
    gzread(fin[0],&len,sizeof(char));
    varname = (char *)calloc( ((int)len+1), sizeof(char) );
    gzread(fin[0],varname,(int)len*sizeof(char));
    (*ColNames)[vidx] = varname;
    gzread(fin[0],&type,sizeof(char));
    (*ColTypes)[vidx] = type;
    gzread(fin[0],&mult,sizeof(float));
    (*ColMults)[vidx] = mult;
    UsedBytes += sizeof(char) + (int)len*sizeof(char) + sizeof(char) + sizeof(float);

    // Update number of bytes per input line
    switch ( type ) {
    case OUT_TYPE_CHAR:
      (*NumBytes) += sizeof(char);
      break;
    case OUT_TYPE_SINT:
    case OUT_TYPE_USINT:
      (*NumBytes) += sizeof(short int);
      break;
    case OUT_TYPE_INT:
      (*NumBytes) += sizeof(int);
      break;
    case OUT_TYPE_FLOAT:
      (*NumBytes) += sizeof(float);
      break;
    case OUT_TYPE_DOUBLE:
      (*NumBytes) += sizeof(double);
      break;
    default:
      fprintf( stderr, "WARNING: Using default variable type - DOUBLE.\n" );
      (*NumBytes) += sizeof(double);
      break;
    }

    // Determine daily from sub-daily estimate type (based on create_output_list.c from VIC model)
    (*ColAggTypes)[vidx] = GetVarAggType( (*ColNames)[vidx] );

    // Update model parameters
    if ( strncmp( varname, "OUT_SOIL_LIQ", 12 ) == 0 ) (*NumLayers)++;
    if ( strncmp( varname, "OUT_SOIL_TNODE", 14 ) == 0 ) (*NumNodes)++;
    if ( strncmp( varname, "OUT_SWE_BAND", 12 ) == 0 ) (*NumBands)++;
    if ( strncmp( varname, "OUT_FDEPTH", 10 ) == 0 ) (*NumFrostFronts)++;
    if ( strncmp( varname, "OUT_LAKE_NODE", 13 ) == 0 ) (*NumLakeNodes)++;

    vidx++;

  }
  (*NumCols) = vidx;
  (*ColNames) = (char **)realloc((*ColNames),vidx*sizeof(char *));
  (*ColTypes) = (char *) realloc((*ColTypes),vidx*sizeof(char));
  (*ColMults) = (double *)realloc((*ColMults),vidx*sizeof(double));  
  (*ColAggTypes) = (char *)realloc((*ColAggTypes),vidx*sizeof(char));  

  return (0);

}

int get_header_ASCII( gzFile **fin, 
		      char  ***ColNames, 
		      char   **ColTypes, 
		      double  **ColMults, 
		      char   **ColAggTypes, 
		      int     *TimeStep, 
		      int     *NumLayers, 
		      int     *NumNodes, 
		      int     *NumBands, 
		      int     *NumFrostFronts, 
		      int     *NumLakeNodes,
		      int     *NumCols,
		      int     *NumBytes ) {
  /************************************************************************
    This routine reads the ASCII VIC model file header (VIC 4.1.0 r4 and later)
    that includes information on the model simulation and format of the output 
    file.
  ************************************************************************/

  int vidx, UsedBytes;

  char           *TmpStr, *ErrStr;
  char            delimiters[] = " \t\n";
  unsigned short int Identifier[4];  // 0xFFFF, repeated 4 times
  unsigned short int Nbytes;         // Number of bytes in the header, including Identifier
  // Part 1: Global Attributes
  unsigned short int Nbytes1;        // Number of bytes in part 1
  int             nrecs;          // Number of records in the file
  //int             dt;             // Time step length in hours
  int             startyear;      // Year of first record
  int             startmonth;     // Month of first record
  int             startday;       // Day of first record
  int             starthour;      // Hour of first record
  char            ALMA_OUTPUT;    // 0 = standard VIC units; 1 = ALMA units
  char            Nvars;          // Number of variables in the file, including date fields
  // Part 2: Variables
  unsigned short    Nbytes2;        // Number of bytes in part 2
  // For each variable, the following fields: { len varname type mult }
  char            len;            // Number of characters in varname
  char           *varname;        // Variable name [len]
  char            type;           // Code identifying variable type
  double           mult;           // Multiplier for variable

  // Initialize model parameters
  (*NumLayers)      = 0;
  (*NumNodes)       = 0;
  (*NumBands)       = 0;
  (*NumFrostFronts) = 0;
  (*NumLakeNodes)   = 0;

  TmpStr = (char *)malloc(MaxCharData*sizeof(char));

  // process global attributes
  ErrStr = gzgets( fin[0],TmpStr, MaxCharData );
  sscanf( TmpStr, "# %*s %i\n", &nrecs ); 
  ErrStr = gzgets( fin[0],TmpStr, MaxCharData );
  sscanf( TmpStr, "# %*s %i\n", TimeStep ); 
  ErrStr = gzgets( fin[0],TmpStr, MaxCharData );
  sscanf( TmpStr, "# %*s %i-%i-%i %i:%*i:%*i\n", &startyear, &startmonth, &startday, &starthour ); 
  ErrStr = gzgets( fin[0], TmpStr, MaxCharData );
  sscanf( TmpStr, "# %*s %i\n", &ALMA_OUTPUT ); 
  ErrStr = gzgets( fin[0],TmpStr, MaxCharData );
  sscanf( TmpStr, "# %*s %i\n", &Nvars ); 

  // initialize variable attribute arrays
  (*NumCols) = (int)Nvars;
  (*ColNames) = (char **)malloc(MaxColumns*sizeof(char *));
  (*ColTypes) = (char *) calloc(MaxColumns,sizeof(char));
  (*ColMults) = (double *)calloc(MaxColumns,sizeof(double));
  (*ColAggTypes) = (char *) calloc(MaxColumns,sizeof(char));

  // Process Variables
  gzgets( fin[0], TmpStr, MaxCharData );
  varname = strtok( TmpStr, delimiters ); // skip # header marker
  vidx = 0;
  varname = strtok( '\0', delimiters ); // get next variable name
  do {
    // store variable name
    (*ColNames)[vidx] = (char *)malloc(50*sizeof(char));
    strcpy( (*ColNames)[vidx], varname );
    // Determine daily from sub-daily estimate type (based on create_output_list.c from VIC model)
    (*ColAggTypes)[vidx] = GetVarAggType( (*ColNames)[vidx] );
    // Update model parameters
    if ( strncmp( varname, "OUT_SOIL_LIQ", 12 ) == 0 ) (*NumLayers)++;
    if ( strncmp( varname, "OUT_SOIL_TNODE", 14 ) == 0 ) (*NumNodes)++;
    if ( strncmp( varname, "OUT_SWE_BAND", 12 ) == 0 ) (*NumBands)++;
    if ( strncmp( varname, "OUT_FDEPTH", 10 ) == 0 ) (*NumFrostFronts)++;
    if ( strncmp( varname, "OUT_LAKE_NODE", 13 ) == 0 ) (*NumLakeNodes)++;

    vidx++;

    varname = strtok( '\0', delimiters ); // get next variable name

  } while ( varname );
  
  (*NumCols) = vidx;
  (*ColNames) = (char **)realloc((*ColNames),vidx*sizeof(char *));
  (*ColTypes) = (char *) realloc((*ColTypes),vidx*sizeof(char));
  (*ColMults) = (double *)realloc((*ColMults),vidx*sizeof(double));  
  (*ColAggTypes) = (char *)realloc((*ColAggTypes),vidx*sizeof(char));  

  free ((char *)TmpStr);

  (*NumBytes) = (*NumCols); // Bytes not used for ASCII, so set equal to number of columns
  
  return (0);
  
}
 
char GetVarAggType( char *VarName ) {

  char AggType;

  // Determine daily from sub-daily estimate type (based on create_output_list.c from VIC model)
  if ( strcmp( VarName, "OUT_LAKE_DEPTH" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_LAKE_ICE" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_LAKE_ICE_FRACT" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_LAKE_ICE_HEIGHT" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_LAKE_MOIST" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_LAKE_SURF_AREA" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_LAKE_VOLUME" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_ROOTMOIST" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SMFROZFRAC" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SMLIQFRAC" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SNOW_CANOPY" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SNOW_COVER" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SNOW_DEPTH" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SOIL_ICE" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SOIL_LIQ" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SOIL_MOIST" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SOIL_WET" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SURFSTOR" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SURF_FROST_FRAC" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SWE" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_WDEW" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_LAKE_ICE_TEMP" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_LAKE_SURF_TEMP" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SNOW_CANOPY_BAND" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SNOW_COVER_BAND" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SNOW_DEPTH_BAND" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SWE_BAND" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SOIL_DEPTH" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_POROSITY" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_ZSUM_NODE" ) == 0 ) 
    AggType = (char)AGG_TYPE_END;
  else if ( strcmp( VarName, "OUT_SUBSIDENCE" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_BASEFLOW" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_DELINTERCEPT" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_DELSOILMOIST" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_DELSWE" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_DELSURFSTOR" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_EVAP" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_EVAP_BARE" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_EVAP_CANOP" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_EVAP_LAKE" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_INFLOW" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_PREC" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_RAINF" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_REFREEZE" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_RUNOFF" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_SNOW_MELT" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_SNOWF" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_SUB_BLOWING" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_SUB_CANOP" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_SUB_SNOW" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_SUB_SURFACE" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_TRANSP_VEG" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_DELTACC_BAND" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else if ( strcmp( VarName, "OUT_SNOW_MELT" ) == 0 ) 
    AggType = (char)AGG_TYPE_SUM;
  else {
    //fprintf( stderr, "WARNING: Unknown aggregation type for variable %s, using AGG_TYPE_AVG by default.\n", VarName );
    AggType = (char)AGG_TYPE_AVG;
  }

  return( AggType );

}

int get_record_PEN( int            BinaryFile,
		    char         **RawData, 
		    int           *RawPtr,
		    char         **ColNames,
		    char          *ColTypes,
		    double        *ColMults,
		    int           *date, 
		    double        *data, 
		    int            NumCols,
                    int            NumOut,
		    int           *DataStart,
		    int            CalcPE,
		    PenInfoStruct  PenInfo,
		    int            CalcTR,
		    PenInfoStruct  TRoffInfo,
		    int            CalcTSM,
		    PenInfoStruct  TSMInfo ) {
  /********************************************************************************
    Makes sure that the correct get_record function is called based on the type
    of data file being processed.

    Modifications:
    16-Jan-2017 Modified to pass flags indicating whether or not PE and 
      TOTAL_RUNOFF should be calculated fromt he current data set.   KAC
    30-May-2017 Added flags for calculating modified growing degree days 
      and chilling hours.                                 KAC

  ********************************************************************************/

  int                 cidx, vidx;
  double Rnet, G, U, Ts, RH, Ta;
  double SVP, vpd, delta, lv, gamma, A, denominator;

  if ( BinaryFile ) 
    get_record_BINARY( RawData, RawPtr, ColNames, ColTypes, ColMults,
		       date, data, NumCols, DataStart );
  else
    get_record_ASCII( RawData, RawPtr, ColNames, ColTypes, ColMults,
		       date, data, NumCols, DataStart );

  /***** Calculate PE and add to data[cidx] *****/
  if ( CalcPE ) {
    Rnet = data[PenInfo.ColNumList[0]]; // net radiation
    G = data[PenInfo.ColNumList[1]]; // ground heat flux
    U = data[PenInfo.ColNumList[2]]; // wind speed
    Ts = data[PenInfo.ColNumList[3]]; // surface temperature
    RH = data[PenInfo.ColNumList[4]]; // relative humidity
    Ta = data[PenInfo.ColNumList[5]]; // air temperature
    
    SVP = 0.6108 * exp((17.27 * Ts)/(237.3+Ts));
    if(Ts<0) SVP *= 1.0 + .00972 * Ts + .000042 * Ts;
    
    vpd = SVP*(1. - RH/100.);
    
    delta = 4098.*SVP / ((237.3 + Ts) * (237.3 + Ts));

    /* calculate latent heat of vaporization (J/kg). Eq. 4.2.1 in Handbook of 
       Hydrology, assume Ts is Tair */
    lv = 2501000 - 2361 * Ta;
  
    /* calculate gamma. Eq. 4.2.28. Handbook of Hydrology (kPa/C) */
    gamma = 1628.6 * 101.3/(lv);
    A = (Rnet-G)*86400./lv;
    denominator = delta + gamma*(1.+.33*U);

    data[PenInfo.ColNumList[6]] = ((delta/denominator)*A + (gamma/denominator)*(900./(Ta+275.))*U*vpd);
  }

  /***** Compute total runoff and add to data[cidx] *****/
  if ( CalcTR )
    data[TRoffInfo.ColNumList[2]] = data[TRoffInfo.ColNumList[0]] + data[TRoffInfo.ColNumList[1]];

  /***** Compute total soil moisture and add to data[cidx] *****/
  if ( CalcTSM )
    data[TSMInfo.ColNumList[3]] = data[TSMInfo.ColNumList[0]] + data[TSMInfo.ColNumList[1]] + data[TSMInfo.ColNumList[2]];

  return (0);

}

int get_record_BINARY( char  **RawData, 
		       int    *RawPtr,
		       char  **ColNames,
		       char   *ColTypes,
		       double  *ColMults,
		       int    *date, 
		       double  *data, 
		       int     NumCols,
		       int    *DataStart ) {
  /***********************************************************************************
    Parses a binary VIC model output file based on information from the header file.
  ************************************************************************************/
  
  extern char *VarType[];

  char               *cptr, IntVal;
  short int          *siptr;
  unsigned short int *usiptr;
  int                *iptr, TmpInt;
  float              *fptr;
  double             *dptr, TmpDouble;
  int                 cidx, vidx;

  cidx = 0;
  for ( vidx = 0; vidx < NumCols; vidx++ ) {
    IntVal = 1;
    // read column variable
    switch ( ColTypes[vidx] ) {
    case OUT_TYPE_CHAR:
      cptr = (char *)&RawData[0][(*RawPtr)];
      TmpInt = (int)cptr[0];
      (*RawPtr) += sizeof(char);
      break;
    case OUT_TYPE_SINT:
      siptr = (short int *)&RawData[0][(*RawPtr)];
      TmpInt = (int)siptr[0];
      (*RawPtr) += sizeof(short int);
      break;
    case OUT_TYPE_USINT:
      usiptr = (unsigned short int *)&RawData[0][(*RawPtr)];
      TmpInt = (int)usiptr[0];
      (*RawPtr) += sizeof(unsigned short int);
      break;
    case  OUT_TYPE_INT:
      iptr = (int *)&RawData[0][(*RawPtr)];
      TmpInt = (int)iptr[0];
      (*RawPtr) += sizeof(int);
      break;
    case OUT_TYPE_FLOAT:
      fptr = (float *)&RawData[0][(*RawPtr)];
      TmpDouble = (float)fptr[0];
      IntVal = 0;
      (*RawPtr) += sizeof(float);
      break;
    case OUT_TYPE_DOUBLE:
      dptr = (double *)&RawData[0][(*RawPtr)];
      TmpDouble = (double)dptr[0];
      IntVal = 0;
      (*RawPtr) += sizeof(double);
      break;
    default:
      fprintf( stderr, "WARNING: Unknown variable type (%i) for variable %i.\n", (int)ColTypes[vidx], vidx );
      return (-1);
    }

    // Handle date information
    if ( strcmp( ColNames[vidx], "YEAR" ) == 0 ) {
      // year
      if ( IntVal ) date[0] = TmpInt;
      else date[0] = (int)TmpDouble;
    }
    else if ( strcmp( ColNames[vidx], "MONTH" ) == 0 ) {
      // month 
      if ( IntVal ) date[1] = TmpInt;
      else date[1] = (int)TmpDouble;
    }
    else if ( strcmp( ColNames[vidx], "DAY" ) == 0 ) {
      // day 
      if ( IntVal ) date[2] = TmpInt;
      else date[2] = (int)TmpDouble;
    }
    else if ( strcmp( ColNames[vidx], "HOUR" ) == 0 ) {
      // hour
      if ( IntVal ) date[3] = TmpInt;
      else date[3] = (int)TmpDouble;
    }

    // Process data variables
    else {
      // Store current variable
      if ( cidx == 0 ) (*DataStart) = vidx;
      if ( IntVal ) data[cidx] = (double)TmpInt/ColMults[vidx];
      else data[cidx] = (double)TmpDouble/ColMults[vidx];
      cidx++;
    }
  }

  return (0);

}

int get_record_ASCII( char  **RawData, 
		      int    *RawPtr,
		      char  **ColNames,
		      char   *ColTypes,
		      double  *ColMults,
		      int    *date, 
		      double  *data, 
		      int     NumCols,
		      int    *DataStart ) {
  /***********************************************************************************
    Parses a ASCII VIC model output file based on information from the header file.
  ************************************************************************************/
  
  extern char *VarType[];

  char           *TmpStr, *ErrStr;
  char            delimiters[] = " \t\n";
  char           *varValue;
  int             cidx, vidx;

  vidx = 0;
  cidx = 0;
  TmpStr = (char *)malloc(MaxCharData*sizeof(char));
  // read next line from file
  strcpy( TmpStr, RawData[(*RawPtr)] );
  // process all variables in line
  varValue = strtok( TmpStr, delimiters ); // skip # header marker
  while ( varValue && vidx < NumCols ) { 

    // Handle date information
    if ( strcmp( ColNames[vidx], "YEAR" ) == 0 ) {
      // year
      date[0] = atoi(varValue);
    }
    else if ( strcmp( ColNames[vidx], "MONTH" ) == 0 ) {
      // month 
      date[1] = atoi(varValue);
    }
    else if ( strcmp( ColNames[vidx], "DAY" ) == 0 ) {
      // day 
      date[2] = atoi(varValue);
    }
    else if ( strcmp( ColNames[vidx], "HOUR" ) == 0 ) {
      // hour
      date[3] = atoi(varValue);
    }

    // Process data variables
    else {
      // Store current variable
      data[cidx] = atof(varValue);
      cidx++;
    }

    // Get next variable
    varValue = strtok( '\0', delimiters ); // get next variable name
    vidx++;

  }

  if ( vidx != NumCols ) {
    fprintf( stderr, "WARNING: Incomplete line in data file.\n" );
    return(-1);
  }

  (*RawPtr) ++;
  free ((char*)TmpStr);

  return (0);

}

char *get_next_string(char *str, int *index, char sepchar) {
  
  /*************************************************************
    This routine extracts items from a string, *str, separated
    by sepchar.  *index points to the current location in the
    string *str.  To extract the first N items for *str, set
    *index to 0 in the calling routine, and call this function
    N times sending the updated value of *index each time.
  *************************************************************/
  int idx, pos=0;
  char *tmpstr;
 
  tmpstr = (char *)calloc(256,sizeof(char));
 
  idx = *index;
 
  while(str[idx]!=sepchar && str[idx]!='\0' && str[idx]!='\n') {
    tmpstr[pos] = str[idx];
    pos++;
    idx++;
  }
  tmpstr[pos] = '\0';
  idx++;
 
  *index = idx;
 
  return(tmpstr);
 
}

char *reset_spaces( char *TmpStr ) {
  // routine converts "_" characters in field and dataspace names back into " "s
  int pos;

  if ( !TmpStr ) return( TmpStr );
  pos = 0;
  while( TmpStr[pos] != '\0' ) {
    if ( TmpStr[pos] == '_' ) 
      TmpStr[pos] = ' ';
    pos++;
  }
  return TmpStr;
}

