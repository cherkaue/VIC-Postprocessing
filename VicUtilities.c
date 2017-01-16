#include <VicStuff.h>

int get_header_NEW( gzFile **fin, 
		char  ***ColNames, 
		char   **ColTypes, 
		float  **ColMults, 
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
  float           mult;           // Multiplier for variable

  // Initialize model parameters
  (*NumLayers)      = 0;
  (*NumNodes)       = 0;
  (*NumBands)       = 0;
  (*NumFrostFronts) = 0;
  (*NumLakeNodes)   = 0;

  // Check for correct header form
  gzread(fin[0],Identifier,4*sizeof(unsigned short int));
  if ( Identifier[0] != 0xFFFF && Identifier[1] != 0xFFFF && Identifier[2] != 0xFFFF && Identifier[3] != 0xFFFF ) {
    fprintf( stderr, "ERROR: Header is not in standard binary form.  This program does not yet work with ASCII headers.\n" );
    return (-1 );
  }

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
  (*ColMults) = (float *)malloc(MaxColumns*sizeof(float));
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
      fprintf( stderr, "WARNING: Using default variable type - FLOAT.\n" );
      (*NumBytes) += sizeof(float);
      break;
    }

    // Determine daily from sub-daily estimate type (based on create_output_list.c from VIC model)
    if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_DEPTH" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_ICE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_ICE_FRACT" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_ICE_HEIGHT" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_MOIST" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_SURF_AREA" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_VOLUME" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_ROOTMOIST" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SMFROZFRAC" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SMLIQFRAC" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOW_CANOPY" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOW_COVER" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOW_DEPTH" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SOIL_ICE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SOIL_LIQ" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SOIL_MOIST" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SOIL_WET" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SURFSTOR" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SURF_FROST_FRAC" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SWE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_WDEW" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_ICE_TEMP" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_LAKE_SURF_TEMP" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOW_CANOPY_BAND" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOW_COVER_BAND" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOW_DEPTH_BAND" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SWE_BAND" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SOIL_DEPTH" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_POROSITY" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_ZSUM_NODE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_END;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SUBSIDENCE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_BASEFLOW" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_DELINTERCEPT" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_DELSOILMOIST" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_DELSWE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_DELSURFSTOR" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_EVAP" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_EVAP_BARE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_EVAP_CANOP" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_EVAP_LAKE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_INFLOW" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_PREC" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_RAINF" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_REFREEZE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_RUNOFF" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOW_MELT" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOWF" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SUB_BLOWING" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SUB_CANOP" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SUB_SNOW" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SUB_SURFACE" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_TRANSP_VEG" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_DELTACC_BAND" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else if ( strcmp( (*ColNames)[vidx], "OUT_SNOW_MELT" ) == 0 ) 
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_SUM;
    else {
      //fprintf( stderr, "WARNING: Unknown aggregation type for variable %s, using AGG_TYPE_AVG by default.\n", (*ColNames)[vidx] );
      (*ColAggTypes)[vidx] = (char)AGG_TYPE_AVG;
    }

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
  (*ColMults) = (float *)realloc((*ColMults),vidx*sizeof(float));  
  (*ColAggTypes) = (char *)realloc((*ColAggTypes),vidx*sizeof(char));  

  //(*NumLayers) = (int)TmpHdr[9];
  //(*NumNodes) = (int)TmpHdr[10];
  //(*NumBands) = (int)TmpHdr[11];
  //(*NumFrostFronts) = (int)TmpHdr[12];
  //(*NumLakeNodes) = (int)TmpHdr[13];

  return (0);

}

int get_record_NEW( char   *RawData, 
		    int    *RawPtr,
		    char  **ColNames,
		    char   *ColTypes,
		    float  *ColMults,
		    int    *date, 
		    float  *data, 
		    int     NumCols,
		    int    *DataStart ) {

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
      cptr = (char *)&RawData[(*RawPtr)];
      TmpInt = (int)cptr[0];
      (*RawPtr) += sizeof(char);
      break;
    case OUT_TYPE_SINT:
      siptr = (short int *)&RawData[(*RawPtr)];
      TmpInt = (int)siptr[0];
      (*RawPtr) += sizeof(short int);
      break;
    case OUT_TYPE_USINT:
      usiptr = (unsigned short int *)&RawData[(*RawPtr)];
      TmpInt = (int)usiptr[0];
      (*RawPtr) += sizeof(unsigned short int);
      break;
    case  OUT_TYPE_INT:
      iptr = (int *)&RawData[(*RawPtr)];
      TmpInt = (int)iptr[0];
      (*RawPtr) += sizeof(int);
      break;
    case OUT_TYPE_FLOAT:
      fptr = (float *)&RawData[(*RawPtr)];
      TmpDouble = (double)fptr[0];
      IntVal = 0;
      (*RawPtr) += sizeof(float);
      break;
    case OUT_TYPE_DOUBLE:
      dptr = (double *)&RawData[(*RawPtr)];
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
      if ( IntVal ) data[cidx] = (float)TmpInt/ColMults[vidx];
      else data[cidx] = (float)TmpDouble/ColMults[vidx];
      cidx++;
    }
  }

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

