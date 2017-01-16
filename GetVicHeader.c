#include <VicStuff.h>

int main(int argc, char *argv[])
/**********************************************************************
  GetVicHeader.c          Keith Cherkauer          February 7, 2008

  This program extracts information from a VIC binary output header 
  of a given file (created with PRT_HEADER = TRUE in the VIC global 
  file).  Information is output to stdout.

**********************************************************************/
{

  gzFile *InFile;
  char    filename[2048], AggText[50];;
  char  **ColNames, *ColTypes, *ColAggTypes;
  float  *ColMults;
  int     TimeStep, NumLayers, NumNodes, NumBands, NumFrostFronts;
  int     NumLakeNodes, NumCols, ErrNum, vidx, NumBytes;

  // check if usage message should be printed
  if ( argc != 2 ) {
    fprintf(stderr,"\nUsage: %s <filename>\n", argv[0]);
    fprintf(stderr,"\n\tThis program extracts information from a VIC model binary header at the start \n\tof a given file (created with PRT_HEADER = TRUE in \n\tthe VIC global file).  Information is output to stdout.\n");
    fprintf(stderr,"\n\t<filename> is the name of anbinary VIC model output file \n\t\twith header information. \n\n");
    exit(0);
  }

  // handle command line arguments
  strcpy(filename,argv[1]);

  // Open file 
  if ( ( InFile  = gzopen(filename, "rb") ) == NULL ) {
    printf("File %s does not exist\n", filename);
    exit(1);
  }

  // Read header information
  ErrNum = get_header_NEW( &InFile, &ColNames, &ColTypes, &ColMults, 
			   &ColAggTypes, &TimeStep, &NumLayers, 
			   &NumNodes, &NumBands, &NumFrostFronts, 
			   &NumLakeNodes, &NumCols, &NumBytes );

  // Write header information to stdout
  fprintf( stdout, "\n" );
  for( vidx = 0; vidx < NumCols; vidx++ ) {
    switch (ColAggTypes[vidx]) {
    case AGG_TYPE_AVG:
      strcpy( AggText, "Average Value" );
      break;
    case AGG_TYPE_BEG:
      strcpy( AggText, "Begin Value" );
      break;
    case AGG_TYPE_END:
      strcpy( AggText, "End Value" );
      break;
    case AGG_TYPE_MAX:
      strcpy( AggText, "Maximum Value" );
      break;
    case AGG_TYPE_MIN:
      strcpy( AggText, "Minimum Value" );
      break;
    case AGG_TYPE_SUM:
      strcpy( AggText, "Sum Value" );
      break;
    default:
      strcpy( AggText, "Unknown" );
    }
    
    fprintf( stdout, "Column %2i:%30s%5i%12g (%s)\n", vidx, ColNames[vidx], ColTypes[vidx], ColMults[vidx], AggText );
  }
  fprintf( stdout, "\nTime Step: %i\n", TimeStep );
  fprintf( stdout, "Number of Layers: %i\n", NumLayers );
  fprintf( stdout, "Number of Soil Thermal Nodes: %i\n", NumNodes );
  fprintf( stdout, "Number of Elevation Bands: %i\n", NumBands );
  fprintf( stdout, "Number of Frost Fronts: %i\n", NumFrostFronts );
  fprintf( stdout, "Number of Lake Nodes: %i\n", NumLakeNodes );
  fprintf( stdout, "Number of Columns: %i\n\n", NumCols );

  return (0);

}
  
