#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include <DateFuncs.h>
#include <VicStats_1d.h>

#ifndef NUM_DATE_VALS
#define NUM_DATE_VALS 4
#endif

#ifndef NoData
#define NoData -999
#endif

#ifndef NumBlocks  
#define NumBlocks 1       // Number of four year blocks to process concurrently
#define NumBlockDays 1461 // Number of days in a four year block
#endif

#ifndef LARGE
#define LARGE_CHK_VAL 1e+20
#define SMALL_CHK_VAL 1e-4
#endif

#ifndef TRUE
#define FALSE 0
#define TRUE !FALSE
#endif

#ifndef FAIL
#define FAIL -1
#endif

#ifndef MaxColumns
#define MaxColumns 1000
#endif

#define MaxCharData 5000 // Maximum number of characters per line of ASCII input file.

/***** Output BINARY format types (from vicNl_def.h) *****/
#define OUT_TYPE_DEFAULT 0 /* Default data type */
#define OUT_TYPE_CHAR    1 /* chaßßr */
#define OUT_TYPE_SINT    2 /* short int */
#define OUT_TYPE_USINT   3 /* unsigned short int */
#define OUT_TYPE_INT     4 /* int */
#define OUT_TYPE_FLOAT   5 /* single-precision floating point */
#define OUT_TYPE_DOUBLE  6 /* double-precision floating point */

/***** Output aggregation method types (from vicNl_def.h) *****/
#define AGG_TYPE_AVG     0 /* average over agg interval */
#define AGG_TYPE_BEG     1 /* value at beginning of agg interval */
#define AGG_TYPE_END     2 /* value at end of agg interval */
#define AGG_TYPE_MAX     3 /* maximum value over agg interval */
#define AGG_TYPE_MIN     4 /* minimum value over agg interval */
#define AGG_TYPE_SUM     5 /* sum over agg interval */

/***** Define data type structures *****/
typedef struct {
  char **ColNameList;
  char **ColStatList;
  char UseColFile;
  int *ColNumList;
  int MaxColNum;
  int Ncols;
  float *Thres;
} StatInfoStruct; // structure for storing information on summary statistics

typedef struct {
  char **ColNameList;
  int *ColNumList;
  int MaxColNum;
  int Ncols;
  int *OutputCols;
  int Noutput;
} PenInfoStruct; // structure for storing column information for additional calculations (PE and Total Runoff)

// subroutine prototypes

int get_header( gzFile **, char ***, char **, double **, char **, int *, int *, 
		int *, int *, int *, int *, int *, int * );
int get_header_ASCII( gzFile **, char ***, char **, double **, char **, int *, int *, 
		int *, int *, int *, int *, int *, int * );
int get_header_BINARY( gzFile **, char ***, char **, double **, char **, int *, int *, 
		int *, int *, int *, int *, int *, int * );
int    get_record_PEN( int, char **, int *, char **, char *, double *, int *, double *, 
		       int, int, int *, int, PenInfoStruct, int, PenInfoStruct );
int    get_record_ASCII( char **, int *, char **, char *, double *, int *, double *, int, 
		       int * );
int    get_record_BINARY( char **, int *, char **, char *, double *, int *, double *, int, 
			  int * );
char GetVarAggType( char * );
/***
int    get_header_NEW( gzFile  **fin, char ***, char **, float **, char **, 
		       int *, int *, int *, int *, int *, int *, int *, int * );
int    get_record_NEW( char *, int *, char **, char *, float *, int *, float *, int, 
		       int * );
***/
char *get_next_string(char *, int *, char );
char *reset_spaces( char * );
