#!/bin/env bash
# Created on March 22, 2017
#  by Keith Cherkauer
#
# This script adds a VIC model output style header to ASCII forcing data.
# By adding the header, you can use post-processing scripts that require
# the header to also process input forcing data files.  
#
# NOTE: The first three columns of the file must be YEAR, MONTH and DAY
# for daily data (add HOUR as the fourth column, if file is subdaily).
# If the forcing file does not have dates, see AddDateColumn.py to add it.
#
# Debugging location:
# /home/cherkaue/ScratchDrive/INCCIA
# Debugging command:
# /depot/phig/apps/source/VIC-Postprocessing/AddHeaderToVicForcingData.sh Indiana/ data TestOut/  36160 24 1915-01-01 4 IN_PREC IN_TMAX IN_TMIN IN_WIND

cd $PBS_O_WORKDIR

if [ $# -lt 8 ]; then
  printf "\nUsage: $0 <Input Dir> <Input Prefix> <Output Dir> <N recs> <DT> <Start Date> <N vars> <Var 1> [<Var 2> [<Var 3> [... <Var N>]]]\n\n"
else
  # read command line arguments
  InputDir=$1
  shift
  InputPrefix=$1
  shift
  OutputDir=$1
  shift
  Nrec=$1
  shift
  DT=$1
  shift
  StartDate=$1
  shift
  Nvar=$1
  shift
  DateCol=0
  TMAX=-1
  TMIN=-1
  for ((idx=1;idx <= $Nvar; idx++)) {
    VarNames[$idx]=$1
    if [ "$1" == "YEAR" ]; then
	DateCol=1
    fi
    if [ "$1" == "IN_TMAX" ]; then
	TMAX=$idx
    fi
    if [ "$1" == "IN_TMIN" ]; then
	TMIN=$idx
    fi
    shift
  }

  # check input path
  if [ ! -e $InputDir ]; then
    printf "ERROR: Input directory, $InputDir, does not exist!\n"
    exit 1
  fi

  # check output path
  if [ ! -e $OutputDir ]; then
    mkdir -p $OutputDir
  fi

  # update variable list for new columns
    # if daily data, add IN_TAVG
  if [ "$TMAX" -ge 0 ] && [ "$TMIN" -ge 0 ] && [ "$DT" -eq 24 ]; then # need to identify columns with IN_TMAX and IN_TMIN
      Nvar=$((Nvar+1))
      VarNames[$Nvar]="IN_TAVG"
  fi

    # if date information is not provided, then add it
  if [ $DateCol == 0 ]; then
      if [ $DT == 24 ]; then
	  for ((idx=$Nvar;idx > 0; idx--)) {
	      VarNames[$idx+3]=${VarNames[$idx]}
	  }
	  VarNames[1]=YEAR
	  VarNames[2]=MONTH
	  VarNames[3]=DAY
	  Nvar=$((Nvar+3))
      else
	  for ((idx=$Nvar;idx > 0; idx--)) {
	      VarNames[$idx+4]=${VarNames[$idx]}
	  }
	  VarNames[1]=YEAR
	  VarNames[2]=MONTH
	  VarNames[3]=DAY
	  VarNames[4]=DAY
	  Nvar=$((Nvar+4))
      fi
  fi

  # create output file header
  TmpHeader=`mktemp`
  printf "# NRECS: $Nrec\n" > $TmpHeader
  printf "# DT: $DT\n" >> $TmpHeader
  printf "# STARTDATE: $StartDate\n" >> $TmpHeader
  printf "# ALMA_OUTPUT: 0\n" >> $TmpHeader
  printf "# NVARS: $Nvar\n" >> $TmpHeader
  printf "#" >> $TmpHeader
  for ((idx=1;idx <= $Nvar; idx++)) {
    printf " ${VarNames[$idx]}"  >> $TmpHeader
  }
  printf "\n"  >> $TmpHeader

  # process all files in the selected directory
  for InFile in "$InputDir"/"$InputPrefix"*
  do
    name=${InFile##*/} # strip InputPath from filename
    OutFile=$OutputDir/$name
    echo Working on $name

    # create temporary working directory
    TmpDir=`mktemp -d ./TMPXXXXX`

    # set temporary file
    TmpFile=$InFile

    # if daily data, add IN_TAVG
    if [ "$TMAX" -ge 0 ] && [ "$TMIN" -ge 0 ] && [ "$DT" -eq 24 ]; then # need to identify columns with IN_TMAX and IN_TMIN
	awk '{ TAVG = ( $'$TMAX' + $'$TMIN' ) / 2.; printf( "%s %s\n", $0, TAVG ) }' $TmpFile > $TmpDir/TAVG
	TmpFile=$TmpDir/TAVG # reset tmp file
    fi

    # if date information is not provided, then add it
    if [ $DateCol == 0 ]; then
	if [ $DT == 24 ]; then
	    AddDateColumn.py $TmpFile $StartDate D "%Y %m %d" > $TmpDir/DATE
	else
	    AddDateColumn.py $TmpFile $StartDate "$DT"H "%Y %m %d %H" > $TmpDir/DATE
	fi
	TmpFile=$TmpDir/DATE # reset tmp file
    fi

    # write final version of file
    cat $TmpHeader $TmpFile > $OutFile

    # clean up temporary directory
    rm -rf $TmpDir

  done

  # clean up our mess
  rm $TmpHeader

fi

exit 0