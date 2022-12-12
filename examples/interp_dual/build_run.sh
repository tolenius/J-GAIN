#!/bin/bash 

set -e 
set -x 


SCRIPT_PATH=$(realpath $0)
SCRIPT_AD=$(dirname $SCRIPT_PATH)

FC=gfortran
FCFLAGS="-g -ffree-line-length-none -cpp -fcheck=all -ffpe-trap=invalid,zero,overflow -O3"



FILE=dual_table
export OBJS_AD=$SCRIPT_AD/build/obj
export BIN_AD=$SCRIPT_AD/build/bin
export INT_LIB_AD=$SCRIPT_AD/../../interpolator/build/lib
export INT_INC_AD=$SCRIPT_AD/../../interpolator/build/include




###############################

 

[ -d $OBJS_AD ] && rm -rf $OBJS_AD
[ -d $BIN_AD ] && rm -rf $BIN_AD
 


mkdir -p $OBJS_AD
mkdir -p $BIN_AD
 
$FC $FCFLAGS  $FILE.f90 -o $BIN_AD/$FILE.exe -L$INT_LIB_AD -I$INT_INC_AD  -J$OBJS_AD -linterp
 
$BIN_AD/$FILE.exe

 




