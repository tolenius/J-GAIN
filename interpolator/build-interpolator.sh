#!/bin/bash 

set -e 
set -x 


SCRIPT_PATH=$(realpath $0)
INT_DIR=$(dirname $SCRIPT_PATH)

FC=gfortran
FCFLAGS="-g -ffree-line-length-none -cpp -fcheck=all -ffpe-trap=invalid,zero,overflow -O3"



export SRC_AD=$INT_DIR/src
export LIB_AD=$INT_DIR/build/lib
export INC_AD=$INT_DIR/build/include
export OBJS_AD=$INT_DIR/build/obj



###############################

 

[ -d $OBJS_AD ] && rm -rf $OBJS_AD
[ -d $LIB_AD ] && rm -rf $LIB_AD
[ -d $INC_AD ] && rm -rf $INC_AD


mkdir -p $OBJS_AD
mkdir -p $LIB_AD
mkdir -p $INC_AD

files=(mo_debug
mo_kind
mo_runtime
mo_table_lookup
mo_serial_lookup
mo_utils
)



for F in ${files[@]}; do 
   $FC $FCFLAGS -c $SRC_AD/$F.f90 -o $OBJS_AD/$F.o  -I$INC_AD  -J$INC_AD 
done 


ar cr $LIB_AD/libinterp.a $OBJS_AD/*.o 






