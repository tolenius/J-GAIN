#!/bin/bash 

build_type=$1  # mpi serial 


if [ "x$build_type" != "xserial" ] && [ "x$build_type" != "xmpi" ]; then 
    echo please define build-type: serial or mpi 
	echo ./build-gen.sh build-type
    exit 	
fi 
   


SCRIPT_PATH=$(realpath $0)
GEN_DIR=$(dirname $SCRIPT_PATH)

export BUILd_DIR=$GEN_DIR/build.$build_type
export SRC_AD=$GEN_DIR/src
export OBJS_AD=$BUILd_DIR/obj
export BIN_AD=$BUILd_DIR/bin
export ACDC_LIB_AD=$GEN_DIR/../acdc/build/lib
export ACDC_INC_AD=$GEN_DIR/../acdc/build/include
MAKE_FILE=$GEN_DIR/scr/makefile.$build_type


###############################



[ -d $BUILd_DIR ] && rm -rf $BUILd_DIR



mkdir -p $OBJS_AD
mkdir -p $BIN_AD



make -f $MAKE_FILE


