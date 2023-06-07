#!/bin/bash

set -x 

gfortran -g arbit_f.f90 -o  arbit_f.exe 
./arbit_f.exe
