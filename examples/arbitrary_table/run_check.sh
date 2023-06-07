#!/bin/bash

set -x 

gfortran -g check.f90 -o  check.exe 
./check.exe
