## Generating tables from arbitrary data sources without the J-GAIN generator

This folder contains simple code examples on formatting and writing out J-GAIN compatible look-up tables from an arbitrary data source for the given descriptor file "lookup_3dep_4indep.desc". Code examples are given in Fortran and C++. The examples produce identical output files "lookup_3dep_4indep_cpp.bin" and "lookup_3dep_4indep_fortran.bin". The size of each file is 302940 bytes.

Note that while formation rate tables by default include one dependent variable (i.e. the formation rate), it is also possible to save more than one dependent variables for the same independent variables. We demonstrate this in the examples by including 3 dependent and 4 independent variables. In this case, the table includes an additional dimension, and must be saved considering the storing order: examples for column-major and row-major languages are provided by the Fortran and C++ codes, respectively.

* To run the Fortran code: bash run_f.sh
* To run the C++ code: bash run_c.sh 
* To check the output: bash run_check.sh | less 
