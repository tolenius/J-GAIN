# Look-up table generator

The look-up table generator requires input (1) for the modeled chemical system, i.e. molecular cluster set, and (2) for the ambient conditions. The molecular cluster input is set in the 'acdc' folder as detailed in the instructions within the folder. The ambient conditions, i.e. the parameter space over which formation rates are computed, are set in the generator namelist file as instructed below.

The generator creates two main output files that are used by the interpolator: the look-up table and a descriptor file that contains information on the parameters and the table properties.

***Note*** While the table generator is primarily designed for formation rates obtained by molecular cluster dynamics modeling, it is possible to use other data sources by saving similar tables in the framework of another model (see 'Other data sources' below).

## Building the generator

After creating the ACDC subroutines, the table generator can be built in either serial or parallel mode:
```console
$ ./build-gen.sh serial
```
or
```console
$ ./build-gen.sh mpi
```

## Setting the ambient conditions and other table properties, and running the generator

Before running the compiled program, the values for the independent parameters are set in namelist.gen. The namelist file contains explanations for the settings, of which the most essential are:

* Number of independent parameters that determine the formation rate (including vapor concentrations and other ambient conditions)
* Names of the parameters: user-defined names for vapor species and 'T', 'CS', 'IPR', and 'RH' for temperature, cluster coagulation sink (as reference sink corresponding to H<sub>2</sub>SO<sub>4</sub>), atmospheric ion production rate, and relative humidity; 'T' and 'CS' are mandatory while 'IPR' and 'RH' are optional and must be left out when charged species are not included and when RH is not considered, respectively
* Minimum and maximum values (in SI units) and number of data points for each parameter
* Logical variables determining if the parameter is a vapor, and if the parameter values are placed evenly on a logarithmic instead of a linear scale; for logarithmic parameters the minimum and maximum values are given as the base-10 logarithm of the actual value

The namelist contains also the following:

* User-defined name of the output look-up table
* Units of all parameters to be saved in the descriptor file; SI units must be used for all parameters

The user needs to ensure that the order of the independent parameters is consistent in all related arrays (i.e. names, dimensions, units, ...; the exact order doesn't matter as long as the arrays are consistent).

Examples of namelists can be found in /examples/gen_serial/ and /examples/gen_mpi/. When all input is set, the generator can be run (in /examples/gen_serial/ and /examples/gen_mpi/) by
```console
$ ./run.sh
```
or by submitting the run as a batch script. The output will be in the directory defined in the namelist.

## Other data sources

Saving tables from other models or data sources can be done as follows:

1. The information on the independent parameters and their ranges need to be saved in a descriptor file, examples of which can be found in 'examples/tables' and used as a template. The file contains the following:
    * Number, names and units of dependent parameters (here trivially one parameter, i.e. the formation rate)
    * Number of independent parameters
    * Lists for the following settings, where the list elements correspond to the independent parameters:
        * Number of values
        * Name
        * Units
        * Minimum and maximum values (note that for parameters defined as logarithmic, the limit values are given as base-10 logarithm of the actual value; see below)
        * Logical telling if the parameter is defined as logarithmic (i.e. its values will be placed evenly on a logarithmic scale)
        * Logical telling if the parameter is a vapor
    * Name of the binary file containing the formation rates
    * Total number of data points in the binary file (i.e. the product of the numbers of values for the independent parameters)

2. The formation rates need to be saved in a binary file as a flattened 1-dimensional array. If the rate depends on *N* independent parameters that are listed in the descriptor file in the order *x*<sub>1</sub>, *x*<sub>2</sub>, ..., *x*<sub>*N*</sub>, the flattened array is created by iterating through the independent parameters in the following order:

$i$=0\
for $i_1$=1 to $i_{\text{max},1}$\
&emsp; ...\
&emsp; for $i_N$=1 to $i_{\text{max},N}$\
&emsp; &emsp; &emsp; $i=i+1$\
&emsp; &emsp; &emsp; Save $f(i)=f(x_{1,i_1}, x_{2,i_2}, ..., x_{N,i_N})$\
&emsp; end\
&emsp; ...\
end

Simple code examples for constructing the loop and saving the table file are provided in examples/arbitrary_table. The examples include Fortran and C++, but a similar structure can be included in any existing model code.
