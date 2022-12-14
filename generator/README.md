# Look-up table generator

The look-up table generator requires input (1) for the modeled chemical system, i.e. molecular cluster set, and (2) for the ambient conditions. The molecular cluster input is set in the 'acdc' folder as detailed in the instructions within the folder. The ambient conditions, i.e. the parameter space over which formation rates are computed, are set in the generator namelist file as instructed below.

The generator creates two main output files that are used by the interpolator: the look-up table and a descriptor file that contains information on the parameters and the table properties.

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
