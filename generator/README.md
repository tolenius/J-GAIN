# NN: The look-up table generator

The look-up table generator requires input (1) for the simulated chemical system, i.e. molecular cluster set, and (2) for the ambient conditions. As the table generator applies the ACDC cluster model, the molecular cluster input is same as for ACDC (https://github.com/tolenius/ACDC). Steps for setting up the input and running the table generator are detailed below.

The generator creates two output files: the look-up table and a descriptor file that contains information on the parameters and the table properties.

**Note**: The codes have been compiled and tested using gcc-9.2.0 and mpich-3.3.2 MPI implementation.

## Setting up the chemical system

### Molecular cluster data

The input files for the molecular cluster set and thermodynamics include the following:

* File listing the clusters included in the simulation
* File containing the cluster formation free energies, i.e. enthalpies and entropies
* In the case that charged species are included, file containing the cluster dipole moments and polarizabilities

These files are described in detail in the ACDC manual at https://github.com/tolenius/ACDC. It is also recommended to read the Quick Guide (in the same ACDC repository) to ensure that the cluster set is adequate and the simulation settings reasonable.

The input files can be saved in the subdirectory of the ACDC program: /src/ACDC_plugin_V02/Perl_input/.

**Note**: File acdc_simulation_setup.f90 includes some additional settings for the ACDC simulation. While they are by default set to be suitable for the table generator, it is good to note the following:

* For the look-up table application, the steady-state setting must always be on: solve_ss = .true.
* The vapor concentrations are by default set constant by using constant concentrations for vapors monomers, i.e. single vapor molecules. However, for very strongly clustering species, such as mixtures of sulfuric acid and strong amines, a significant fraction of vapor molecules may be clustered with one or more molecules of another compound. In this case, it is reasonable to define the vapor concentration as a sum of concentrations of free and clustered molecules in subroutine 'sources_and_constants' to avoid artificial overprediction of the formation rate. This setting can be used also for weakly clustering species for which it gives same result as using monomer concentrations, but note that the setting in 'sources_and_constants' is system-specific as it involves species names and thus less flexible.

### Generating the ACDC files

Before running the table generator, the ACDC subroutines need to be created by running run_perl.sh. The user needs to define the following variables in this script:

* The thermodynamics file and the possible dipole moment and polarizability file
* Names of the vapor species (as they appear in the ACDC input files)
* If the cluster set includes charged species, the 'l_incl_ions' variable must be set to 1 to automatically include ion-related settings.
* By default, the cluster set file name is assumed to be of the format 'input_AB_neutral.inp' or 'input_AB_neutral_neg_pos.inp', where A and B are the species names. The user can modify this by e.g. simply overwriting the 'cluster_file' variable.

When the variables are set, running

```console
$ cd /src/ACDC_plugin_V02/
$ ./run_perl.sh
$ cd ../..
```
creates the Fortran files that will be used by the table generator. The script also creates an ACDC makefile that will be automatically included when compiling the codes for the table generator (see below).

**Note**: By default, the temperature range is set only in the namelist for the table generator (see below), and 'variable_temperature' is set to 1 in run_perl.sh. A special case are tables where relative humidity (RH) is included in the parameters determining the formation rate. In this case, temperature and RH must be set in run_perl.sh as indicated in the script.

## Compilation

After creating the subroutines for the given cluster set, the table generator can be compiled in either serial or parallel (MPI) mode by using one of the two available makefiles: **makefile.serial** and **makefile.mpi**. The essential variables in the makefiles include the following:

**makefile.serial**: for serial compilation
```console
PROJECT_HOME = # Path to the working directory
```
**makefile.mpi**: for parallel compilation
```console
PROJECT_HOME = # Path to the working directory
PROC_NUM = # Number of MPI processes
```

The program is compiled as follows in the folder where the codes are located:

### Serial compilation

```console
$ ln -sf makefile.serial
$ make
```

### Parallel compilation by MPI

```console
$ ln -sf makefile.mpi
$ make
```

**Note**: This step will create all files in the path *$PROJECT_HOME/build* and the executable *lookup_gen.exe* in *$PROJECT_HOME*.

**Note**: OpenMP is not supported in this version.

### Cleaning

To clean all files generated in the compilation step, run the command below:
```console
$ make cleanall
```

## Setting the ambient conditions and other table properties in namelist.gen

Before running the compiled program, the values for the independent parameters are set in namelist.gen. The namelist file contains explanations for the settings, of which the most essential are:

* Number of independent parameters that determine the formation rate (including vapor concentrations and other ambient conditions)
* Names of the parameters: user-defined names for vapor species and 'T', 'CS', 'IPR', and 'RH' for temperature, cluster coagulation sink (as reference sink corresponding to H<sub>2</sub>SO<sub>4</sub>), atmospheric ion production rate, and relative humidity; 'T' and 'CS' are mandatory while 'IPR' and 'RH' are optional and must be left out when charged species are not included and when RH is not considered, respectively
* Minimum and maximum values (in SI units) and number of data points for each parameter
* Logical variables determining if the parameter is a vapor, and if the parameter values are placed evenly on a logarithmic instead of a linear scale; for logarithmic parameters the minimum and maximum values are given as the base-10 logarithm of the actual value

The namelist contains also the following:

* User-defined name and filepath of the output look-up table
* Units of all parameters to be saved in the descriptor file; SI units must be used for all parameters

## Running the code

After the compilation step the generator is run by
```console
$ make run
```
When the run is finished, the output will be in the directory defined in namelist.gen.
