# Molecular model set-up for generating tables from given input data

The table generator applies automatic molecular cluster modeling, implemented through the cluster model ACDC, as the combination of cluster dynamics simulations and qusntum chmistry input is a standard approach for formation rate calculations. In general, the cluster thermodynamics input can be requested from quantum chemistry modelers or found in publications. Questions on the input files and formats can be addressed to the ACDC or J-GAIN contacts.

## Setting up the chemical system with molecular cluster thermodynamics data

The input files for the molecular cluster set and thermodynamics include the following:

* File listing the clusters included in the simulation
* File containing the cluster formation free energies, i.e. enthalpies and entropies
* In the case that charged species are included, file containing the cluster dipole moments and polarizabilities

These files are described in detail in the ACDC manual at https://github.com/tolenius/ACDC. It is also recommended to read the Quick Guide (in the same ACDC repository) to ensure that the cluster set is adequate and the simulation settings reasonable.

The input files can be saved in the subdirectory of the ACDC implementation: /acdc/input/.

**Note**: File acdc_simulation_setup.f90 in /acdc/src/ includes some additional settings for the ACDC simulation. While they are by default set to be suitable for the table generator, it is good to note the following:

* For the look-up table application, the steady-state setting must always be on: solve_ss = .true.
* The vapor concentrations are by default set constant by using constant concentrations for vapors monomers, i.e. single vapor molecules. However, for very strongly clustering species, such as mixtures of sulfuric acid and strong amines, a significant fraction of vapor molecules may be clustered with one or more molecules of another compound. In this case, it is reasonable to define the vapor concentration as a sum of concentrations of free and clustered molecules in subroutine 'sources_and_constants' to avoid artificial overprediction of the formation rate. This setting can be used also for weakly clustering species for which it gives same result as using monomer concentrations, but note that the setting in 'sources_and_constants' is system-specific as it involves species names and thus less flexible.

## Generating the ACDC files

Before running the table generator, the ACDC subroutines need to be created by build-acdc.sh. The user needs to define the following variables in the beginning of the script:

* Names of the vapor species (as they appear in the ACDC input files)
* The thermodynamics file and the possible dipole moment and polarizability file
* If the cluster set includes charged species, the 'l_incl_ions' variable must be set to 1 to automatically include ion-related settings.
* By default, the cluster set file name is assumed to be of the format 'AB_neutral_suffix.inp' or 'AB_neutral_neg_pos_suffix.inp', where 'A' and 'B' are the species names, and '_suffix' is a user-defined suffix.

When the variables are set, running
```console
$ ./build-acdc.sh
```
creates the Fortran files that will be used by the table generator.

**Note**: By default, the temperature range is set only in the namelist for the table generator (see the instructions for the generator), and 'variable_temperature' is set to 1 in the ACDC build script. A special case are tables where relative humidity (RH) is included in the parameters determining the formation rate. In this case, temperature and RH must be set in the ACDC script as indicated within the script.
