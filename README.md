# Daniel
first build ACDC by runing 
```
# set the parameters  like 
# dip_file="dip_pol_298.15K_426clusters2016Apr25.txt"
# hs_file="HS298.15K_426clusters2016Apr25.txt"
# in build-acdc.sh the run
bash build-acdc.sh
```
build the generator 
```
bash build-gen.sh serial  # select between serial and mpi 
```
To run it see the example in examples/gen_mpi or examples/gen_serial


build the interpolator by running 
```
bash build-interpolator.sh
```
in the examples folder there are many examples:
* gen_mpi: run the generator in parallel mode usning mpi 
* gen_serial: run the generator in serial
* interp_dual: interolation using tow tables as illustrated in the paper
* interp_single: simple interoplation using one table and using the actual table without modification
* interp_single_log: interoplation using one table and using the actual table without modification but in the log space for both table values and log parameters.
* tables: two samples of pre-generated tables to be used in the examples.



# NN &ndash; New-particle formation rate look-up table generator and interpolator

NN is a tool for incorporating aerosol formation rate predictions from molecular modeling to large-scale atmospheric models using look-up tables. NN consists of two parts: a table generator and a table interpolator.

* NN can be applied for any chemical system with arbitrary number of species, and with user-defined selection of ambient parameters that determine the formation rate.
* The table generator takes quantum chemical molecular cluster thermodynamics data as input, and applies the cluster population model ACDC to solve the new-particle formation rate for given ambient conditions.
* Parameter ranges and table resolution are defined by the user.
* The table interpolator uses the given look-up tables and multivariate interpolation to determine the formation rate at wanted conditions.
* The interpolator is easily coupled to a host model, and can be applied for several alternate or additive tables for particle formation for different chemical systems.

Instructions for applying the generator and the interpolator can be found in their respective folders.

## Prerequisites

In order to use NN (Unix or Windows), the following software must be installed:

* Fortran
* Perl (only for the table generator)

## Citation

If you use the codes provided in this repository in any study, please cite the following paper:

Yazgi, D., and Olenius, T.: TITLE, Journal XX, YY-YY (2022), doi

## License

This project is licensed under the terms of the GNU General Public License v3.0, as provided with this repository.

The table generator also includes the following open-source routines: the molecular cluster population dynamics program Atmospheric Cluster Dynamics Code (ACDC; GPL 3), and the Fortran VODE solver (In the Public Domain; IPD). For these programs, the original license or acknowledgement is found in the same folder where the program source code is located and must not be separated from the source code.

## Authors

Code development: **Daniel Yazgi** (daniel.yazgi@smhi.se)

Project PI: **Tinja Olenius** (tinja.olenius@alumni.helsinki.fi)
