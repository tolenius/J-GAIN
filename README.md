# J-GAIN: Aerosol particle formation rate look-up table generator and interpolator

J-GAIN is a tool for incorporating aerosol formation rate predictions from molecular modeling to large-scale atmospheric models using look-up tables. J-GAIN consists of two parts: a table generator and a table interpolator.

* J-GAIN can be applied for any chemical system with arbitrary number of species, and with user-defined selection of ambient parameters that determine the formation rate.
* The table generator takes as input quantum chemical molecular cluster thermodynamics data, and applies the cluster population model ACDC to solve the new-particle formation rate for given ambient conditions.
* Parameter ranges and table resolution are defined by the user.
* The table interpolator uses the given look-up tables and multivariate interpolation to determine the formation rate at wanted conditions.
* The interpolator is easily coupled to a host model, and can be applied for several alternate or additive tables for particle formation for different chemical systems.

This repository contains 4 main folders:

1. acdc: &emsp; &emsp; &emsp; Molecular model set-up for generating new tables from given input data
2. generator: &ensp; &ensp; Generation of tables based on the model set-up in 'acdc'
3. interpolator: &ensp; Interpolation routines
4. examples: &emsp; &nbsp; Examples for table generation settings, and for calling the interpolator for arbitrary tables

Instructions for each step can be found in the respective folders.

## Prerequisites

In order to use J-GAIN, the following software must be installed:

* Fortran
* Perl (only for the table generator)

The current versions have been tested with gfortran of gcc (GCC) 6.4.0. and openmpi (for parallel compilation of the table generator).

## Citation

If you use the codes provided in this repository in any study, please cite the following works:

* This repository (https://github.com/tolenius/J-GAIN)
* Yazgi, D., and Olenius, T.: J-GAIN v1.0: A flexible tool to incorporate aerosol formation rates obtained by molecular models into large-scale models, submitted to Geosci. Model Dev. Discuss. (2022)

## License

This project is licensed under the terms of the GNU General Public License v3.0, as provided with this repository.

The table generator also applies the following open-source routines: the molecular cluster population dynamics program Atmospheric Cluster Dynamics Code (ACDC; GPL 3), and the Fortran VODE solver (In the Public Domain; IPD). For these programs, the original license or acknowledgement is found in the same folder where the program source code is located and must not be separated from the source code.

## Authors

Code development: **Daniel Yazgi** (daniel.yazgi@smhi.se)

Project PI: **Tinja Olenius** (tinja.olenius@alumni.helsinki.fi)
