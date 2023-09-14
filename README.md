# J-GAIN: Aerosol particle formation rate look-up table generator and interpolator

J-GAIN is a tool for incorporating aerosol formation rate predictions from molecular modeling to large-scale atmospheric models using look-up tables. J-GAIN consists of two parts: a table generator and a table interpolator. By default, the table generator performs automatic molecular cluster dynamics modeling with quantum chemical input data, as this the standard method for formation rate predictions, although also other data sources are possible.

* J-GAIN can be applied for any chemical system with arbitrary number of species, and with user-defined selection of ambient parameters that determine the formation rate.
* The table generator and interpolator are two separate parts that can be used independently:
    * The table generator saves multidimensional formation rate data together with table information.
        * The generator takes as input quantum chemical molecular cluster thermodynamics data, and applies the ACDC cluster dynamics solver to determine the new-particle formation rate for given ambient conditions.
        * Parameter ranges and table resolution are defined by the user.
        * Also other modeling approaches or data sources can be used by saving the data in the given table format.
    * The table interpolator uses the given look-up tables and multivariate interpolation to determine the formation rate at wanted conditions.
        * The interpolator is easily coupled to a host model, and can be applied for several alternate or additive tables for particle formation for different chemical systems.
        * Tables can be trivially added or removed, as long as the table format is correct.

This repository contains 5 main folders:

1. acdc\
    Molecular model set-up for generating new tables from given input data by molecular cluster dynamics modeling
2. generator\
    Generation of tables based on the model set-up in 'acdc', and instructions for data sources other than cluster modeling
3. interpolator\
    Interpolation routines
4. examples\
    Examples for table generation settings, and for calling the interpolator for arbitrary tables
5. GMD_example\
    Input and step-by-step instructions for reproducing a test case for H<sub>2</sub>SO<sub>4</sub>-NH<sub>3</sub> chemistry data presented by Yazgi and Olenius, 2023, doi:10.5194/egusphere-2022-1464

## Usage

Information on each part can be found in the respective folders. The molecular cluster thermodynamics input can be requested from quantum chemistry data providers, if it is not directly available in publications.

Briefly, the whole procedure from table generation to retrieval of interpolated values is conducted as follows:

* Table generation by the default approach (otherwise follow the instructions for other data sources in the 'generator' folder):
    1. In the 'acdc' folder, set the molecular cluster input as instructed, and execute the script that builds the cluster model set-up
    2. In the 'generator' folder, set the ranges of ambient conditions and execute the scripts that run the table generator
* Interpolation:
    1. In the 'interpolator' folder, execute the script that builds the interpolator
    2. In the 'examples' folder, use the example scripts, e.g. interp_dual/dual_table.f90, for implementing the interpolator
    3. In the script, set the table path(s) and the input for the ambient conditions, and execute the program

The example case in folder 'GMD_example' provides a detailed guide for the usage.

## Prerequisites

In order to use J-GAIN, the following software must be installed:

* Fortran
* Perl (only for table generation by molecular cluster modeling)

The current versions have been tested with gfortran of gcc (GCC) 6.4.0. and openmpi (for parallel compilation of the table generator).

## Citation

If you use the codes provided in this repository in any study, please cite the following works:

* This repository (https://github.com/tolenius/J-GAIN)
* Yazgi, D., and Olenius, T.: J-GAIN v1.1: a flexible tool to incorporate aerosol formation rates obtained by molecular models into large-scale models, *Geosci. Model Dev.* 16, 5237–5249 (2023), https://doi.org/10.5194/gmd-16-5237-2023

## License

This project is licensed under the terms of the GNU General Public License v3.0, as provided with this repository.

The table generator also applies the following open-source routines: the molecular cluster population dynamics program Atmospheric Cluster Dynamics Code (ACDC; GPL 3), and the Fortran VODE solver (In the Public Domain; IPD). For these programs, the original license or acknowledgement (if not the same as for J-GAIN) is found in the same folder where the program source code is located and must not be separated from the source code.

## Authors

Code development: **Daniel Yazgi** (daniel.yazgi@smhi.se)

Project PI: **Tinja Olenius** (tinja.olenius@alumni.helsinki.fi)

