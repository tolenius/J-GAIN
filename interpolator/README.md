# NN: The look-up table interpolator

The look-up table interpolator can be readily coupled to a host model through a provided example routine. The user needs to (1) define the used table(s), and (2) set the input parameters related to the table(s).

**Note**: The codes have been compiled and tested using **xxx: WHICH GFORTRAN VERSION?**

## Selecting the look-up tables

The tables are selected by filling in the directory path and names of the table and descriptor files in namelist.nml.

The namelist also includes a logical variable for determining if the formation rate is interpolated on a logarithmic instead of a linear scale; by default this is set to  logarithmic (.true.) as this generally gives more accurate results. (An independent parameter is interpolated on a logarithmic scale if the parameter is defined as 'logarithmic'; see the instructions for the table generator.)

**xxx: How are several tables defined?**

**xxx: How can alternate or additive tables be defined? We need to include an example of this.**

## Setting the input parameters

The example program demonstrates loading the table(s) and calling the interpolator for an example case of particle formation from H<sub>2</sub>SO<sub>4</sub> and NH<sub>3</sub>. The program contains an example of a subroutine with the independent parameters as input and the interpolated formation rate as output. Here, the input parameters 'A', 'N', 'CS', 'T', and 'IPR' correspond to [H<sub>2</sub>SO<sub>4</sub>], [NH<sub>3</sub>], cluster scavenging sink, temperature, and ion pair production rate.

The user needs to modify the input parameters to correspond to the used table(s). All input and output are in SI units.

## Treating values outside the table

The treatment of cases with independent parameter values that fall outside the range of the look-up table is determined by variable 'beyond_limits_policy' in the example program. This array contains a value for each input parameter.

Currently implemented values are:

* 0 - The value of the independent parameter is set to the minimum/maximum limit of the table when the input value is below/over the table limits
* 1 - Otherwise as 0, but zero J is returned when the parameter is below the lower limit of the table; this option should be used for vapor concentrations and the lower limits set so that effectively no particle formation is expected at concentrations below the limits
