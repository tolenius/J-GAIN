# Look-up table interpolator

The interpolator is built by
```console
$ ./build-interpolator.sh
```
Examples of routines for calling the interpolator within a host model are given in the 'examples' folder. The user needs to (1) define the used table(s), (2) determine whether linear/logarithmic axes are used, and (3) set the input parameters related to the table(s) and decide how out-of-range values are treated. These settings are clarified below, and demonstrated in the example routines.

## Selecting the look-up tables and defining axes scales

The tables are loaded by an initialization routine (inputting the path and names of the table and descriptor files), as shown in the example programs. The user can define if parameters are interpolated along linear or logarithmic axes. This choice can improve the interpolation accuracy: we recommend using logarithmic scale for the formation rate, and for the independent parameters that are generally defined as logarithmic in the table generator namelist, namely the vapor concentrations, the cluster sink, and the ion production rate.

In the examples, linear/logarithmic interpolation is chosen according to 2 logical variables: one for the formation rate, and one for the (logarithmic) independent parameters.

**Note**: When parameters that are defined as logarithmic upon generation (i.e. placed evenly on a logarithmic scale) are set to be interpolated along logarithmic axis, their input values must be given as base10-logarithms of the actual values, instead of using the actual values.

## Setting the input parameters and treating out-of-range values

The interpolation routine is simply called by the values of the independent parameters related to the table(s), that need to be in the same order as they appear in the table descriptor file. By default, the molecular model implementation of the table generator assumes SI units for all parameters.

The interpolation returns the interpolated value for the formation rate, as well as information on the variables. Most importantly, the output includes a logical variable that tells if the input value of an independent parameter is within the table range, or below or above the limits. Also information on if which parameters correspond to vapor concentrations is included. These can be used to determine how to treat cases with value(s) outside the range. We recommend the following for out-of-range input:

* If vapor concentration value(s) is/are below the lower limit(s) of the table, set the formation rate to zero (the lower limits should also be set so that effectively no particle formation is expected at concentrations below them)
* Otherwise, use the interpolated value (obtained by fixing the value of the independent parameter to the minimum/maximum limit of the table)

## Example routines

All example programs in /examples/ apply simple, pre-generated (dummy) tables in /examples/tables/, and can be tested by running in the build script in each folder.

### Example of default implementation

The case in /examples/interp_dual/ demonstrates the application of two additive tables (that would in practice correspond to different chemical pathways). The case includes the treatment of axes scales and out-of-range values, as discussed above.

### Examples of linear/logarithmic interpolation

Folders /examples/interp_single/ and /examples/interp_single_log/ give examples of the details in choosing linear/logarithmic interpolation axes applying one table. Note that in actual applications, the recommended settings need to be used (see above).

## Note
When using the Intel compiler `ifort`, add the flag `-assume byterec` to ensure the code runs as expected.
