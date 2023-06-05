# Examples for applying the generator and the interpolator

This folder contains the following examples:

* gen_serial: run the table generator in serial
* gen_serial_rh: run the table generator in serial with RH included as an independent variable. Build acdc with variable_temperature=0
* gen_mpi: run the table generator in parallel using mpi
* interp_dual: example of a practical implementation of the interpolator using two separate, additive tables (see also the documentation in 'interpolator')
* interp_single: additional example of simple interpolation using linear interpolation axes
* interp_single_log: additional example of simple interpolation using logarithmic interpolation axes for both the table values and the input parameters defined as logarithmic
* tables: two simple, pre-generated tables to be used in the examples
