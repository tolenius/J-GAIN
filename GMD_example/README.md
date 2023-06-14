# Step-by-step instructions with an H<sub>2</sub>SO<sub>4</sub>-NH<sub>3</sub> example case

The detailed instructions below clarify the whole procedure from preparing and running the table generator with given molecular cluster thermochemistry input data to retrieving and evaluating interpolated values. The H<sub>2</sub>SO<sub>4</sub>-NH<sub>3</sub> example case presented in the work by Yazgi and Olenius (*Geosci. Model Dev.* 2023) is used for demonstration, and the same instructions can be applied for other molecular cluster input data and/or other look-up tables.

The first section below summarizes the usage of an automatic template for completing the steps. In addition, to clarify the essential procedures, the next section explicitly lists the commands for the main steps. For these instructions, the J-GAIN repository must be downloaded as is including all subfolders, with no changes in the structure.

## Automatic template

The script run_me.sh can be used to complete all steps from table generation to interpolation and evaluation of the interpolated values. The script generates two tables according to the molecular cluster input data in the 'acdc' folder and the namelists in the 'GMD_example/nam' subfolder: a coarser table (experiment case) and a higher-resolution table that gives the exact reference values (reference case), and produces a figure corresponding to Figure 4 by Yazgi and Olenius (*Geosci. Model Dev.* 2023). The figure shows the error in the values interpolated from the coarser table, obtained by comparison to the reference table.

The script is run by
```console
$ ./run_me.sh [option]
```
where the option refers to the steps to be executed. Using option '8' executes all steps; to see other options run the script without any input for the option.

The template is intended as an example for creating and comparing tables. While the interpolation is here performed over the \[H<sub>2</sub>SO<sub>4</sub>\] and \[NH<sub>3</sub>\] axes for the H<sub>2</sub>SO<sub>4</sub>-NH<sub>3</sub> example case, similar tests can be conducted for other model cases by examining and modifying run_me.sh, the namelists in the 'nam' subfolder and the scripts in the 'src' subfolder. The user can explore the interpolation and plotting routines in 'src' and adapt the scripts for their needs.

## Commands for the main steps

The commands for each main step are listed below; for more information, see the documentation in the respective folders.

1. For given molecular cluster thermochemistry input files, the first step is to create the cluster model. For this, go to the 'acdc' folder and place the input in the 'input' subfolder; the example input files included in the repository correspond to the H<sub>2</sub>SO<sub>4</sub>-NH<sub>3</sub> example case. In 'acdc', run
```console
$ ./build-acdc.sh
```

2. To build the table generator (here in the serial mode), go to 'generator' and run
```console
$ ./build-gen.sh serial
```

3. To generate tables for the H<sub>2</sub>SO<sub>4</sub>-NH<sub>3</sub> case, use the namelist files in 'GMD_example/nam' (filename extension '.gen'; for evaluating interpolated values, two tables can be generated as in the automatic template discussed above). Go to 'examples/gen_serial', copy or link the namelist to the file namelist.gen, and run
```console
$ ./run.sh
```

4. To build the table interpolator, go to 'interpolator' and run
```console
$ ./build-interpolator.sh
```

5. Both the tables and the interpolation routines are now ready. For an example of retrieving interpolated values from arbitrary tables, see e.g. 'examples/interp_dual'. Insert correct look-up table file names and conditions at which interpolated values will be determined in the example script dual_table.f90, and compile and run by
```console
$ ./build_run.sh
```

6. The user can also construct own templates for iterating over given values. For example, for table comparisons, the values can be set according to the values of a reference table (i.e. by starting from the minimum values and proceeding according to the reference table intervals). This is done also in the automatic script 'GMD_example/src/table_to_table_interpolator.f90'.
