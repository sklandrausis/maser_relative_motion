# MRM - Maser relative motion

## Capabilities of MRM
MRM is complex of python scripts that allow analysis of proper motion and evolution of long-surviving cloudlets using EVN maser observation data. As input, it takes two types of files data files created by AIPS task ... and ISPEC files AIPS task ...

## Configuration of MRM
MRM consist of two configuration files: 1) config.cfg, 2) plot.cfg, both are located in the directory config. Configuration file plot.cfg have only one section main that contain matplotlib configuration see more in https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html. Configuration file config.cfg have these sections paths, parameters, groups. The paths section contain all the data input and output paths. The parameters section are a collection of hardcoded parameters used in the data processing. The groups it is used for Gauss fitting.

| **Paths** | **Description** |
| --- | --- |
| dataFiles | AIPS output|
| relative_motion | input and output files for relative motion scripts|
| ispec | ISPEC|
| cloudlet |  cloudlet Gauss fit statistics | 
| groups |  group files for  Gauss fit| 

| **Parameters** | **Description** |
| --- | --- |
| fileOrder | order of input files|
| dates | dates of observation|
|gauss | which groups double gauss are used in which observation|

## Analysis of proper motion

| **Scripts** | **Description** |
| --- | --- |
| global_groups.py | Takes AIPS output and create global groups and save result in file relative_motion/relative_motion_group.dat |
| mean_motion.py | Takes script global_groups.py output and compute mean motion result is stored in two files output_mean_motion.dat and position_angle_motion_linearity.dat|
| linear_errors.py | Takes script mean_motion.py output and compute linear errors result is stored in three files linearity_errors_fitted_cm.dat, linearity_errors_fitted_tex_cm.dat, linearity_errors_fitted_tex_sort.dat|
| position_angle.py | Takes script linear_errors.py output and compute position angle result is stored|

## Evolution of long-surviving cloudlets
| **Scripts** | **Description** |
| --- | --- |
|total_spectrum.py | Takes AIPS output dataFiles and ispec file and plot total spectrum and map for all observations |
| select_groups_manual.py | Allow to create cloudlet groups for observation and save result in groups directory|
| gauss_fit_for_all_groups.py | Fits Gauss components for all groups and prints LATEX table for Gauss fit parameters |
| gauss_fit_for_multiple_groups.py | Fits Gauss components for selected groups. Take a input list of groups |
| gauss_fit_for_single_group.py | Fits Gauss components for selected group. Take a input group number |
| group_3d.py | Plot 3D plot for all groups x axis - RA, y axis - DEC, z axis - velocity  |
| gauss_fit_and_position_angle.py | Fits Gauss components for selected group and compute position angle. Store results in two files cloudlet_<group number>._coords.csv and cloudlet_<group number>._sats.csv"  |
| gauss_fit_and_position_angle2.py | Same as gauss_fit_and_position_angle.py |
| gauss_fit_and_position_angle_for_single_epoch.py | Same as gauss_fit_and_position_angle.py but for single observation |
|vel_fit_vs_time.py | Plot velocity fit vs time |


## Getting Help

Bug reports, feature requests and make contributions (e.g. code patches) can be reported by opening a &quot;new issue&quot; ticket on GitHub. Please give as much information (e.g. the software component, version) as you can in the ticket. For bugs, it is extremely useful if a small self-contained code snippet that reproduces the problem is provided.

## Acknowledgements
This software was written by Jānis Šteinbergs under the supervision of Artis Aberfelds. If you make use of this software to get results that appear in a publication or presentation please include this acknowledgement: &quot;We have made use of MDPS, a tool developed by Janis Steinbergs.&quot;

This work was supported by Latvian Council of Science Project “Research of Galactic Masers” Nr.: lzp-2018/1-0291.
