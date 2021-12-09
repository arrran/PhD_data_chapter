# info on the radar map files.


## process with highres and with line 8 too

01 (35 originally) interpolate line2line
02 grid it with background into a tiff
03 subtrack ice thickness from surface topo to get ice base

## for MITgcm smooth roof

40
41
42
43 are copies of 35 36 37 39, but with smoothed interpolation to get output for MITgcm

## for hydrostatic equilibrium

44_where_float

## key files

are all the ones loaded in DATA/Jupyter/RADAR/RADAR_MAP
look at the plots to see whats what


### there has been some moving of files onto volumes, the later notebooks should have the accurate file locations


