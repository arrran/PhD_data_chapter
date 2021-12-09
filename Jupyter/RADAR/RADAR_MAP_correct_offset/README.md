# info on the radar map files.


## process with highres and with line 8 too

01 (35 originally) interpolate line2line
02 grid it with background into a tiff
03 subtrack ice thickness from surface topo to get ice base

## for MITgcm smooth roof

05A
05B
05C
05D are copies of 35 36 37 3 but with smoothed interpolation to get output for MITgcm

## for hydrostatic equilibrium

04_where_float

## key files

are all the ones loaded in DATA/Jupyter/RADAR/RADAR_MAP
look at the plots to see whats what


### there has been some moving of files onto volumes, the later notebooks should have the accurate file locations


I redid the offset, earlier notebooks are at /Users/home/whitefar/DATA/Jupyter/RADAR/RADAR_MAP

I may have accidentally written over all the files so the wrong offset files may not even be in /Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP_wrong_offset
