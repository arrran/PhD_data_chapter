# info on the radar map files.

## radar data organisation

45 redo radar correct offset. Was originally 0_RADAR_bed_conglo...

## interpolation
26_RADAR_interpolate_line2line_by_distance_downchan is the money interpolation

it is then gridded at  20_RADAR_interpolate_channel_background_separate_wholechannel.ipynb

surface is interpolated later at  /31_RADAR_grid_make_tiff_of_ice_base.ipynb and used to find base of ice

## redo as highres and with line 8 too

35 interpolate line2line
36 grid it with background into a tiff
37 subtrack ice thickness from surface topo to get ice base

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


