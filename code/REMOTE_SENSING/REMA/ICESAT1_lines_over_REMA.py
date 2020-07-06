#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 09:54:26 2020

@author: arran

This script gets the elevations from REMA lines for points along the icesatlines

icesat lines are read in as shp, elevations are added, then outputs shp files

the shapefiles_of_icesat1_over_channel files of the lines are made at ICESAT1/ICESAT1_fixing_blocky_residual_problems.ipynb
these use the coordinates made using save_zps() from icesat1_save_df_of_lines_residuals.py

NB extrapolation is done by x for 0211 and y for 0099, v important, must redo REMA projection if changing this
"""

import rasterio as rio
import rasterio.mask
import fiona
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np
import os
import sys
import glob
import matplotlib.pyplot as plt
from scipy import interpolate

from shapely.geometry import LineString


gis_path ="/Users/home/whitefar/DATA/REMOTE_SENSING/ICESAT1/shapefiles_of_icesat1_over_channel/"

REMA_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/'
#REMA_files_paths = glob.glob(os.path.join(REMA_filepath,"**.tif"),recursive=True)

#indicies of REMA strips which intersect the study area + buffer, ie larger than channel area.
#this file is made in find_crop_REMA, 
indicies_which_intersect = np.loadtxt("/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/indicies_which_intersect.txt").astype(int).tolist()



lines_files_paths = glob.glob(os.path.join(gis_path,"**.shp"),recursive=True)
lines_names = [os.path.splitext(os.path.split(line_file_path )[1])[0] for line_file_path in lines_files_paths]


gdf = gpd.read_file('/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1.shp',crs="EPSG:3031")
gdf['nid']=gdf.index


REMA_shapes_df =gdf.iloc[indicies_which_intersect]
REMA_shapes_df.reset_index(drop=True,inplace=True)


# Check for intersection of line and REMA, then if they intersect, write the elevations of the REMA to the icesatline shp file
for i, line_file_path in enumerate(lines_files_paths):
    
    icesat_line = gpd.read_file(line_file_path)
    
    #print("sampling elevations on "+lines_names[i])
        
    for s, REMA_shape in enumerate(REMA_shapes_df.name):
        
        # check for intersection
        if not REMA_shapes_df.geometry.iloc[s].intersects( LineString(icesat_line.geometry.tolist()) ):
            #print("no intersection with "+REMA_shape)
            continue
        
        #print("yes, intersection with "+REMA_shape)
                
        tiff_stripe_fname = REMA_shape + "_dem.tif"

        #write REMA elevations to radarline                
        with rio.open(REMA_filepath + tiff_stripe_fname) as src:
            
            coords = [(x,y) for x, y in zip(icesat_line.geometry.x, icesat_line.geometry.y)]
            
            elevations = [elevation[0] for elevation in src.sample(coords)]
        
        #column_name =f"i{indicies_which_intersect[s]}date{REMA_shape.split('_')[2]}"
        
        column_name =f"RE_{REMA_shapes_df.iloc[s].acquisitio}" #this can be changed to "nid_{REMA_shapes_df.iloc[s].nid}" and then use df to get more info on the REMA strip
        icesat_line[column_name] = pd.Series(elevations).replace(-9999.0, np.nan)
        
        #print(f"elevations printed to line for REMA on {REMA_shapes_df.acquisit_1.iloc[s]}")
        
        print(f"{s}/{len(REMA_shapes_df)} of way through REMA strip")
        
        del tiff_stripe_fname, coords, elevations
    
    print(f"{i}/{len(lines_files_paths)} of way through lines")
    
    icesat_line.to_file(gis_path+lines_names[i]+".shp")
    
    
    del icesat_line
    
# =============================================================================

# write a dictionary which associates each line with REMA strips, 

lines_dict_name = {}
lines_dict_nid = {}
lines_dict_date = {}

for i, line_file_path in enumerate(lines_files_paths):
    
    icesat_line = gpd.read_file(line_file_path)
    
    REMAnid = []
    REMAname = []
    REMAdate = []
        
    for s, REMA_shape in enumerate(REMA_shapes_df.name):
        
        if not REMA_shapes_df.geometry.iloc[s].intersects( LineString(icesat_line.geometry.tolist()) ):
            continue
        
        REMAnid.append(f"nid_{REMA_shapes_df.iloc[s].nid}")    
        REMAname.append(REMA_shape)
        REMAdate.append(REMA_shapes_df.iloc[s].acquisitio)
        
    
    print(f"{i}/{len(lines_files_paths)} of way through lines")
    
    lines_dict_nid[lines_names[i]] = REMAnid
    lines_dict_name[lines_names[i]] = REMAname
    lines_dict_date[lines_names[i]] = REMAdate

with open(gis_path+'REMAnid_over_icesat1lines.txt','w') as f:
    f.write(str(lines_dict_nid))
with open(gis_path+'REMAdate_over_icesat1lines.txt','w') as f:
    f.write(str(lines_dict_date))
with open(gis_path+'REMAname_over_icesat1lines.txt','w') as f:
    f.write(str(lines_dict_name))
    
    
# =============================================================================

#PLOT

# #open the dictionary associating each line with REMA strips
# with open(gis_path+'REMAnid_over_icesatlines.txt','r') as f:
#     ld = eval(f.read())


# line_index = [x for x in zip(lines_names,range(len(lines_names)))]

# #print the line names with an index
# print(list(zip(range(len(lines_names)),lines_names)))

# #with line from first to last point removed
# i=14  #which radarline



# def plot_line(line_name, df=df,legend='date', remove_trend=True,lr=False):
#     """
#     Plots the icesat with REMA elevations.
#     each elevation line has a linear line between first and last points removed, to try and reduce tide effects.
#     only use this on lines perpendicular to channel
#     """
#     leg = []
#     rl = gpd.read_file(gis_path+line_name+".shp")
#     plt.figure(figsize=(10,7))
#     for REMA in ld[line_name]:
        
#         nid = int(REMA.split('_')[1])
#         date = df.loc[nid].acquisitio
#         #whether to flip line lr
#         if (rl.geometry.x.iloc[0] - rl.geometry.x.iloc[-1]) < 0:
#             lr=False
#         else:
#             lr=True
        
#         if lr==True:
#             remaline =  rl[REMA].iloc[::-1]
#         else:
#             remaline =  rl[REMA]
#         if remove_trend==True:
#             #make line from first to last point
#             f = interpolate.interp1d( [rl.distan_cum.iloc[0],rl.distan_cum.iloc[-1]], [rl[REMA].iloc[0],rl[REMA].iloc[-1]])
#             trendline = f(rl.distan_cum)
#             if lr==True:
#                 plt.plot(rl.distan_cum,remaline-trendline[::-1])
#             else:
#                 plt.plot(rl.distan_cum,remaline-trendline)
                    
#         else:
#             plt.plot(rl.distan_cum,remaline)
#         plt.title(line_name)
#         if legend=='date':
#             leg.append(date)
#         elif legend=='nid':
#             leg.append(nid)
        
        
#     plt.legend(leg) 
#     plt.grid(True)
#     plt.xlabel("distance, 'm'")
# #     plt.xlim([-200,4500])
# #     plt.ylim([-18,2])
#     plt.show()
  








# #THIS MAKES A GEODATAFRAME WITH A ROW FOR EACH RADARLINE
# lines = []
# names = []

# for i, line_file_path in enumerate(lines_files_paths):
#     df = gpd.read_file(line_file_path)
#     lines.append( LineString(df.geometry.tolist()) )
#     names.append(lines_names[i])
#     del df
#     print(i / len(lines_files_paths))


# lines_gdf = gpd.GeoDataFrame({'shortname': names}, geometry=lines,crs="EPSG:3031")
# #lines_gdf.to_file(gis_path+"lines_linestrings_gpkg")





