#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 09:54:26 2020

@author: arran

This script gets the elevations from REMA lines for points along the radarlines

Radarlines are read in as gpkgs, elevations are added, then outputs shp files
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


gis_path ="/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RES/PROCESSED_LINES_GISFILE/"

REMA_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/'
#REMA_files_paths = glob.glob(os.path.join(REMA_filepath,"**.tif"),recursive=True)

indicies_which_intersect = np.loadtxt("/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/indicies_which_intersect.txt").astype(int).tolist()



lines_files_paths = glob.glob(os.path.join(gis_path,"**.gpkg"),recursive=True)
lines_names = [os.path.splitext(os.path.split(line_file_path )[1])[0] for line_file_path in lines_files_paths]


gdf = gpd.read_file('/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1.shp',crs="EPSG:3031")
gdf['nid']=gdf.index


REMA_shapes_df =gdf.iloc[indicies_which_intersect]
REMA_shapes_df.reset_index(drop=True,inplace=True)


# Check for intersection of line and REMA, then if they intersect, write the elevations of the REMA to the radarline
for i, line_file_path in enumerate(lines_files_paths):
    
    radar_line = gpd.read_file(line_file_path)
    
    #print("sampling elevations on "+lines_names[i])
        
    for s, REMA_shape in enumerate(REMA_shapes_df.name):
        
        # check for intersection
        if not REMA_shapes_df.geometry.iloc[s].intersects( LineString(radar_line.geometry.tolist()) ):
            #print("no intersection with "+REMA_shape)
            continue
        
        #print("yes, intersection with "+REMA_shape)
                
        tiff_stripe_fname = REMA_shape + "_dem.tif"

        #write REMA elevations to radarline                
        with rio.open(REMA_filepath + tiff_stripe_fname) as src:
            
            coords = [(x,y) for x, y in zip(radar_line.geometry.x, radar_line.geometry.y)]
            
            elevations = [elevation[0] for elevation in src.sample(coords)]
        
        #column_name =f"i{indicies_which_intersect[s]}date{REMA_shape.split('_')[2]}"
        
        column_name =f"nid_{REMA_shapes_df.iloc[s].nid}"
        radar_line[column_name] = pd.Series(elevations).replace(-9999.0, np.nan)
        
        #print(f"elevations printed to line for REMA on {REMA_shapes_df.acquisit_1.iloc[s]}")
        
        print(f"{s}/{len(REMA_shapes_df)} of way through REMA strip")
        
        del tiff_stripe_fname, coords, elevations
    
    print(f"{i}/{len(lines_files_paths)} of way through lines")
    
    radar_line.to_file(gis_path+lines_names[i]+".shp")
    
    
    del radar_line
    
# =============================================================================

# write a dictionary which associates each line with REMA strips, 

lines_dict_name = {}
lines_dict_nid = {}

for i, line_file_path in enumerate(lines_files_paths):
    
    radar_line = gpd.read_file(line_file_path)
    
    REMAnid = []
    REMAname = []
        
    for s, REMA_shape in enumerate(REMA_shapes_df.name):
        
        if not REMA_shapes_df.geometry.iloc[s].intersects( LineString(radar_line.geometry.tolist()) ):
            continue
        
        REMAnid.append(f"nid_{REMA_shapes_df.iloc[s].nid}")    
        REMAname.append(REMA_shape)    
    
    print(f"{i}/{len(lines_files_paths)} of way through lines")
    
    lines_dict_nid[lines_names[i]] = REMAnid
    lines_dict_name[lines_names[i]] = REMAname

with open(gis_path+'REMAnid_over_radarlines.txt','w') as f:
    f.write(str(lines_dict_nid))
with open(gis_path+'REMAdate_over_radarlines.txt','w') as f:
    f.write(str(lines_dict_date))
with open(gis_path+'REMAname_over_radarlines.txt','w') as f:
    f.write(str(lines_dict_name))
    
    
# =============================================================================

#PLOT

#open the dictionary associating each line with REMA strips
with open(gis_path+'REMAnid_over_radarlines.txt','r') as f:
    ld = eval(f.read())


line_index = [x for x in zip(lines_names,range(len(lines_names)))]

#print the line names with an index
print(list(zip(range(len(lines_names)),lines_names)))

#with line from first to last point removed
i=14  #which radarline


def plot_line_rm_tide(line_name):
    """
    Plots the radarline with REMA elevations.
    each elevation line has a linear line between first and last points removed, to try and reduce tide effects.
    only use this on lines perpendicular to channel
    """
    
    leg = []
    
    rl = gpd.read_file(gis_path+line_name+".shp")
    
    for REMA in ld[lines_names[i]]:
        #make line from first to last point
        f = interpolate.interp1d( [rl.distan_cum.iloc[0],rl.distan_cum.iloc[-1]], [rl[REMA].iloc[0],rl[REMA].iloc[-1]])
        trendline = f(rl.distan_cum)
    
        plt.plot(rl.distan_cum,rl[REMA]-trendline)
        plt.title(lines_names[i][1:5])
        leg.append(REMA)
    plt.legend(leg) 
    plt.show()


def plot_line(line_name):
      """
    Plots the radarline with REMA elevations.
    no tideline removed.
    """ 
    
    
    #with no trendline removed
    
    rl = gpd.read_file(gis_path+lines_names[i]+".shp")
    
    leg = []
    
    for REMA in ld[lines_names[i]]:
        
        plt.plot(rl.distan_cum,rl[REMA])
        plt.title(lines_names[i][1:5])
        leg.append(REMA)
    plt.legend(leg) 
    plt.show()








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





