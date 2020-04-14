#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 09:54:26 2020

@author: arran
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

from shapely.geometry import LineString

REMA_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/'
#REMA_files_paths = glob.glob(os.path.join(REMA_filepath,"**.tif"),recursive=True)

indicies_which_intersect = np.loadtxt("/Users/home/whitefar/DATA/REMA_2m_strips/indicies_which_intersect.txt").astype(int).tolist()

REMA_shapes_df = gpd.read_file('/Users/home/whitefar/DATA/REMA_2m_strips/REMA_Strip_Index_Rel1.shp').iloc[indicies_which_intersect]


gis_path ="/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RES/PROCESSED_LINES_GISFILE/"

lines_files_paths = glob.glob(os.path.join(gis_path,"**.gpkg"),recursive=True)
lines_names = [os.path.splitext(os.path.split(line_file_path )[1])[0] for line_file_path in lines_files_paths]


# =============================================================================
# 
###
##  1 reprocess line0
# 2 # take a look #REMOVE TRENDLINE and plot


# def sample_tiff():
#     """
#     """
    
    
# =============================================================================

i=13
s=6
# line_file_path = lines_files_paths[i]
# REMA_shape = REMA_shapes_df.name.iloc[s]


#PUT THESE IN FUNCTIONS
#check for interesection
lines_files_paths1 = [lines_files_paths[i]]
REMA_shapes_df1 = REMA_shapes_df

for i, line_file_path in enumerate(lines_files_paths):
    
    radar_line = gpd.read_file(line_file_path).drop(['level_0', 'index'],axis=1)
    
    print("sampling elevations on "+lines_names[i])
        
    for s, REMA_shape in enumerate(REMA_shapes_df.name):
        
        if not REMA_shapes_df.geometry.iloc[s].intersects( LineString(radar_line.geometry.tolist()) ):
            print("no intersection with "+REMA_shape)
            continue
        
        print("yes, intersection with "+REMA_shape)
                
        tiff_stripe_fname = REMA_shape + "_dem.tif"
                
        with rio.open(REMA_filepath + tiff_stripe_fname) as src:
            
            coords = [(x,y) for x, y in zip(radar_line.geometry.x, radar_line.geometry.y)]
            
            elevations = [elevation[0] for elevation in src.sample(coords)]
        
        #column_name =f"i{indicies_which_intersect[s]}date{REMA_shape.split('_')[2]}"
        column_name =f"d{REMA_shape.split('_')[2]}"
        radar_line[column_name] = pd.Series(elevations).replace(-9999.0, np.nan)
        
        print(f"elevations printed to line for REMA on {REMA_shapes_df.acquisit_1.iloc[s]}")
        
        print(f"{s}/{len(REMA_shapes_df)} of way through REMA strip")
        
        del tiff_stripe_fname, coords, elevations
    
    print(f"{i}/{len(lines_files_paths)} of way through lines")
    
    radar_line.to_file(gis_path+lines_names[i]+".shp")
    
    del radar_line
    
# =============================================================================

# =============================================================================

i=13
rl = gpd.read_file(gis_path+lines_names[i]+".shp")

for REMA in ['d20151001', 'd20161109', 'd20141209', 'd20170114','d20170925', 'd20150104', 'd20151010', 'd20121224']:
    plt.plot(rl.distance_a,rl[REMA]-rl[REMA].iloc[-1])
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





