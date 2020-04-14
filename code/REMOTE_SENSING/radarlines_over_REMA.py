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
## reprocess line0
#run this file for KIS2 lines only


# def sample_tiff():
#     """
#     """
    
    
# =============================================================================

# line_file_path = lines_files_paths[0]
# REMA_shape = REMA_shapes_df.name.iloc[0]
# i=0 

#PUT THESE IN FUNCTIONS
#check for interesection
lines_files_paths1 = [lines_files_paths[0]]
REMA_shapes_df1 = REMA_shapes_df[:4]

for i, line_file_path in enumerate(lines_files_paths):
    
    radar_line = gpd.read_file(line_file_path).drop(['level_0', 'index', 'year', 'day', 'hour',
                                             'minute', 'second', 'x', 'y'],axis=1)
    
    print("sampling elevations on "+lines_names[i])
        
    for s, REMA_shape in enumerate(REMA_shapes_df.name):
        
        if not REMA_shapes_df.geometry.iloc[s].intersects( LineString(radar_line.geometry.tolist()) ):
            #print("no intersection with "+REMA_shape)
            continue
        
        #print("yes, intersection with "+REMA_shape)
                
        tiff_stripe_fname = REMA_shape + "_dem.tif"
                
        with rio.open(REMA_filepath + tiff_stripe_fname) as src:
            
            coords = [(x,y) for x, y in zip(radar_line.geometry.x, radar_line.geometry.y)]
            
            elevations = [elevation[0] for elevation in src.sample(coords) for coords in radar_line.geometry.tolist()]
                
        radar_line[REMA_shape.replace('.','o')] = pd.Series(elevations).replace(-9999.0, np.nan)
        
        #print(f"elevations printed to line for REMA on {REMA_shapes_df.acquisit_1.iloc[s]}")
        
        print(f"{s}/{len(REMA_shapes_df)} of way through REMA strip")
        
        del tiff_stripe_fname, coords, elevations
    
    print(f"{i}/{len(lines_files_paths)} of way through lines")
    
    radar_line.to_file(gis_path+lines_names[i]+".shp")
    
    del radar_line
    
    
    
    









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





