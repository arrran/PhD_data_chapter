#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 09:54:26 2020

@author: arran

This script plots differences in REMA strips
"""

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

REMA_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/'
#REMA_files_paths = glob.glob(os.path.join(REMA_filepath,"**.tif"),recursive=True)

indicies_which_intersect = np.loadtxt("/Users/home/whitefar/DATA/REMA_2m_strips/indicies_which_intersect.txt").astype(int).tolist()

df = gpd.read_file('/home/arran/PHD/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1/REMA_Strip_Index_Rel1.shp',crs="EPSG:3031")


df = gpd.read_file('/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1.shp')

field_area_df = gpd.read_file("/Users/home/whitefar/DATA/REMA_2m_strips/study_area_buffer_geo.shp",crs="EPSG:3031")



# =============================================================================
#on laptop

#find all REMAs which intersect visible channel


channel_df = gpd.read_file("/home/arran/PHD/DATA/REMOTE_SENSING/fieldwork_shapefiles/visible_channel.shp")

channel = channel_df.geometry.to_crs(epsg=3031).iloc[0]

intersects_list = []
for i in range(df.shape[0]):
    
    if df.geometry.iloc[i].intersects(channel):
        intersects_list.append(i)

        print( intersects_list)
    
np.savetxt("/home/arran/PHD/DATA/REMOTE_SENSING/fIeldwork_shapefiles/indicies_which_intersect_channel.txt",np.array(intersects_list))

# =============================================================================


# intersects_list = [122083, 122087, 122088, 122089, 131225, 131226, 131228, 145068, 145073, 145074, 150097, 159199, 159200, 159202]


#find which rema strips intersect which other ones and output in a dictionary

REMA_shapes_channel = df.iloc[intersects_list]

intersects = {}

for i in range(REMA_shapes_channel.shape[0]):
    strip0 = REMA_shapes_channel.geometry.iloc[i]
    striplist = []
    for j in range(REMA_shapes_channel.shape[0]):
        strip1 = REMA_shapes_channel.geometry.iloc[j]
        if strip1.intersects(strip0):
            striplist.append(intersects_list[j])
    print(len(striplist))
    
    intersects[str(intersects_list[i])] = striplist
    del striplist
    
    



def difference_rema(strip1_index,strip2_index):
    """
    Input: two rema strips, tiffs,
    input tiffs is as an index of the shapefile REMA_Strip_Index_Rel1.shp. This has all the info
    
    This script
    1. Finds the intersection polygon using the shapefiles. REMA_Strip_Index_Rel1.shp
    2. cuts both strips to the shapefile and then takes the differnce
    ----
    output difference tiff
    
    """
    
    tiff_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/'
    
    
    df =   gpd.read_file('/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1.shp')         
    
    strip1_ = df.iloc[strip1_index]
    strip2_ = df.iloc[strip2_index]
    
    intersection = strip1_.geometry.intersection(strip2_.geometry)
    
    #crop the rema strip .tiff    
    with rio.open(tiff_filepath + strip1_.name) as src:
        #print(src.transform)
        out_image1, out_transform = rasterio.mask.mask(src, intersection,crop=True)
        data_type = out_image.dtype
        out_meta2 = src.meta.copy()
    with rio.open(temp_directory + strip2_.name) as src:
        #print(src.transform)
        out_image2, out_transform = rasterio.mask.mask(src, intersection,crop=True)
        data_type = out_image.dtype
        out_meta2 = src.meta.copy()   
        
    diff_image = out_image1 - out_image2
    
    #difference the two tiffs
    
    out_meta.update({"driver": "GTiff",
                 "height": diff_image.shape[1],
                 "width": diff_image.shape[2],
                 "transform": out_transform,
                 "dtype" : data_type})
    
    print(f"i = {i}, tiff cropped")
        
        
        
    #write the cropped tiff to file
    with rio.open(output_filepath + tiff_stripe_fname, "w", **out_meta) as dest:
        dest.write(diff_image)
        
    print("cropped tiff written to " + output_filepath + tiff_stripe_fname)









        
    
        
    #print(out_transform)
    
    out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform,
                 "dtype" : data_type})
    
    print(f"i = {i}, tiff cropped")
    
    
    
    
    
    
for i,strip in df.iloc[intersects_list].itterows():
    
    #first 
    
    stripe_name = df.name.iloc[j]
    tiff_stripe_fname = stripe_name + "_dem.tif"
    
    

