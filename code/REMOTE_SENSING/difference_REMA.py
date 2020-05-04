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

import rasterio as rio
import rasterio.mask

from shapely.geometry import LineString

REMA_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/'
#REMA_files_paths = glob.glob(os.path.join(REMA_filepath,"**.tif"),recursive=True)

# indicies which intersect channel
indicies_which_intersect = np.loadtxt("/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/indicies_which_intersect.txt").astype(int).tolist()

# df = gpd.read_file('/home/arran/PHD/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1/REMA_Strip_Index_Rel1.shp',crs="EPSG:3031")


df = gpd.read_file('/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1.shp')



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




#all REMAs which intersect visible channel
intersects_list = [122083, 122087, 122088, 122089, 131225, 131226, 131228, 145068, 145073, 145074, 150097, 159199, 159200, 159202]


#get dataframe of only the REMA strips over the channel
REMA_shapes_channel = df.iloc[intersects_list]
REMA_shapes_channel = REMA_shapes_channel.assign(stripid=REMA_shapes_channel.index.to_series())
# REMA_shapes_channel.reset_index(drop=True,inplace=True)

#write new columns with intersection polygons
for i,strip in REMA_shapes_channel.iterrows(): 
    REMA_shapes_channel['intersects'+str(strip.stripid)] = REMA_shapes_channel.geometry.intersects(strip.geometry)
    
    print(i)
    
    

        





def difference_rema(df,strip1_index,strip2_index):
    """
    Input: two rema strips, 
    input tiffs is as an index of the shapefile REMA_Strip_Index_Rel1.shp. This has all the info
    
    This script
    1. Finds the intersection polygon using the shapefiles. REMA_Strip_Index_Rel1.shp
    2. cuts both strips to the shapefile and then takes the differnce
    ----
    output difference tiff
    
    """
    
    tiff_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/'
    output_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/DIFFERENCES/'
  
    strip1_ = df.iloc[strip1_index]
    strip2_ = df.iloc[strip2_index]
    
  
    intersection = strip1_.geometry.intersection(strip2_.geometry)
    
    #crop the rema strip .tiff    
    with rio.open(tiff_filepath + strip1_['name'] + "_dem.tif") as src:
        #print(src.transform)
        out_image1, out_transform = rasterio.mask.mask(src, gpd.GeoSeries(intersection),crop=True)
        data_type = out_image1.dtype
        out_meta = src.meta.copy()
    with rio.open(tiff_filepath + strip2_['name'] + "_dem.tif") as src:
        #print(src.transform)
        out_image2, _ = rasterio.mask.mask(src, gpd.GeoSeries(intersection),crop=True)
        
   
    #difference the two tiffs
    diff_image = out_image1 - out_image2
    diff_image[diff_image==0.] = -9999.
    
    out_meta.update({"driver": "GTiff",
                 "height": diff_image.shape[1],
                 "width": diff_image.shape[2],
                 "transform": out_transform,
                 "dtype" : data_type})
    
    print(f"{strip1_index}-{strip2_index}, tiffs cropped")
        
        
        
    #write the cropped tiff to file
    with rio.open(output_filepath + f"REMA_{strip1_index}-{strip2_index}_diff.tif", "w", **out_meta) as dest:
        dest.write(diff_image)
        
    diff_date_str = strip1_.acquisitio[2:4] + strip2_.acquisitio[2:4]
    
        
    print("cropped tiff written to " + output_filepath + f"REMA{diff_date_str}_{strip1_index}-{strip2_index}dif.tif")





# =============================================================================
# MAKE ALL THE DIFFERENCE TIFFS by iterating over the 

#all REMAs which intersect visible channel
intersects_list = [122083, 122087, 122088, 122089, 131225, 131226, 131228, 145068, 145073, 145074, 150097, 159199, 159200, 159202]
df = gpd.read_file('/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1.shp')



# =============================================================================
# which strips intersect which?
REMA_shapes_channel = df.iloc[intersects_list]
REMA_shapes_channel = REMA_shapes_channel.assign(stripid=REMA_shapes_channel.index.to_series())
# REMA_shapes_channel.reset_index(drop=True,inplace=True)

#write new columns with intersection polygons
for i,strip in REMA_shapes_channel.iterrows(): 
    REMA_shapes_channel['intersects'+str(strip.stripid)] = REMA_shapes_channel.geometry.intersects(strip.geometry)
    
    print(i)

#     
#     
# =============================================================================
# Make REMA difference tiffs

for stripid1 in intersects_list[5:]:
    
    print(f"{intersects_list.index(stripid1)} / {len(intersects_list)}")
    
    for stripid2 in intersects_list:
    
        if REMA_shapes_channel[f'intersects{stripid2}'].loc[stripid1] == False:
            continue
        if stripid1 == stripid2:
            continue
        
        difference_rema(df,stripid1,stripid2)
# =============================================================================


    
    
    

