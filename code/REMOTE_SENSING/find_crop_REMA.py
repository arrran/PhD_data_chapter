#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:53:44 2020

@author: arran
"""
#List of indicies in gpd.read_file('REMA_Strip_Index_Rel1.shp') that overlap fieldsite
[122083, 122087, 122088, 122089, 125681, 125682, 131225, 131226, 131228, 131733, 131734, 145068, 145069, 145073, 145074, 150097, 159199, 159200, 159202, 159203, 159536, 159538, 172558, 172560, 172561, 178272]

#from https://www.pgc.umn.edu/data/rema/ scroll down and get REMA Stip Index (Esri Shapefile)


import tarfile

import os
import sys

sys.path.append(os.path.abspath('/Users/home/whitefar/DATA/code/'))


#from download_data import download_to_path copy this so use old python
import rasterio as rio
import rasterio.mask
import fiona
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np

target_area = gpd.read_file("/Users/home/whitefar/DATA/REMA_2m_strips/study_area_buffer_geo.shp")




# =============================================================================
    
def crop_REMA(i, target_area, output_filepath = '/Volumes/arc_04/whitefar/DATA/REMA_STRIPES'):
    
    """
    """
    
    stripe_name = df.name.iloc[i]
    
    temp_directory = '/Users/home/whitefar/DATA/tmp/'
    
    download_to_path(temp_directory + stripe_name + '.tar.gz', df.fileurl.iloc[i])
    print(f'i = {i}, tar.gz downloaded')
        
    zipped_stripe_path =  temp_directory + stripe_name + '.tar.gz'
    
    #path for output tiff
    tiff_stripe_fname = stripe_name + "_dem.tif"
    shape_stripe_fname = "index/"+stripe_name + "_index.shp"
    
    #extract the .tar.gz and save 2 files; .tif and .shp
    with tarfile.open(zipped_stripe_path) as tar:
        stripe_tiff = tar.extract(member=tiff_stripe_fname, path=temp_directory)
        stripe_shape = tar.extract(member=shape_stripe_fname, path= temp_directory)
        stripe_shape = tar.extract(member=shape_stripe_fname[:-4]+'.shx', path= temp_directory)
        stripe_shape = tar.extract(member=shape_stripe_fname[:-4]+'.prj', path= temp_directory)
        stripe_shape = tar.extract(member=shape_stripe_fname[:-4]+'.dbf', path= temp_directory)
    print(f'i = {i}, tar.gz extracted')
    
    #if the tar.gz rema strip intersects the target_area
    if not gpd.read_file(temp_directory + shape_stripe_fname).geometry.intersects(target_area.geometry).iloc[0]:
        print(f'i = {i}, No intersection between target area and tiff {df.name.iloc[i]}')
        
        os.remove(zipped_stripe_path)
        os.remove(temp_directory + tiff_stripe_fname)
        os.remove(temp_directory + shape_stripe_fname)
        os.remove(temp_directory + shape_stripe_fname[:-4]+'.shx')
        os.remove(temp_directory + shape_stripe_fname[:-4]+'.prj')
        os.remove(temp_directory + shape_stripe_fname[:-4]+'.dbf')
        del tiff_stripe_fname, shape_stripe_fname, zipped_stripe_path, temp_directory, stripe_name
        
        #False = no overlap
        return False
    else:
        print(f'i = {i}, Yes, intersection between target area and tiff {df.name.iloc[i]}')
        

    #crop the rema strip .tiff    
    #tiff is in UTM 3031, so convert target_area to_crs 3031 in the process
    with rio.open(temp_directory + tiff_stripe_fname) as src:
        #print(src.transform)
        out_image, out_transform = rasterio.mask.mask(src, target_area.geometry.to_crs(epsg=3031),crop=True)
        data_type = out_image.dtype
        out_meta = src.meta.copy()
    
        
    #print(out_transform)
    
    out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform,
                 "dtype" : data_type})
    
    print(f"i = {i}, tiff cropped")
        
    #write the cropped tiff to file
    with rio.open(output_filepath + tiff_stripe_fname, "w", **out_meta) as dest:
        dest.write(out_image)
        
    print("cropped tiff written to " + output_filepath + tiff_stripe_fname)
    
    #remove the uncropped tiff and tarfile
    os.remove(zipped_stripe_path)
    os.remove(temp_directory + tiff_stripe_fname)
    os.remove(temp_directory + shape_stripe_fname)
    os.remove(temp_directory + shape_stripe_fname[:-4]+'.shx')
    os.remove(temp_directory + shape_stripe_fname[:-4]+'.prj')
    os.remove(temp_directory + shape_stripe_fname[:-4]+'.dbf')
    del (tiff_stripe_fname, stripe_tiff, stripe_shape, shape_stripe_fname, zipped_stripe_path,
         temp_directory, stripe_name, out_image, out_transform, data_type,
         out_meta)
    
    #True = true that it overlapped the area and outputted a cropped tiff
    return True
    
# =============================================================================


# =============================================================================
# Check intersection on home laptop

df = gpd.read_file('/home/arran/PHD/DATA/REMOTE_SENSING/REMA_2m_strips/REMA_Strip_Index_Rel1/REMA_Strip_Index_Rel1.shp')

field_area_df = gpd.read_file('/home/arran/PHD/DATA/REMOTE_SENSING/REMA_2m_strips/study_area_buffer_geo.shp')

field_area = field_area_df.geometry.to_crs(epsg=3031).iloc[0]

intersects_list = []
for i in range(df.shape[0]):
    
    if df.geometry.iloc[i].intersects(field_area):
        intersects_list.append(i)
        print( intersects_list)
    

# =============================================================================



df = pd.read_csv('/Users/home/whitefar/DATA/REMA_2m_strips/KAMB_CHANNEL/attribute_table_stripes_over_channel.txt',delimiter='\t')



intersects_list = []
for i in range(df.shape[0]):
    intersects = crop_REMA(i,target_area)
    if intersects:
        intersects_list.append(i)
        
np.savetxt('/Volumes/arc_04/whitefar/DATA/REMA_STRIPES/indicies_which_intersect.txt',np.array(intersects_list))