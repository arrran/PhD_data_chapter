#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:53:44 2020

@author: arran
"""
# i= 166 "SETSM_WV01_20161109_1020010058134D00_10200100576C9100_seg1_2m_v1.0",
# i= 93 "SETSM_WV01_20161027_1020010056BA2E00_1020010058DFD800_seg1_2m_v1.0",
# i = 23"SETSM_WV02_20161220_1030010061866000_103001006007BA00_seg1_2m_v1.0",
# i= 15 "SETSM_WV03_20161220_1040010026391800_1040010026480500_seg1_2m_v1.0"]

# for goo in good:
#     for i,item in enumerate(df.name.to_list()):
#         if goo==item:
#             print(i)

import tarfile

import os
import sys

sys.path.append(os.path.abspath('/Users/home/whitefar/DATA/code/'))


from download_data import download_to_path
import rasterio as rio
import rasterio.mask
import fiona
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon

target_area = gpd.read_file("/Users/home/whitefar/DATA/REMA_2m_strips/study_area_buffer_geo.shp")

# with fiona.open("/Users/home/whitefar/DATA/REMA_2m_strips/study_area_buffer_geo.shp", "r") as shapefile:
#     field_area = [feature["geometry"] for feature in shapefile]

i= 166    

# CROSSES AINT WERKING!
    
def crop_REMA(i, output_filepath = '/Volumes/arc_04/whitefar/DATA/REMA_STRIPES'):
    
    """
    """
    
    stripe_name = df.name.iloc[i]
    
    temp_directory = '/Volumes/arc_04/whitefar/DATA/REMA_STRIPES/tmp/'
    
    download_to_path(temp_directory + stripe_name + '.tar.gz', df.fileurl.iloc[i])
        
    zipped_stripe_path =  temp_directory + stripe_name + '.tar.gz'
    
    #path for output tiff
    tiff_stripe_fname = stripe_name + "_dem.tif"
    shape_stripe_fname = "index/"+stripe_name + "_index.shp"
    
    #extract the .tar.gz and save 2 files; .tif and .shp
    with tarfile.open(zipped_stripe_path) as tar:
        stripe_tiff = tar.extract(member=tiff_stripe_fname, path=temp_directory)
        stripe_shape = tar.extract(member=shape_stripe_fname, path= temp_directory)
        stripe_shape1 = tar.extract(member=shape_stripe_fname[:-4]+'.shx', path= temp_directory)
        stripe_shape1 = tar.extract(member=shape_stripe_fname[:-4]+'.prj', path= temp_directory)
        stripe_shape1 = tar.extract(member=shape_stripe_fname[:-4]+'.dbf', path= temp_directory)
    
    #if the tar.gz rema strip intersects the target_area
    if not gpd.read_file(temp_directory + shape_stripe_fname).geometry.crosses(target_area.geometry).iloc[0]:
        print(f'No intersection between target area and tiff {df.name.iloc[i]}')
    
        os.remove(tiff_stripe_fname)
        os.remove(shape_stripe_fname)
        os.remove(shape_stripe_fname[:-4]+'.shx')
        os.remove(shape_stripe_fname[:-4]+'.prj')
        os.remove(shape_stripe_fname[:-4]+'.dbf')
        del tiff_stripe_fname, shape_stripe_fname, zipped_stripe_path, temp_directory, stripe_name
        
        return
        
        
    
    
    #crop the .tiff    
    with rio.open(temp_directory + tiff_stripe_fname) as src:
        #print(src.transform)
        out_image, out_transform = rasterio.mask.mask(src, target_area.geometry,crop=True)
        data_type = out_image.dtype
        out_meta = src.meta.copy()
        
    #print(out_transform)
    
    out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform,
                 "dtype" : data_type})
    
    cropped_stripe_path = output_filepath + tiff_stripe_fname
    
    #write the cropped tiff to file
    with rio.open(cropped_stripe_path, "w", **out_meta) as dest:
        dest.write(out_image)
    
    #remove the uncropped tiff
    os.remove(temp_directory + tiff_stripe_fname)
    os.remove(temp_directory + stripe_name + '.tar.gz')
    



df = pd.read_csv('/Users/home/whitefar/DATA/REMA_2m_strips/KAMB_CHANNEL/attribute_table_stripes_over_channel.txt',delimiter='\t')

for i in df.shape[0]:
    crop_REMA(i, output_filepath = '/Volumes/arc_04/whitefar/DATA/REMA_STRIPES')