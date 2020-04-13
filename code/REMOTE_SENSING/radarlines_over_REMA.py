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

from shapely import LineString

REMA_filepath = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/REMA_STRIPES/'

gis_path ="/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RES/PROCESSED_LINES_GISFILE/"

lines_files_paths = glob.glob(os.path.join(gis_path,"**.gpkg"),recursive=True)
lines_names = [os.path.splitext(os.path.split(line_file_path )[1])[0] for line_file_path in lines_files_paths]












#THIS MAKES A GEODATAFRAME WITH A ROW FOR EACH RADARLINE
lines = []
names = []

for i, line_file_path in enumerate(lines_files_paths):
    df = gpd.read_file(line_file_path)
    lines.append( LineString(df.geometry.tolist()) )
    names.append(lines_names[i])
    del df
    print(i / len(lines_files_paths))


lines_gdf = gpd.GeoDataFrame({'shortname': names}, geometry=lines,crs="EPSG:3031")













LineString(df.geometry.tolist())

for line_file_path in lines_files_paths[1:]:
    df.append ( gpd.read_file(line_file_path) )

for line in df:
        
    with rio.open(temp_directory + tiff_stripe_fname) as src:
        src.sample(line.geometry)
            #print(src.transform)
    


