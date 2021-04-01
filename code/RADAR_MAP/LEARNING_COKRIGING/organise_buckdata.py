,#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 10:48:41 2020

@author: whitefar
"""

import glob
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point, LineString
import geopandas as gpd
from shapely.ops import nearest_points
import numpy as np
from tqdm import tqdm
from scipy import interpolate
import datetime as dt
from scipy.signal import savgol_filter
from shapely.affinity import scale

buck_allpoints = gpd.read_file('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_allpoints.shp').drop(['Latitude (', 'Longitude'],axis=1).rename(columns={'Height (Gl':'h_gps'})

buck_allpoints['x'] = buck_allpoints.geometry.x.copy()
buck_allpoints['y'] = buck_allpoints.geometry.y.copy()

buck_allpoints = buck_allpoints.drop(['geometry'],axis=1).reindex(columns=['x','y','h_gps','ID']).copy()

buck_allpoints.to_csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_allpoints.csv',index=False)

xmin = 1749306.0530032176
xmax = 1749358.4058849094
ymin = 5422466.273025215
ymax = 5422487.752541591

!cd /Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/
!gdal_translate -projwin 1749305 5422489 1749360 5422465 -of GTiff DEM_Buck_small.tif DEM_Buck_overpoints.tif #overpoints is the DEM restricted to a box over the points
NB is  xmin ymax xmax ymin
!cd /Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/
!gdal_translate -of XYZ DEM_Buck_overpoints.tif DEM_Buck_overpoints.xyz

buck_DEM = pd.read_csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/DEM_Buck_small.xyz',
                       header=None,names = ['x','y','h_DEM'],sep=' ')
buck_DEM['ID'] = 'DEM'

buck_DEM.to_csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_DEM.csv',index=False)


buck_DEM = pd.read_csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/DEM_Buck_overpoints.xyz',
                       header=None,names = ['x','y','h_DEM'],sep=' ')
buck_DEM['ID'] = 'DEM'

buck_DEM.to_csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_DEM_overpoints.csv',index=False)

,
#sample points from RASTER
# Check for intersection of line and REMA, then if they intersect, write the elevations of the REMA to the radarline


buck_allpoints = gpd.read_file('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_allpoints.shp').drop(['Latitude (', 'Longitude'],axis=1).rename(columns={'Height (Gl':'h_gps'})

#write REMA elevations to radarline                
with rio.open('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/DEM_Buck_small.tif') as src:
    
    coords = [(x,y) for x, y in zip(buck_allpoints.geometry.x, buck_allpoints.geometry.y)]
    
    elevations = [elevation[0] for elevation in src.sample(coords)]
    
    #column_name =f"i{indicies_which_intersect[s]}date{REMA_shape.split('_')[2]}"
    
    column_name =f"DEM_sampled_points"
    buck_allpoints[column_name] = pd.Series(elevations).replace(-9999.0, np.nan)
    
    #print(f"elevations printed to line for REMA on {REMA_shapes_df.acquisit_1.iloc[s]}")

buck_allpoints['x'] = buck_allpoints.geometry.x.copy()
buck_allpoints['y'] = buck_allpoints.geometry.y.copy()

buck_allpoints = buck_allpoints.drop(['geometry'],axis=1).reindex(columns=['x','y','h_gps',"DEM_sampled_points",'ID']).copy()
    
buck_allpoints.to_csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_allpoints_DEMsampled.csv',index=False)


