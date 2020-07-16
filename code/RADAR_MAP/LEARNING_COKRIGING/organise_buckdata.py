#!/usr/bin/env python3
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





!cd /Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/
!gdal_translate -of XYZ DEM_Buck_small.tif DEM_Buck_small.xyz

buck_DEM = pd.read_csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/DEM_Buck_small.xyz',
                       header=None,names = ['x','y','h_DEM'],sep=' ')
buck_DEM['ID'] = 'DEM'

buck_DEM.to_csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_DEM.csv',index=False)
