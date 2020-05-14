#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 16:09:22 2020

@author: arran
"""
import matplotlib.pyplot as plt
import numpy as np
import glob
from numpy import linalg as LA
from functools import reduce
import os
import sys
import time
import datetime as dt
import pandas as pd
import geopandas as gpd
import scipy as sp
from scipy import signal

from geopandas import GeoDataFrame
from shapely.geometry import Point
import fiona
from tqdm import tqdm


gdf = gpd.read_file('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/radar_gnss_track_KIS2.shp')
gdf.crs = "EPSG:3031"

# temporary distance from previous
tmp_dfp = [Point.distance(gdf.geometry.iloc[i]) for i,Point in enumerate(gdf.geometry.iloc[1:])]
tmp_dfp[:0] = [0]

gdf['dx'] = pd.Series(tmp_dfp)
gdf['distan_cum'] = gdf.dx.cumsum()

to_keep = []
d=10

#Resample by dropping points less than 10m apart
for i,row in tqdm(gdf.iterrows()):
    
    d += row.dx
    
    if d < 1:
        to_keep.append(False)
    else:
        to_keep.append(True)
        d = 0
    
rs = gdf[to_keep]
rs.reset_index(drop=True,inplace=True)

# =============================================================================
# =============================================================================

# =============================================================================
# =============================================================================

buffer_size = 1

#find all the points with other points close by, add them to a list
point_intersects_all = []
index_all =[]

to_skip = []

for i in tqdm(range(0,rs.shape[0])):
    
    #if the inex has already been counted, then skip.
    if i in to_skip:
        continue
    
    # get points which are close in space, and far in time
    intx = rs[ (rs.geometry.intersects(rs.iloc[i].geometry.buffer(buffer_size)) == True) #the 8 is 8m away
              &  (  (rs.timestamp-rs.timestamp.iloc[i])>3600 )
              ].index1.tolist()     
    
    # add these points to a list
    if len(intx) > 0:
        point_intersects_all.append(intx)
        index_all.append(rs.index1.iloc[i])   
    
    to_skip += intx
    
    

# =============================================================================
# =============================================================================

#organise so that each cluster of points is represented by 1 row = 1 point, and the separate height chnges are added

repeated = gdf.iloc[index_all] #points which are repeated
repeated['intx_points'] = point_intersects_all


time_ele_list_all = []
for i,row in tqdm(repeated.iterrows()):
    
    time_ele_list = [] # a list of tuples of time and elevations for each point
    
    time_ele_list.append( [gdf.iloc[i].timestamp,
                           gdf.iloc[i]['HGT(m)'],
                           gdf.iloc[i]['HGT(m)']  - gdf.iloc[i]['HGT(m)'] ])
    
    for index in row.intx_points:
        time_ele_list.append(  [gdf.iloc[index].timestamp,
                                gdf.iloc[index]['HGT(m)'],
                                gdf.iloc[index]['HGT(m)']  - gdf.iloc[i]['HGT(m)'] ])
                                
    
    time_ele_list_all.append(np.array(time_ele_list))
    
                                
repeated['time_ele'] = time_ele_list_all

max_dh_f = lambda row : row[:,2].max()
repeated['max_dh'] =  repeated.time_ele.apply(max_dh_f) #max change in height 
min_dh_f = lambda row : row[:,2].min()
repeated['min_dh'] =  repeated.time_ele.apply(min_dh_f)
repeated['tide'] =  repeated.max_dh - repeated.min_dh


repeated.drop(['time_ele','intx_points'],1).to_file(f'/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/tides_1msample_{buffer_size}mbuffer.shp')
repeated.to_pickle(f'/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/tides_1msample_{buffer_size}mbuffer.pkl')

print(f"done for buffer {buffer_size}m")

# # =============================================================================
# =============================================================================

# Try and work out which points with no tidal signal are like that due to the
# fact the timing was at same time in tidal cycle

repeated = pd.read_pickle('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/tides_10msample_50mbuffer.pkl')

time_diff = []
for i, row in repeated.iterrows():
    #find time difference from min to max
    time_diff.append( row.time_ele[np.argmax(row.time_ele[:,2]),0] - row.time_ele[np.argmin(row.time_ele[:,2]),0])
repeated['time_diff'] = (np.array(time_diff)/44712)%1

repeated.to_pickle('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/tides_10msample_50mbuffer_timediff.pkl')





