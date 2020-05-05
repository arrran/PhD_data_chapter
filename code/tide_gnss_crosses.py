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

# =============================================================================
# CONCATENATE ALL GNSS DATA INTO ONE GEODATAFRAME

files_paths = glob.glob('/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/PPP/**/*.pos',recursive=True)
dfs = []
for file_path in files_paths: 
    #date = os.path.splitext(os.path.split(file_path)[1])[0]
    dfs.append( pd.read_csv(file_path,header=5,delim_whitespace=True,
                            usecols=['YEAR-MM-DD','HR:MN:SS.SS','LATDD','LATMN','LATSS','LONDD','LONMN','LONSS','HGT(m)'])  )
df = pd.concat(dfs)

#add lat lon in decimal degrees
if df.LATDD.to_numpy()[0] > 0:
    raise ValueError('need to recode, is done for negative latitude')
df["Latitude"] = df.LATDD.to_numpy() - df.LATMN.to_numpy()/60 - df.LATSS.to_numpy()/3600
df["Longitude"] = df.LONDD.to_numpy() - df.LONMN.to_numpy()/60 - df.LONSS.to_numpy()/3600

#remove non decimal lat lon
df = df.drop('LATDD', 1)
df = df.drop('LATMN', 1)
df = df.drop('LATSS', 1)
df = df.drop('LONDD', 1)
df = df.drop('LONMN', 1)
df = df.drop('LONSS', 1)

#add datetime object

df['datetime'] = pd.to_datetime(df['YEAR-MM-DD'] + 'T' + df['HR:MN:SS.SS'])

#remove original date and time
df = df.drop('YEAR-MM-DD', 1)
df = df.drop('HR:MN:SS.SS', 1)

#another column with timestamp
ts_func = lambda t : t.timestamp()
df['timestamp'] = df.datetime.apply(ts_func)

geometry = [Point(xy) for xy in zip(df.Longitude, df.Latitude)]
gdf = GeoDataFrame(df, crs={'init': 'epsg:4326'}, geometry=geometry,)  
gdf = gdf.rename(columns={'geometry': 'Points'}).set_geometry('Points').to_crs(epsg=3031)

# 
# =============================================================================
# CUT TO KIS2, everything x < -370000

gdf = gdf[gdf.geometry.x < -370000]

gdf.reset_index(drop=True)

gdf['index1']=gdf.index
# =============================================================================
# 
def intersects_with(row,gdf):
    
    i=row.index1
    
    if i-6<0:
        f = 0
    else:
        f=i-6
    if i+6 > len(gdf):
        t=len(gdf)
    else:
        t = i+6
        
    intx = gdf.drop(range(f,t))[gdf.drop(range(f,t)).geometry.intersects(row.Points.buffer(10))
                                == True].index.tolist()        
        
    if len(intx) != 0:
        #print('found intersection')
        return intx
    else:
        return False

gdf_all = []

point_intersects_all = []
index_all =[]

for i in tqdm(range(gdf.shape[0])):
    
    point_intersects = intersects_with(gdf.iloc[i],gdf)

    if point_intersects != False:
        point_intersects_all.append(point_intersects)
        index_all.append(gdf.index1.iloc[i])
    
    if (i%5000 == 0) & (i != 0):
        print(f'cleaning double ups at i = {i}')
        df_intersections = pd.DataFrame({'intersections': point_intersects_all,
                                         'index1': index_all})
        df_intersections = df_intersections[df_intersections.index1.diff() != 1 ]
        df_intersections.to_pickle(f'/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/intersects_{i-5000}-{i}.pkl')
        point_intersects_all = []
        index_all =[]  
        del df_intersections
        print(f'i = {i-100} till {i} saved')
    
df_intersections = pd.read_pickle('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/intersects_{i-1000}-{i}.pkl')

# =============================================================================

    


np.save( '/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/intersect_list.npy', gdf['intersects_points'].to_numpy() )




df = pd.DataFrame({'a': np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]),
                   'b': np.array([3, 3, 3, 3, 3, 3, 3, 3, 3])})
df['index1']=df.index

def intersects_with(row):
    print(row.index)
    

df.apply( intersects_with,axis=1)  

ts_func = lambda row :
    
    
    
    
for i in tqdm(range(10000)):
    e=i+i*2
    
# =============================================================================
# =============================================================================
# # 
# =============================================================================
# =============================================================================

# df2['b'] = np.array([False, False, False, True, True, True, False, True, False])
# df2.drop(range(3,6))[df2.drop(range(3,6)).b==True]
intersect_df = pd.DataFrame({'point_2': [i[1] for i in intersect_list],
                       'point_1': [i[0] for i in intersect_list]})

intersect_df.drop_duplicates('point_2',inplace=True)
    
# =============================================================================
#     

#workings to remove double ups
# a = [0,1,2,2,4,13,3,4,10,11,12,13]
# list(set(a)) #removes double ups
# b = list(zip(range(len(a)),a))
# list(dict.fromkeys(a))

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

duplicates_indicies=[]

for item in intersect_list1:
    duplicates_indicies.append(list_duplicates_of(intersect_list1,item)) # fix this [:][1]

# =============================================================================

source = "ABABDBAAEDSBQEWBAFLSAFB"
print(list_duplicates_of(source, 'B'))

for l in a:
    print(l)
    print(list_duplicates_of(a,l))


a = [[i,[1,2,3]] for i in range(5)]


b = [i[1] for i in a]
#workings to remove the points around the point    
