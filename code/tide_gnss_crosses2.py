#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 16:09:22 2020

@author: arran

Output a shapefile with tracks from GNSS around KIS2 doing radarlines etc
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
# CONCATENATE ALL GNSS DATA recorded on the radar sled INTO ONE GEODATAFRAME (rough copy of Load_ppp)

files_paths = glob.glob('/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/PPP/**/*.pos',recursive=True)
dfs = []
for file_path in files_paths: 
    #date = os.path.splitext(os.path.split(file_path)[1])[0]
    dfs.append( pd.read_csv(file_path,header=5,delim_whitespace=True,
                            usecols=['YEAR-MM-DD','HR:MN:SS.SS','LATDD','LATMN','LATSS','LONDD','LONMN','LONSS','HGT(m)'])  )
df = pd.concat(dfs)
df.reset_index(drop=True,inplace=True)


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

gdf.reset_index(drop=True,inplace=True)

gdf['index1']=gdf.index


gdf = gdf.drop("datetime",1)
gdf.to_file('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/radar_gnss_track_KIS2.shp')
