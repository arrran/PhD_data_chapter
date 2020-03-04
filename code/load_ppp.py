#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 10:20:44 2020

@author: whitefar
"""

#
#SAMPLE DATA FILE /Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/PPP/full_output_12-11/2019-12-11.pos
#HDR GRP CANADIAN GEODETIC SURVEY, SURVEYOR GENERAL BRANCH, NATURAL RESOURCES CANADA
#HDR ADR GOVERNMENT OF CANADA, 588 BOOTH STREET ROOM 334, OTTAWA ONTARIO K1A 0Y7
#HDR TEL 343-292-6617
#HDR EMA nrcan.geodeticinformation-informationgeodesique.rncan@canada.ca
#NOTE: Estimated positions are at the epoch of data
#DIR FRAME  STN   DAYofYEAR YEAR-MM-DD HR:MN:SS.SS NSV GDOP RMSC(m) RMSP(m)       DLAT(m)       DLON(m)       DHGT(m)          CLK(ns)  TZD(m) SDLAT(95%) SDLON(95%) SDHGT(95%) SDCLK(95%) SDTZD(95%) LATDD LATMN    LATSS LONDD LONMN    LONSS     HGT(m) UTMZONE    UTM_EASTING   UTM_NORTHING UTM_SCLPNT UTM_SCLCBN
#BWD NAD83 2019  345.000000 2019-12-11 00:00:00.00  12  2.0   0.210  0.0042        0.2390        0.3555       -0.5027      144460.0397  2.2740     0.0249     0.0215     0.0876     0.1869     0.0029   -82    40 16.27626  -156     1 11.59742    26.9273      -1         0.0000         0.0000   0.000000   0.000000
#BWD NAD83 2019  345.000012 2019-12-11 00:00:01.00  12  1.5   0.233  0.0040        0.2401        0.3553       -0.5002      142822.0407  2.2740     0.0247     0.0215     0.0864     0.1838     0.0029   -82    40 16.27623  -156     1 11.59746    26.9297      -1         0.0000         0.0000   0.000000   0.000000
#BWD NAD83 2019  345.000023 2019-12-11 00:00:02.00  12  1.5   0.254  0.0040        0.2389        0.3551       -0.4948      141183.9453  2.2740     0.0249     0.0217     0.0869     0.1847     0.0029   -82    40 16.27627  -156     1 11.59751    26.9351      -1         0.0000         0.0000   0.000000   0.000000
#BWD NAD83 2019  345.000035 2019-12-11 00:00:03.00  12  1.5   0.262  0.0036        0.2394        0.3554       -0.4945      139545.8715  2.2740     0.0250     0.0219     0.0874     0.1854     0.0029   -82    40 16.27625  -156     1 11.59744    26.9354      -1         0.0000         0.0000   0.000000   0.000000
#BWD NAD83 2019  345.000046 2019-12-11 00:00:04.00  12  1.5   0.262  0.0035        0.2407        0.3553       -0.4935      137907.7596  2.2740     0.0251     0.0221     0.0879     0.1862     0.0029   -82    40 16.27621  -156     1 11.59747    26.9364      -1         0.0000         0.0000   0.000000   0.000000
#BWD NAD83 2019  345.000058 2019-12-11 00:00:05.00  12  1.5   0.262  0.0036        0.2420        0.3555       -0.4924      136269.4901  2.2740     0.0253     0.0223     0.0883     0.1868     0.0029   -82    40 16.27617  -156     1 11.59741    26.9375      -1         0.0000         0.0000   0.000000   0.000000
#BWD NAD83 2019  345.000069 2019-12-11 00:00:06.00  12  1.5   0.281  0.0036        0.2423        0.3556       -0.4941      134631.1636  2.2740     0.0254     0.0225     0.0887     0.1875     0.0029   -82    40 16.27616  -156     1 11.59738    26.9358      -1         0.0000         0.0000   0.000000   0.000000
#

import numpy as np
import glob
import os
import time
import datetime as dt
import subprocess
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
from geopandas import GeoDataFrame

def load_ppp_date(date):
    """
    INPUT: a date e.g. '2019-12-29'
        """
    
    input_files = f"/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/PPP/**/"+date+".pos"
    
    #find all the input files
    files_paths = glob.glob(input_files)
    
    dfs = []
    
    for file_path in files_paths:
            
        file_path = files_paths[0]
        #date = os.path.splitext(os.path.split(file_path)[1])[0]
    
        dfs.append( pd.read_csv(file_path,header=5,delim_whitespace=True,
                                usecols=['YEAR-MM-DD','HR:MN:SS.SS','LATDD','LATMN','LATSS','LONDD','LONMN','LONSS','HGT(m)'])  )
        
    df = pd.concat(dfs)
    
    #add lat lon in decimal degrees
    df["Latitude"] = df.LATDD.to_numpy() + df.LATMN.to_numpy()/60 + df.LATSS.to_numpy()/3600
    df["Longitude"] = df.LONDD.to_numpy() + df.LONMN.to_numpy()/60 + df.LONSS.to_numpy()/3600
    
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
    gdf = GeoDataFrame(df, crs={'init': 'epsg:3031'}, geometry=geometry,)  
    gdf = gdf.rename(columns={'geometry': 'Points'}).set_geometry('Points').to_crs(epsg=3031)
    
    
    return gdf
