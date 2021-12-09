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

# % program   : RTKPOST ver.2.4.3 b33
# % inp file  : \\staff\Home\SCIFAC\whitefa1\DocumentsRedir\GNSS_buck\processed_rinex_rtk\buck83.obs
# % inp file  : \\staff\Home\SCIFAC\whitefa1\DocumentsRedir\GNSS_buck\processed_rinex_files\wgtn_airport\WGTN_165_11to17
# % inp file  : \\staff\Home\SCIFAC\whitefa1\DocumentsRedir\GNSS_buck\processed_rinex_rtk\buck83.nav
# % inp file  : \\staff\Home\SCIFAC\whitefa1\DocumentsRedir\GNSS_buck\ephemeris_gps\brdc1650.20g.Z
# % inp file  : \\staff\Home\SCIFAC\whitefa1\DocumentsRedir\GNSS_buck\ephemeris_gps\brdc1650.20n.Z
# % obs start : 2020/06/13 00:30:57.0 GPST (week2109 520257.0s)
# % obs end   : 2020/06/13 05:13:26.0 GPST (week2109 537206.0s)
# % pos mode  : moving-base
# % freqs     : L1+L2
# % solution  : forward
# % elev mask : 0.0 deg
# % dynamics  : off
# % tidecorr  : off
# % ionos opt : broadcast
# % tropo opt : saastamoinen
# % ephemeris : broadcast
# % navi sys  : gps glonass
# % amb res   : continuous
# % amb glo   : on
# % val thres : 3.0
# % antenna1  :                       ( 0.0000  0.0000  0.0000)
# % antenna2  :                       ( 0.0000  0.0000  0.0000)
# %
# % (lat/lon/height=WGS84/ellipsoidal,Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)
# %  GPST                  latitude(deg) longitude(deg)  height(m)   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio
# 2020/06/13 00:32:38.000  -41.335287074  174.784435390   149.5281   1  11   0.0071   0.0052   0.0119  -0.0024  -0.0021  -0.0017  -0.00    3.0
# 2020/06/13 00:32:39.000  -41.335290494  174.784438816   148.8965   2  11   0.4951   0.3653   0.8312  -0.1660  -0.1472  -0.1129  -0.00    2.9
# 2020/06/13 00:32:40.000  -41.335289028  174.784438162   148.9176   2  11   0.4057   0.2991   0.6805  -0.1361  -0.1203  -0.0940  -0.00    2.6
# 2020/06/13 00:32:41.000  -41.335287200  174.784439318   149.0411   2  11   0.3520   0.2594   0.5901  -0.1181  -0.1042  -0.0819  -0.00    2.8
# 2020/06/13 00:32:42.000  -41.335288416  174.784436740   148.9743   2  11 

def load_ppp(files_paths = ['/Users/home/whitefar/Documents/Stuff/GNSS_buck/outputs/buck_day1.pos']):
    """
    INPUT: a list of file paths
        """
    
    dfs = []
    
    for file_path in files_paths:
            
        #date = os.path.splitext(os.path.split(file_path)[1])[0]
    
        dfs.append( pd.read_csv(file_path,header=25,delim_whitespace=True,) )
                    # usecols=['GPST','latitude(deg)', 'longitude(deg)', 'height(m)'   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio])  
    
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
    
    
    return gdf
