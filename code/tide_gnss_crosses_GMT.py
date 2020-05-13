#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 16:09:22 2020

@author: arran

In this script
"""


# horganhu@arc-gljbrv2:/mnt/hdrive/PROJECTS$ more MELT_GLINE/bin/cross_over.sh
# #!/bin/sh
# # Cross over estimates
# #Create TAG
# x2sys_init TIDES_RADAR -Dgeoz -Ellz -F -Gg
# ls -1 /Users/seismic/PROJECTS/MELT_GLINE/RADAR/RAW/*/MERGED/*exportpicks.txt > radar.lis
# #Convert radar picks to geoz format
# n=1
# max=`wc -l radar.lis | awk '{print $1}'`
# while (test $n -le $max)
# do
# infile=`awk '{if(NR==num) print($0)}' num=$n radar.lis`
# outdir=`dirname $infile`
# outfile=$outdir/`basename $infile _exportpicks.txt`_picks.llz
# awk '{if(substr($1,1,1)!=";") printf("%10.7f %10.7f %5.2f\n", 360+$2, $1, $7)}' $infile > $outfile
# echo converting $infile to $outfile
# n=`expr $n + 1`
# done

# #Convert bedmap picks to geoz format
# infile=/Users/seismic/PROJECTS/MELT_GLINE/BEDMAP/bedmap_subset.txt
# outfile=/Users/seismic/PROJECTS/MELT_GLINE/BEDMAP/bedmap_icethk.llz
# awk -F, '{if((substr($0,1,1)!=";") && ($4 != 0) && ($4!=-9999)) printf("%10.7f %10.7f %5.2f\n", $3, $2, $4)}' $infile > $outfile

# #Create lists
# ls -1 /Users/seismic/PROJECTS/MELT_GLINE/RADAR/RAW/*/MERGED/*_picks.llz > radar.lis
# ls -1 /Users/seismic/PROJECTS/MELT_GLINE/BEDMAP/bedmap_icethk.llz >> radar.lis
# x2sys_cross `cat radar.lis` -TTIDES_RADAR -2 -Qe -Wd0.5k > ../BEDMAP/crossings.txt

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
from shapely.geometry import MultiPoint

gdf = gpd.read_file('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/radar_gnss_track_KIS2.shp')
gdf.crs = "EPSG:3031"
gdf = gdf.round(5)
gdf.drop_duplicates()

out_array = np.hstack((gdf.Longitude.to_numpy().reshape(gdf.Longitude.shape[0],1), gdf.Latitude.to_numpy().reshape(gdf.Latitude.shape[0],1)))
np.savetxt('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/radar_gnss_track_KIS2.xy',out_array)

! export X2

! x2sys_init KIS2 -Dgeo -Exy -F -Gg

! x2sys_cross -TKIS2 radar_gnss_track_KIS2.xy > intersections.xy


! X2SYS_HOME="/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/"

# =============================================================================

# =============================================================================

intersections_raw = pd.read_csv("/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/INTERSECTIONS/intersections.xy",
                                delimiter='\t', skiprows=1,header=0,names=['lon','lat'])

intersections_raw['points'] = [Point(xy) for xy in zip(intersections_raw.lon, intersections_raw.lat)]
intersections_gdf = gpd.GeoDataFrame(intersections_raw,geometry=intersections_raw.points,crs="EPSG:4326")

intersections_points = MultiPoint(intersections_gdf.geometry.tolist())


# .to_crs(3031)
# refine_timesync_func = lambda t : t + pd.Timedelta(Dt)

# self.radata.datetime =  self.radata.datetime.apply(refine_timesync_func)
        
# dist_from_intersections_line = lambda row :print(row.index1)

# gdf.apply(dist_from_intersections_line)

common_points = gdf[gdf.to_crs(4326).geometry.distance(intersections_points) < 0.1]
common_points =common_points.reset_index(drop=True)

common = common_points.groupby(['Latitude', 'Longitude']).geometry
for i,row in common_points.itterows():
    
    
    