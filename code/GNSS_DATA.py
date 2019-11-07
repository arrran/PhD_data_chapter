#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:05:01 2019

@author: whitefar
"""

import matplotlib.pyplot as pl
import numpy as np
import glob
from numpy import linalg as LA
from functools import reduce
import os
import time
import datetime 
import pandas as pd
import geopandas as gpd

from geopandas import GeoDataFrame
from shapely.geometry import Point
import fiona

#datadir = "/Volumes/arc_02/horganhu/GPS_PROCESSING/TAS2016_KINEMATIC/tac2/ARCHIVE_EXPT1/"  # tac2 tac3 both have data






def load_GNSSdata_for_site(site,start,finish):
    """
    Load txt file from server with GNSS data
    INPUT: 
        - site
        - start, datetime.timestamp()
        - finish, datetime.timestamp()
    
    
    OUTPUT: geodataframe 
    
    Text file looks like this:
        
    #*YY  DOY        Seconds        Latitude     Longitude      Height   SigN  SigE  SigH   RMS  #     Atm       +-       Fract DOY     Epoch  #BF NotF
    #*                                (deg)       (deg)          (m)     (cm)  (cm)  (cm)  (mm)  DD    (mm)     (mm)
    # 2016  250       0.000000  -43.642764669  170.198796821    832.5181   2.1   1.3   3.6   4.5  7  1973.29    43.03   250.00000000000      1  16   0 K
    # 2016  250      29.999999  -43.642764704  170.198796877    832.5074   2.1   1.3   3.6   4.6  7  2029.71    36.85   250.00034722222      2  16   0 S
    # 2016  250      60.000000  -43.642764573  170.198796785    832.4986   2.1   1.3   3.6   6.6  7  2029.71    36.85   250.00069444445      3  16   0 S
    # 2016  250      90.000000  -43.642764776  170.198796723    832.5178   2.1   1.3   3.6   8.3  7  2029.70    36.85   250.00104166666      4  16   0 S

    
    """
    # attribute sites to GNSS units 
    units = { 
    "tal1":"arc2",
    "tal2":"arc1",
    "tac1":"arc1",
    "tac2":"arc5",
    "tac3":"arc4",
    "tar1":"arc8",
    "tar2":"arc6"
    }

    unit=units.get(site)
    print(site,unit)
    #look up all file paths for the GNSS site
    files_paths = glob.glob("/Volumes/arc_02/horganhu/GPS_PROCESSING/TAS2016_KINEMATIC/"  
                            +site+"/ARCHIVE_EXPT1/"+"*.GEOD."+unit+".LC")
    
    if len(files_paths) == 0:
        print(f"No files for {unit} {site}")
        return
    
    #function which converts year day-of-year second to seconds since Jesus
    GNSStime2datetime = lambda gdf : (datetime.datetime(gdf.Year, 1, 1)
                                                 + datetime.timedelta(int(gdf.DOY)-1) +
                                                 datetime.timedelta(seconds=int(gdf.Seconds))).timestamp()
    
    #convertor rounds the seconds to load directly as integer
    converters = {"Seconds": lambda s : int(round(float(s)))} 
        
    dtype = [('*YY',int), ('DOY',int), ("Seconds",float), ('Latitude',float),
                 ('Longitude',float), ('Height',float), ('SigN',float),
                 ('SigE',float), ('SigH',float), ('RMS',float), ('#',float),
                 ('Atm',float), ('+-',float), ('Fract',float), ('DOY.1',float),
                 ('Epoch',int),('#BF',int), ('NotF',str)]
    
    
    
    
     #load each file (one file is made per day)
    for ndayfile, path in enumerate(files_paths): # ndayfile = 0 path = files_paths[0]
        print(ndayfile,"/",len(files_paths),path)  
        
        # First read the timestamp
                
        time_df = pd.read_csv(path, header=0, skiprows=[1], delim_whitespace=True,
                              dtype=dtype, usecols=["*YY","DOY","Seconds"], converters=converters)
        
        time_df.rename(columns={'*YY': 'Year'}, inplace=True)
        time_df["Timestamp"] = time_df.apply(GNSStime2datetime,axis=1)  #convert the time to a new timestamp column
        
        time_df.sort_values(by=['Timestamp'], inplace=True) #sort by the new column
        time_df.reset_index(drop=True, inplace=True)
        
        
        
        #see if there are there any data in the time period
        if not any(np.argwhere( np.logical_and( (time_df.Timestamp.to_numpy() >=start ),( finish>=time_df.Timestamp.to_numpy() ))).flatten()+2):
            continue
            
        #get the indicies for the time period
        time_indicies_skip = np.argwhere( np.logical_or( (time_df.Timestamp.to_numpy() <start ),( finish<time_df.Timestamp.to_numpy() ))).flatten()+2
        
        #now we load the whole dataset with just the time period specified            
        df = pd.read_csv(path, header=0, skiprows=np.concatenate(([1],time_indicies_skip)), delim_whitespace=True, 
                         dtype=dtype, converters=converters)   
        
        df.rename(columns={'*YY': 'Year'}, inplace=True)
        
        df["Timestamp"] = df.apply(GNSStime2datetime,axis=1) 
        
        df["site"] = site
        df["unit"] = unit
        
        if 'df_all' not in locals():
            df_all = df
            print('making new dataframe')
        else:
            df_all = df_all.append(df)
            print('appending to dataframe')
        del df, time_df    
#        
    
    df_all.sort_values(by=['Timestamp', 'site'], inplace=True)
    df_all.reset_index(drop=True, inplace=True)
    
    
    df_all = df_all.reindex(columns=['Timestamp','site', 'unit','Year', 'DOY', 'Seconds', 'Latitude', 'Longitude', 'Height', 'SigN',
       'SigE', 'SigH', 'RMS', 'DDiff', 'Atm', '+-', 'Fract', 'DOY.1', 'Epoch',
       '#BF', 'NotF', 'geometry' ])
    
    return df_all

# =============================================================================

#    
def GNSStime2datetime(year,doy,seconds):
    return (datetime.datetime(year, 1, 1) + datetime.timedelta(int(doy)-1) +\
            datetime.timedelta(seconds=int(seconds))).timestamp()

start = GNSStime2datetime(2016,2,0)
finish = GNSStime2datetime(2016,9,0)

site = "tac2"

gdf = load_GNSSdata_for_site(site,start,finish)
print("still going")

output_folder = "/Volumes/arc_02/whitefar/DATA/TASMAN/GNSS_ABSOLUTE/geodataframe_allGNSS"


units = {
"tal1":"arc2",
"tal2":"arc1",
"tac1":"arc1",
"tac2":"arc5",
"tac3":"arc4",
"tar1":"arc8",
"tar2":"arc6"
}

for site in units:
    site_gdf = load_GNSSdata_for_site(site,start,finish)
    
    if 'geo_df' not in locals():
        geo_df = site_gdf
    else:
        geo_df = geo_df.append(site_gdf, ignore_index=True)
        
    del site_gdf

geo_df.sort_values(by=['Timestamp'], inplace=True)
geo_df.reset_index(drop=True, inplace=True)


geo_df.to_file(output_folder)

del geo_df

# =============================================================================
 #
df = geo_df
geometry = [Point(xy) for xy in zip(df.Latitude, df.Longitude)]
crs = {'init': 'epsg:4326'} 

if 'gdf' not in locals():
    gdf = GeoDataFrame(df, crs=crs, geometry=geometry)
    print('making new dataframe')
else:
    gdf = gdf.append(GeoDataFrame(df, crs=crs, geometry=geometry), ignore_index=True)
    print('appending to dataframe')
    del df, geometry    




#
## =============================================================================
#def test():
#    a =22
#    b = 33
#    print (locals())
#
#
#def load_GNSSdata_time_period()


# =============================================================================
# =============================================================================
# 
# 
# #   converter = {'date_time':lambda s: dt.datetime.strptime(s,'"%Y-%m-%d %H:%M:00"')}
# #   log=np.genfromtxt(data_dir + fname,\
# #                           dtype=[('date_time', 'O'), ('rec_num', '<i8'),('rain', '<f8'), ('level', '<f8')],
# #                           skip_header=header_lines[ind],delimiter=',',\
# #                           names=['date_time','rec_num','rain','level'],\
# #                           converters=converter)
# 
# 
# 
# positionfiles = glob.glob("/Volumes/arc_02/horganhu/GPS_PROCESSING/TAS2016_KINEMATIC/"+site+"/ARCHIVE_EXPT1/"+"*.GEOD."+unit+".LC")
# path = positionfiles[0]
# 
# 
# 
#     
# df = pd.read_csv(path, header=0, skiprows=[1], delim_whitespace=True,nrows=10)
# geometry = [Point(xy) for xy in zip(df.Latitude, df.Longitude)]
# crs = {'init': 'epsg:2263'} #change this
# geo_df = GeoDataFrame(df, crs=crs, geometry=geometry)
# unit=units.get(site)
# 
# =============================================================================
