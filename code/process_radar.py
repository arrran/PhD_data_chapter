#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:22:58 2020

@author: whitefar
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
from numpy import linalg as LA
from functools import reduce
import os
import time
import datetime as dt
import pandas as pd
import geopandas as gpd

from geopandas import GeoDataFrame
from shapely.geometry import Point
import fiona

# =============================================================================
# #check all radar files are in folder
# 
 folder = "/Volumes/arc_04/FIELD_DATA/K8621920/RES"
 
 #do this later
 #converters = {"date_nzst": lambda t : int(round(float(s)))\
 #              "time_nzst": lambda t : } 
 
 dtype = [('line_name',"U50"), ("date_nzst","U12"), ("started_file_nzdt","U12"), ("stopped_file_nzdt","U12"),\
          ("waveforms","int32"), ("filecode","U12")]
 
 meta = np.genfromtxt("/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata",skip_header=1,dtype=dtype)
 
 for filecode in meta['filecode'][:20]:
     if len(glob.glob(os.path.join(folder,"**",str(filecode)+"*"),recursive=True)) ==0:
         print("files are missing for {}".format(filecode))
     else:
         print("file exists for {}".format(filecode))
 
# 
# =============================================================================
#class radarline:
#    
#    def set_filecode(self,filecode):
#        """
#        filecode is the code given to the radar line eg 06348013011
#        """
#        self.filecode = filecode
#        self.ch1_filename = filecode + "ch1"
#        self.ch0_filename = filecode + "ch0"
#        self.info_filename = filecode + "info.txt"
#        self.time_filename = filecode + "time.txt"
#        self.filenames = [self.ch0_filename,self.ch1_filename,self.info_filename, self.time_filename]
#        
#    def load_data(self,path):
#        """
#        load radar data from path and all directories beneath it
#        """
#        

#        self.path = path
#        files_paths = [glob.glob(os.path.join(path,"**",filename),recursive=True)[0] for filename in  self.filenames]
#        
#        print("loading files from:",files_paths)
#        
#        self.ch0 = np.fromfile( files_paths[0] ).reshape(2500,-1) #two binary files
#        self.ch1 = np.fromfile( files_paths[1] ).reshape(2500,-1)
#        self.info = np.genfromtxt( files_paths[2] ) # two text files
#        self.time = np.genfromtxt( files_paths[3] )
#        
# 
#def set_metadata()
#        
#    metadata_path = glob.glob(os.path.join(path,"**","radar_metadata")[0]       

filecode = "06364035101"
#filecode = "06364020457"


path = "/Volumes/arc_04/FIELD_DATA/K8621920/RES"

              
def set_filecode(filecode):
        """
        filecode is the code given to the radar line eg 06348013011
        """
        filecode = filecode
        ch1_filename = filecode + "ch1"
        ch0_filename = filecode + "ch0"
        info_filename = filecode + "info.txt"
        time_filename = filecode + "time.txt"
        filenames = [ch0_filename,ch1_filename,info_filename, time_filename]
        
def load_radar_data(path):
        """
        load radar data from path and all directories beneath it
        
        ...could recode this to import with pandas...
        """
        
        files_paths = [glob.glob(os.path.join(path,"**",filename),recursive=True)[0] for filename in  filenames]
        
        print("loading files from:",files_paths)
        

        
        #two binary files in big endian
        ch0 = np.fromfile( files_paths[0],dtype=">f8", count=-1).reshape(-1,2500)
        ch1 = np.fromfile( files_paths[1],dtype=">f8", count=-1).reshape(-1,2500)
        info = np.genfromtxt( files_paths[2] ,delimiter=',') # two text files
        datetime = [dt.datetime(2019,1,1) - dt.timedelta(days=1) + dt.timedelta(t) for t in np.genfromtxt( files_paths[3] )]
        time_str = np.array([t.strftime("%H:%M:%S %d%b%y") for t in time])
        
        radar = pd.DataFrame({'time':time})
        radar['ch0'] = list(ch0)
        radar['ch1'] = list(ch1)


with open(Filename_ch0, 'rb') as fid:
      data=(fromfile(fid,float64).reshape((-1,int(Samples))).T).byteswap()
def load_gps_data():
        """
        """
        gps_path = "/Users/home/whitefar/DATA/ANT_DATA_1920/RES_GPS/2019-12-30 181325.gpx"
        track_points = gpd.read_file(gps_path,layer='track_points')
        track_points.keys()
        
        track_points['datetime'] = np.array([dt.datetime.strptime(t,"%Y-%m-%dT%H:%M:%S") for t in track_points.time])
        
        radar_to_gps_index = [np.argwhere(abs(track_points.datetime - t)==abs(track_points.datetime - t).min())[0][0] for t in time]
        
        radar['geometry'] = track_points.geometry[radar_to_gps_index].array
        radar['geometry_datetime'] = track_points.datetime[radar_to_gps_index].array
        






    time[140]
    track_points.datetime[radar_to_gps_index[140]]


    date = time[0]
    dates = np.track_points.datetime)
    np.argwhere(abs(dates - date)==abs(dates - date).min())[0][0]




    diff = dates - date
    diff[13666]

    
    gps_data = np.hstack(track_points.geometry.x,track_points.geometry.y,
    
    fiona.listlayers(fname)
    track_points = fiona.open(fname, layer='track_points')
    
    

    layer.crs, layer.bounds
    len(list(track_points.items()))
    track_points


        
        

