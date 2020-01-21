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
# folder = "/Volumes/arc_04/FIELD_DATA/K8621920/RES"
# 
# #do this later
# #converters = {"date_nzst": lambda t : int(round(float(s)))\
# #              "time_nzst": lambda t : } 
# 
# dtype = [('line_name',"U50"), ("date_nzst","U12"), ("started_file_nzdt","U12"), ("stopped_file_nzdt","U12"),\
#          ("waveforms","int32"), ("filecode","U12")]
# 
# meta = np.genfromtxt("/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata",skip_header=1,dtype=dtype)
# 
# for filecode in meta['filecode'][:20]:
#     if len(glob.glob(os.path.join(folder,"**",str(filecode)+"*"),recursive=True)) ==0:
#         print("files are missing for {}".format(filecode))
#     else:
#         print("file exists for {}".format(filecode))
 
# 
# =============================================================================

# 
def metadata_func(fc):
    """
    This function reads the metadata for all radar lines.
    
    input: a filecode (fc) 
    output: a callable dataframe with entries [line_name date_nzdt started_file_nzdt stopped_file_nzdt waveforms filecode]
    
    """
        
    metadata_path = "/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata"
   
    #not sure the convertors work, it seems to use pandas time format, which seem compatible with datetime
    converters1 = {"date_nzdt" :lambda t : dt.datetime.strptime(t,"%Y-%m-%d"),\
                "started_file_nzdt": lambda t : dt.datetime.strptime(t,"%H:%M"),\
                "stopped_file_nzdt" :lambda t : dt.datetime.strptime(t,"%H:%M")}
   
   
    metadata = pd.read_csv(metadata_path,delimiter=' ',converters=converters1)
   
    #date was in a separate column to time, add the date to the time series to make datetimes
    metadata.started_file_nzdt = [dt.datetime.combine(metadata.date_nzdt[i].date(),metadata.started_file_nzdt[i].time()) for i,_ in enumerate(metadata.started_file_nzdt)]
    metadata.stopped_file_nzdt = [dt.datetime.combine(metadata.date_nzdt[i].date(),metadata.stopped_file_nzdt[i].time()) for i,_ in enumerate(metadata.stopped_file_nzdt)]
   
    #list of filecodes
    filecodes = np.loadtxt("/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata",dtype=str,skiprows=1,usecols=5)
   
    #dictionary which returns the row in metadata given filecode
    filecode2metarow = {filecode:row for filecode,row in zip(filecodes,np.arange(0,len(filecodes)))}
   
    return metadata.loc[filecode2metarow[fc]]

def set_timesync(date_in):
    """
    This reads the timesync data
    
    input: a date string, format "%Y-%m-%d"
    output: a time delta, the time difference between pixie_time and utc_time written on the file "time_sync"
    """
    converters = {"exact_nzdt" :lambda t : dt.datetime.strptime(str(t),"%d-%m-%YT%H:%M:%S"),\
                "pixie_time": lambda t : dt.datetime.strptime(str(t),"%d-%m-%YT%H:%M:%S"),\
                "exact_utc_time" :lambda t : dt.datetime.strptime(str(t),"%d-%m-%YT%H:%M:%S")}
    
        
    timesync = pd.read_csv("/Volumes/arc_04/FIELD_DATA/K8621920/RES/time_sync",delimiter=' ',converters=converters)
    #get a time delta
    timesync["dt"] = timesync.pixie_time.array - timesync.exact_utc_time.array
    #make a dictionary which returns the time delta given a date
    timesync_dict = {date:dt for date,dt in zip([D.date().strftime("%Y-%m-%d") for D in timesync.exact_nzdt], timesync.dt)}
        
    return timesync_dict[date_in]
    
        
        
        

class radarline:
    
    def __init__(self,filecode):
            """
            filecode is the code given to the radar line eg 06348013011
            """
            self.filecode = filecode
            self.ch1_filename = filecode + "ch1"
            self.ch0_filename = filecode + "ch0"
            self.info_filename = filecode + "info.txt"
            self.time_filename = filecode + "time.txt"
            self.filenames = [self.ch0_filename,self.ch1_filename,self.info_filename, self.time_filename]
            
            
            self.metadata = metadata_func(filecode)
    
    
                  
    def set_filecode(self,filecode):
            """
            filecode is the code given to the radar line eg 06348013011
            """
            self.filecode = filecode
            self.ch1_filename = filecode + "ch1"
            self.ch0_filename = filecode + "ch0"
            self.info_filename = filecode + "info.txt"
            self.time_filename = filecode + "time.txt"
            self.filenames = [self.ch0_filename,self.ch1_filename,self.info_filename, self.time_filename]

            self.metadata = metadata_func(filecode)
            
    def load_radar_data(self,path = "/Volumes/arc_04/FIELD_DATA/K8621920/RES/"):
            """
            load radar data from path and all directories beneath it
            
            ...could recode this to import with pandas...
            """
            
            if len(glob.glob(os.path.join(path,"**",self.filenames[0]),recursive=True)) == 0:
                raise ValueError("Can't find files with that filecode")
            
            self.files_paths = [glob.glob(os.path.join(path,"**",filename),recursive=True)[0] for filename in  self.filenames]
            
            print("loading files from:")
            for file_path in self.files_paths:
                print(file_path)
            
    
            
            #two binary files in big endian
            #ch0 = np.fromfile( files_paths[0],dtype=">f8", count=-1).reshape(-1,2500)
            #ch1 = np.fromfile( files_paths[1],dtype=">f8", count=-1).reshape(-1,2500)
            self.info = np.genfromtxt( self.files_paths[2] ,delimiter=',') # two text files
            self.pixietimes = [dt.datetime(2019,1,1) - dt.timedelta(days=1) + dt.timedelta(t) for t in np.genfromtxt( self.files_paths[3] )]
            
            self.timesync = set_timesync(self.metadata.date_nzdt.strftime("%Y-%m-%d"))
            
            self.datetime = [pixietime + self.timesync for pixietime in self.pixietimes]
            self.time_str = np.array([t.strftime("%H:%M:%S %d%b%y") for t in self.datetime])
            self.pixie_time_str = np.array([t.strftime("%H:%M:%S %d%b%y") for t in self.pixietimes])
            
            self.radata = pd.DataFrame({'time':self.datetime})
            self.radata['ch0'] = list( np.fromfile( self.files_paths[0],dtype=">f8", count=-1).reshape(-1,2500) )
            self.radata['ch1'] = list( np.fromfile( self.files_paths[1],dtype=">f8", count=-1).reshape(-1,2500) )
    
    
    def load_gps_data(self,gps_path = "/Users/home/whitefar/DATA/ANT_DATA_1920/RES_GPS/2019-12-30 181325.gpx"):
            """
            """
            self.track_points = gpd.read_file(gps_path,layer='track_points')
                        
            self.track_points['datetime'] = np.array([dt.datetime.strptime(t,"%Y-%m-%dT%H:%M:%S") for t in self.track_points.time])
            
            self.radar_to_gps_index = [np.argwhere(abs(self.track_points.datetime - t)==abs(self.track_points.datetime - t).min())[0][0] for t in self.datetime]
            
            self.radata['geometry'] = self.track_points.geometry[self.radar_to_gps_index].array
            self.radata['geometry_datetime'] = self.track_points.datetime[self.radar_to_gps_index].array


#filecode = "06364035101"
#filecode = "06364020457"
#filecode = "06001000411"
#filecode = '06001001502'
            

#1-1-2020
linescamp = radarline("06001001502")
#linescamp.set_filecode("06001001502")
#linescamp.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#linescamp.load_gps_data()

##31-1-2020
#lat_apres = radarline()
#lat_apres.set_filecode("06001000411")
#lat_apres.load_radar_data()
#lat_apres.load_gps_data()
#
#up_chan = radarline()
#up_chan.set_filecode("06001000235")
#up_chan.load_radar_data()
#up_chan.load_gps_data()
#
###30-12-2019
#line5 = radarline()
#line5.set_filecode("06364035101")
#line5.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#line5.load_gps_data()
##
##29-12-2019
#line14 = radarline()
#line14.set_filecode("06363041031")
#line14.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#line14.load_gps_data()
#
##28-12-2019
#line11 = radarline()
#line11.set_filecode("06362023503")
#line11.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#line11.load_gps_data()