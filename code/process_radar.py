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
import datetime 
import pandas as pd
import geopandas as gpd

from geopandas import GeoDataFrame
from shapely.geometry import Point
import fiona

#class radarline:
#    
#    def set_linecode(self,linecode):
#        """
#        linecode is the code given to the radar line eg 06348013011
#        """
#        self.linecode = linecode
#        self.ch1_filename = linecode + "ch1"
#        self.ch0_filename = linecode + "ch0"
#        self.info_filename = linecode + "info.txt"
#        self.time_filename = linecode + "time.txt"
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

linecode = "06364035101"
path = "/Volumes/arc_04/FIELD_DATA/K8621920/RES"

              
def set_linecode(linecode):
        """
        linecode is the code given to the radar line eg 06348013011
        """
        linecode = linecode
        ch1_filename = linecode + "ch1"
        ch0_filename = linecode + "ch0"
        info_filename = linecode + "info.txt"
        time_filename = linecode + "time.txt"
        filenames = [ch0_filename,ch1_filename,info_filename, time_filename]
        
def load_data(path):
        """
        load radar data from path and all directories beneath it
        """
        
        files_paths = [glob.glob(os.path.join(path,"**",filename),recursive=True)[0] for filename in  filenames]
        
        print("loading files from:",files_paths)
        
        ch0 = np.fromfile( files_paths[0] ,dtype="float64").reshape(2500,-1) #two binary files
        ch1 = np.fromfile( files_paths[1],dtype="float64" ).reshape(2500,-1)
        info = np.genfromtxt( files_paths[2] ,delimiter=',') # two text files
        time = [dt.datetime(2019,1,1) - dt.timedelta(days=1) + dt.timedelta(t) for t in np.genfromtxt( files_paths[3] )]
        time_str = [t.strftime("%H:%M:%S %d%b%y") for t in time]



        
        

