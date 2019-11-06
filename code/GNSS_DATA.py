#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:05:01 2019

@author: huw
"""

import matplotlib.pyplot as pl
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

#datadir = "/Volumes/arc_02/horganhu/GPS_PROCESSING/TAS2016_KINEMATIC/tac2/ARCHIVE_EXPT1/"  # tac2 tac3 both have data


working_directory = "/Users/home/whitefar/DATA/TASMAN/ABSOLUTE/"

os.cwd(working_directory)

units = {
    "tal1":"arc2",
    "tal2":"arc1",
    "tac1":"arc1",
    "tac2":"arc5",
    "tac3":"arc4",
    "tar1":"arc8",
    "tar2":"arc6"
    }

#   converter = {'date_time':lambda s: dt.datetime.strptime(s,'"%Y-%m-%d %H:%M:00"')}
#    log=np.genfromtxt(data_dir + fname,\
#                           dtype=[('date_time', 'O'), ('rec_num', '<i8'),('rain', '<f8'), ('level', '<f8')],
#                           skip_header=header_lines[ind],delimiter=',',\
#                           names=['date_time','rec_num','rain','level'],\
#                           converters=converter)

site = "tac2"

positionfiles = glob.glob("/Volumes/arc_02/horganhu/GPS_PROCESSING/TAS2016_KINEMATIC/"+site+"/ARCHIVE_EXPT1/"+"*.GEOD."+unit+".LC")
path = positionfiles[0]



    
df = pd.read_csv(path, header=0, skiprows=[1], delim_whitespace=True,nrows=10)
geometry = [Point(xy) for xy in zip(df.Latitude, df.Longitude)]
crs = {'init': 'epsg:2263'} #change this
geo_df = GeoDataFrame(df, crs=crs, geometry=geometry)
unit=units.get(site)




def pdloadata(site):
    """
    INPUT: Load txt file from server with GNSS data
    OUTPUT: geodataframe 
    
    Text file looks like this:
        
#*YY  DOY        Seconds        Latitude     Longitude      Height   SigN  SigE  SigH   RMS  #     Atm       +-       Fract DOY     Epoch  #BF NotF
#*                                (deg)       (deg)          (m)     (cm)  (cm)  (cm)  (mm)  DD    (mm)     (mm)
# 2016  250       0.000000  -43.642764669  170.198796821    832.5181   2.1   1.3   3.6   4.5  7  1973.29    43.03   250.00000000000      1  16   0 K
# 2016  250      29.999999  -43.642764704  170.198796877    832.5074   2.1   1.3   3.6   4.6  7  2029.71    36.85   250.00034722222      2  16   0 S
# 2016  250      60.000000  -43.642764573  170.198796785    832.4986   2.1   1.3   3.6   6.6  7  2029.71    36.85   250.00069444445      3  16   0 S
# 2016  250      90.000000  -43.642764776  170.198796723    832.5178   2.1   1.3   3.6   8.3  7  2029.70    36.85   250.00104166666      4  16   0 S

    
    """
    unit=units.get(site)
    print(site,unit)
    positionfiles = glob.glob("/Volumes/arc_02/horganhu/GPS_PROCESSING/TAS2016_KINEMATIC/"
                              +site+"/ARCHIVE_EXPT1/"+"*.GEOD."+unit+".LC")
    
    for ndayfile, positionfile in enumerate(positionfiles):
        print(ndayfile,positionfile)  
        dayfile = np.genfromtxt(positionfile,dtype=datatypes,delimiter=None,
                                skip_header=2,usecols=cols,names=fieldnames)
    
      
    
def loaddata(site):
    
    unit=units.get(site)
    print(site,unit)
    positionfiles = glob.glob("/Volumes/arc_02/horganhu/GPS_PROCESSING/TAS2016_KINEMATIC/"
                              +site+"/ARCHIVE_EXPT1/"+"*.GEOD."+unit+".LC")
    #*YY  DOY        Seconds        Latitude     Longitude      Height   SigN  SigE  SigH   RMS  #     Atm       +-       Fract DOY     Epoch  #BF NotF
    #*                                (deg)       (deg)          (m)     (cm)  (cm)  (cm)  (mm)  DD    (mm)     (mm)
    # 2016  335       0.000000  -43.624278447  170.205003007    911.3857   2.0   1.5   3.9   4.9  6  1986.32    43.91   335.00000000000      1  14   0 K    
    fieldnames=['YY','DOY','SS','LAT','LON','ALT','NE','EE','ALTE','RMS',
                'DDIFF','ATM','ATME','DOYFRAC']
    cols = (0,1,2,3,4,5,6,7,8,9,10,11,12,13)
    datatypes=[('YY',int),('DOY',float),('SS',float),('LAT',float),
               ('LON',float),('ALT',float),('NE',float),('EE',float),
               ('ALTE',float),('RMS',float),('DDIFF',float),('ATM',float),
               ('ATME',float),('DOYFRAC',float)]
    for ndayfile, positionfile in enumerate(positionfiles):
        print(ndayfile,positionfile)  
        dayfile = np.genfromtxt(positionfile,dtype=datatypes,delimiter=None,
                                skip_header=2,usecols=cols,names=fieldnames)
        print("ndayfile:",ndayfile)
        print("dayfilesize:",dayfile.size)
        if ndayfile == 0:
            print("Dayfile 0")
            yy = dayfile['YY']
            doy = dayfile['DOY']
            ss = dayfile['SS']
            lat = dayfile['LAT']
            lon = dayfile['LON']
            alt = dayfile['ALT']
            ne = dayfile['NE']
            ee = dayfile['EE']
            alte = dayfile['ALTE']
            rms = dayfile['RMS']
            ddiff = dayfile['DDIFF']
            atm = dayfile['ATM']
            atme = dayfile['ATME']
            doyfrac = dayfile['DOYFRAC']
        else:
            yy = np.append(yy,dayfile['YY'])
            doy = np.append(doy,dayfile['DOY'])
            ss = np.append(ss,dayfile['SS'])
            lat = np.append(lat,dayfile['LAT'])
            lon = np.append(lon,dayfile['LON'])
            alt = np.append(alt,dayfile['ALT'])
            ne = np.append(ne,dayfile['NE'])
            ee = np.append(ee,dayfile['EE'])
            alte = np.append(alte,dayfile['ALTE'])
            rms = np.append(rms,dayfile['RMS'])
            ddiff = np.append(ddiff,dayfile['DDIFF'])
            atm = np.append(atm,dayfile['ATM'])
            atme = np.append(atme,dayfile['ATME'])
            doyfrac = np.append(doyfrac,dayfile['DOYFRAC'])
    return yy, doy, ss, lat, lon, alt, ne, ee, alte, rms, ddiff, atm, atme, doyfrac


# =============================================================================
# =============================================================================

year, day, seconds, lat, lon, height, sigN, sigE, sigH, rms, ddiff, atm, atme, doyfrac = loaddata("tac2")

import matplotlib.pyplot as plt

distance = np.sqrt(lat**2+lon**2)
velocity = np.append(np.diff(distance),np.diff(distance)[-1])

pos = np.vstack((doyfrac,distance, height, sigN, sigE, sigH, rms, ddiff, atm, atme, )).T
sort_by_colum = 0
posi = pos[pos[:,0].argsort()]

plt.plot(doyfrac,distance,'o')

# =============================================================================
# tm_year=2018, tm_mon=12, tm_mday=27, tm_hour=6, tm_min=35, tm_sec=17, tm_wday=3, tm_yday=361, tm_isdst=0

from_date = (3,20,6,5,30) #(month, day, hour,min,sec)
f_tup = (2016,from_date[0],from_date[1],from_date[2],from_date[3],from_date[4],0,0,0 )
f_secs = time.mktime(f_tup)
f_struc = time.localtime(f_secs)
time.strftime("%m/%d/%Y, %H:%M:%S",f_tup)

from_doy = f_struc[7]
from_index = np.searchsorted(posi[:,0],from_doy)

to_date = (3,20,7,5,30) #(month, day, hour,min,sec)
t_tup = (2016,to_date[0],to_date[1],to_date[2],to_date[3],to_date[4],0,0,0 )
t_secs = time.mktime(t_tup)
t_struc = time.localtime(t_secs)
time.strftime("%m/%d/%Y, %H:%M:%S",t_tup)

to_doy = t_struc[7]
to_index = np.searchsorted(posi[:,0],to_doy)


plt.plot(posi[from_index:to_index,0],posi[from_index:to_index,1],'o')

# =============================================================================
day[2700:2900]
day[5000:6000]
np.sqrt(np.arange(8)**2 +2*np.arange(8)**2)


# Program to filter out only the even items from a list

my_list = [1, 5, 4, 6, 8, 11, 3, 12]

"a": lambda x: x+0.5
a([3,2])



# Output: [4, 6, 8, 12]
print(new_list)