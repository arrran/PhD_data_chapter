#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 11:54:01 2020

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
import scipy as sp
from scipy import signal
import subprocess
import gnsscal


files_paths = glob.glob(f"/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/ARC09/192.168.3.1/DSK1/SSN/LOG4_KIS1Hz/**/KIS*")

#del files_paths[100] #this has a zero entry weird

filenames = [os.path.basename(file_path) for file_path in files_paths]

dates = []
datetime_strings = []

for i,filename in enumerate(filenames):
    
    
    
    if filename[4:7] = '0000':
        print(filename + " timestamp is zero, will not process")
        del files_paths[i]
        del filenames[i]
    
   daymonth = dt.date.fromordinal(int(filename[4:7]))
   year = int(dt.datetime.strptime(filename[-3:-1], '%y').strftime('%Y'))   #reads the date eg '19', outputs an integer 2019
   date = daymonth.replace(year=year)
   
   
    if filename[-3:-1] == '19':
        year = 2019
        daymonth = daymonth.replace(2019)
        dates.append(daymonth)
    elif filename[-3:-1] == '20':
        year = 2020
        daymonth = daymonth.replace(2020)
        dates.append(  daymonth)
    else:
        
    
    time = dt.time(ord(filename[7])-97)
    datetime = dt.datetime.combine(date,time)
    datetime_string = datetime.isoformat()[:-6] # a string that looks like this 2020-12-10T21
    datetime_strings.append(datetime_string)
    
out_filenames = ["/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/"+datetime_string for datetime_string in datetime_strings]
 
gpsweeks = [gnsscal.date2gpswd(date)[0] for date in dates] #gpsweek is a number since some epoch, needed for converting to rinex

gpsinterval=30 # GPS sampling interval
antenna_height = 1.55

commands = []
for i,file_path in enumerate(files_paths): 
    command1 = (f"teqc -sep sbf -E -S -C -J -I -week {gpsweeks[i]} -O.int {gpsinterval} "
                "-O.at SEPPOLANT_X_MF "
                "-O.rn 4501365 "
                '-O.rt "SEPT POLARX5" '
                '-O.mn "ARC09" '
                '-O.rv 1.0 '
                '-O.an SEPPOLANT_X_MF '
                f'-O.pe {antenna_height} '
                f'0 0 {file_path} > {out_filenames[i]}')   #info on command at https://www.unavco.org/software/data-processing/teqc/tutorial/tutorial.html
    commands.append(command1)
    
commands_array = np.array(commands)    

with open("/Users/home/whitefar/DATA/code/septentrio_to_rinex_GNSS_commands.sh","w") as f:
    for command in commands:
        f.write(command + "\n")
    

    
command2 = f"teqc +qc -report {out_filenames[i]} > {out_filenames[i]}.qc"
commands.append("/Users/home/whitefar/DATA/code/septentrio_to_rinex_GNSS_commands.sh",command2) 


#for i,file_path in enumerate(files_paths[:2]): 
#    
#    command1 = f'echo -sep sbf -E -S -C -J -I -week {gpsweeks[i]}\
#            -O.int {gpsinterval}\
#            -O.at SEPPOLANT_X_MF\
#            -O.rn 4501365 \
#            -O.rt "SEPT POLARX5" \
#            -O.mn "ARC09" \
#            -O.rv 1.0 \
#            -O.an SEPPOLANT_X_MF \
#            -O.pe {antenna_height}\
#            0 0 {file_path} > {out_filenames[i]}'    #info on command at https://www.unavco.org/software/data-processing/teqc/tutorial/tutorial.html
#    subprocess.call(command1,shell=True) 
#    command2 = f"teqc +qc -report {out_filenames[i]} > {out_filenames[i]}.qc"
#    subprocess.call(command2,shell=True) 
