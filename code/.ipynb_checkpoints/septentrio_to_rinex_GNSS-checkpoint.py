#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 11:54:01 2020

@author: whitefar
"""

import numpy as np
import glob
import os
import time
import datetime as dt
import subprocess
import gnsscal


# 

# =============================================================================
#INPUT

input_files = f"/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/ARC09/192.168.3.1/DSK1/SSN/LOG4_KIS1Hz/**/KIS*"

#Sample filename = KIS1352a.19_

output_file_location = "/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/"

gpsinterval = 30 # GPS sampling interval in secs
antenna_height = 1.55 # in metres

# =============================================================================
# Get information for processing

#find all the input files
files_paths = glob.glob(input_files)

filenames = [os.path.basename(file_path) for file_path in files_paths]

dates = []
datetime_strings = []
dudfile_index = []

#get datetime info out of each filename
for i,filename in enumerate(filenames):
    
    #find files with 0 entrys for a name
    if filename[4:8] == '0000':
        print(filename + " timestamp is zero, will skip this file")
        dudfile_index.append(i)
        continue
    
    daymonth = dt.date.fromordinal(int(filename[4:7]))
    year = int(dt.datetime.strptime(filename[-3:-1], '%y').strftime('%Y'))   #reads the date eg '19', outputs an integer 2019
    date = daymonth.replace(year=year)
    dates.append(date)
    
    time = dt.time(ord(filename[7])-97) #this is weird because the time is given as a letter ie a = 00:00-01:00 b = next hour
    datetime = dt.datetime.combine(date,time)
    datetime_string = datetime.isoformat()[:-6] # a string that looks like this 2020-12-10T21
    datetime_strings.append(datetime_string)

#delete the dud files from the list  
for i in dudfile_index:
    del files_paths[i], filenames[i]
  
#make a list of filenames to output, with the datetime info in a decent format    
out_filenames = [output_file_location + datetime_string for datetime_string in datetime_strings]
 
gpsweeks = [gnsscal.date2gpswd(date)[0] for date in dates] #gpsweek is a number since some epoch, needed for converting to rinex



# =============================================================================

#Write rinex files
for i,file_path in enumerate(files_paths): 
    subprocess.call(f"teqc -sep sbf -E -S -C -J -I -week {gpsweeks[i]} -O.int {gpsinterval} "
                "-O.at SEPPOLANT_X_MF "
                "-O.rn 4501365 "
                '-O.rt "SEPT POLARX5" '
                '-O.mn "ARC09" '
                '-O.rv 1.0 '
                '-O.s M '
                '-O.an SEPPOLANT_X_MF '
                f'-O.pe {antenna_height} '
                f'0 0 {file_path} > {out_filenames[i]}',shell=True)   #info on command at https://www.unavco.org/software/data-processing/teqc/tutorial/tutorial.html

# =============================================================================

#concatenate rinex files
for date in list(set(dates)):
    print(f"trying for {date}")
    allfiles_for_date = [file for file in glob.glob(f"/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/"+ date.isoformat()+'T*')]
    
    #get rid of the files with .qc ending (or any other ending)
    files_for_date = [file for file in allfiles_for_date if os.path.splitext(file)[1]=='']
    
    #to put in the right order of time, they need to be in time order to concatenate
    time_of_file = np.array([int(file[-2:]) for file in  files_for_date])
    index_to_sort = [int(time) for time in time_of_file.argsort()] 
    
    sorted_files_for_date =  [files_for_date[i] for i in index_to_sort]
    
    string_of_files_for_date = ' '.join(sorted_files_for_date)
    
#    subprocess.call("teqc " + string_of_files_for_date + " > "
#                    "/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/" + date.isoformat(),shell=True)
    print(f"Rinex file written for {date}")
    print("teqc " + string_of_files_for_date + " > "
          "/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/" + date.isoformat())

# =============================================================================



#remove all the non concatenated files written in previous step.
for filename in out_filenames:
    os.remove(filename)
    
    # =============================================================================



#GET
#
#! Error ! 2019 Dec 10 20:59:58.000: session sampling interval decreased to 30.000000 seconds
#! Error ! 2019 Dec 10 20:59:59.000: session sampling interval decreased to 30.000000 seconds
