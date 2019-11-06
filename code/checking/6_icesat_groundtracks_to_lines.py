#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 15:21:15 2019

@author: whitefar
"""

#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import zipfile

import os
import sys
import geopandas as gpd
import shapely as sp
import matplotlib.pyplot as plt
import fiona
import shutil

fiona.supported_drivers['libkml'] = 'r'
fiona.supported_drivers['LIBKML'] = 'r'

#working_dir = '/Volumes/arc_02/REMOTE_SENSING/'
working_dir = "/Users/home/whitefar/DATA/Icesat2/Reference_ground_tracks/"

sys.path.append(os.path.abspath(working_dir))
sys.path.append(os.path.abspath('/Users/home/whitefar/DATA/code/'))
os.chdir(working_dir)


# In[4]:


zipped_kml_path = "/Volumes/arc_02/REMOTE_SENSING/ICESAT2/REFERENCE_GROUND_TRACKS/ORIGINALS/"

zipped_kml_fnames = [f for f in os.listdir(zipped_kml_path) if os.path.isfile(os.path.join(zipped_kml_path, f))]

temp_kml_location = "/Volumes/arc_02/REMOTE_SENSING/ICESAT2/REFERENCE_GROUND_TRACKS/TEMP/" #save them temporararily

output_folder = "/Volumes/arc_02/REMOTE_SENSING/ICESAT2/REFERENCE_GROUND_TRACKS/SHAPEFILES/"


zipped_kml_fnames


# In[18]:


zipped_kml_fname = zipped_kml_fnames[2]
zipped_kml_fname


# In[20]:


def make_reference_tracks(zipped_kml_path, zipped_kml_fname,output_folder,temp_kml_location):
    """starts with zip file, outputs  shapefile lines"""


    
    # Unzip the folder of kmls
    zip_fname = zipped_kml_fname
    unzipped_folder = zip_fname[:-4]
    
    print("zip_fname =",zip_fname)

    with zipfile.ZipFile(zipped_kml_path + zipped_kml_fname,"r") as zip_ref:
        zip_ref.extractall(temp_kml_location)
        
    kml_folder = unzipped_folder
    if not os.path.exists(temp_kml_location + kml_folder):
        kml_folder = os.listdir(temp_kml_location)[0]
        if kml_folder == "__MACOSX":
            kml_folder = os.listdir(temp_kml_location)[1]
    if os.path.exists(temp_kml_location + kml_folder):
        print("folder for kmls is " + temp_kml_location + kml_folder)
    else:
        print("cant find folder for kmls ", os.listdir(temp_kml_folder))
    #print("path for kmls is " + temp_kml_path)
    
    temp_kml_path = temp_kml_location + kml_folder  + '/'

    #get list of track names from that folder
    track_kml_fnames = [f for f in os.listdir(temp_kml_path) if os.path.isfile(os.path.join(temp_kml_path, f)) and f[-4:]==".kml"]
    print("length of track_kml_fnames is",len(track_kml_fnames))
    
    output_path = output_folder + kml_folder  + '/'
    
    #make a directory to put output files in
    os.mkdir(output_path)
    print("output is ",output_path)

    #iterate over that folder, making shp files
    for track in track_kml_fnames:        
        gdf = gpd.read_file(temp_kml_path  + track)
        gdf.to_file( output_path  + track[:-4] + ".shp")
        print("written to ",output_path  + track[:-4] + ".shp")
        
    shutil.rmtree(temp_kml_location)
    print("temp files deleted")
    os.mkdir(temp_kml_location)
    
    


# In[14]:



    


# In[15]:


for fname in zipped_kml_fnames:
    
    make_reference_tracks(zipped_kml_path,fname,output_folder,temp_kml_location)


# In[19]:


shutil.rmtree("/Volumes/arc_02/REMOTE_SENSING/ICESAT2/REFERENCE_GROUND_TRACKS/TEMP")
shutil.rmtree("/Volumes/arc_02/REMOTE_SENSING/ICESAT2/REFERENCE_GROUND_TRACKS/SHAPEFILES")
os.mkdir("/Volumes/arc_02/REMOTE_SENSING/ICESAT2/REFERENCE_GROUND_TRACKS/TEMP")
os.mkdir("/Volumes/arc_02/REMOTE_SENSING/ICESAT2/REFERENCE_GROUND_TRACKS/SHAPEFILES")


# In[ ]:


# to reset


# In[ ]:





# In[ ]:





# In[ ]:




