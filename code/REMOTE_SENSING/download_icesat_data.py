#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 15:27:44 2020

@author: arran

This follows the jupyter notebook 
https://github.com/icesat2py/icepyx/blob/master/doc/examples/ICESat-2_DEM_comparison_Colombia_working.ipynb
"""

#from icepyx import is2class as ipd
import os
import shutil
import h5py
import xarray as xr 
# depedencies
import getpass
#from topolib.subsetDat import subsetBBox;
from topolib import icesat2_data 
import glob
import rasterio
from topolib import gda_lib
#from topolib import dwnldArctic
import numpy as np
import geopandas as gpd
from multiprocessing import Pool
import contextily as ctx
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry.polygon import orient


get_ipython().run_line_magic('load_ext', 'autoreload')
from icepyx import is2class as ipd
get_ipython().run_line_magic('autoreload', '2')
#in order to use "as ipd", you have to use autoreload 2, which will automatically reload any module not excluded by being imported with %aimport -[module]



get_ipython().run_line_magic('cd', '~/software/icepyx/dev-notebooks')

shp_filepath = '~/study_area_buffer_geo.shp'



region_a = ipd.Icesat2Data('ATL06', shp_filepath,['2019-02-22','2019-02-28'],start_time='00:00:00', end_time='23:59:59')


earthdata_uid = 'arran'
email = 'arran.whiteford@vuw.ac.nz'
sessionr=region_a.earthdata_login(earthdata_uid, email)

region_a.avail_granules()

region_a.reqparams

region_a.granule_info

region_a.orderIDs

region_a.granules

path = '~/Test1/'

region_a.download_granules(session, path)

#Clean up Outputs folder by removing individual granule folders 
for root, dirs, files in os.walk(path, topdown=False):
    for file in files:
        try:
            shutil.move(os.path.join(root, file), path)
        except OSError:
            pass
        
for root, dirs, files in os.walk(path):
    for name in dirs:
        os.rmdir(os.path.join(root, name))
        
#preprocess        
        
#Convert data into geopandas dataframe, which allows for doing basing geospatial opertaions        
# glob to list of files (run block of code creating wd and path variables if starting processing here)
ATL06_list = sorted(glob.glob(path+'/*.h5'))

filename = ATL06_list[5]
with h5py.File(filename, 'r') as f:
    # List all groups
    pairs=[1, 2, 3]
    beams=['l','r']
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]
#

ATL06_list



        
        
