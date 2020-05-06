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

shp_filepath = '/Users/home/whitefar/DATA/REMOTE_SENSING/REMA_2m_strips/study_area_buffer_geo.shp'

study_area = gpd.read_file(shp_filepath)


region_a = ipd.Icesat2Data('ATL06',shp_filepath ,['2000-10-14','2019-11-15'],start_time='00:00:00', end_time='23:59:59')

# (-411054.19240523444,
#  -739741.7702261859,
#  -365489.6822096751,
#  -699564.516934089)
# (-153.4003855308609, -82.65774797172867, -150.10403144247186, -82.34332209648137)



earthdata_uid = 'whitefar'
email = 'arran.whiteford@vuw.ac.nz'
# pswd = 'Whitefar44'

session = region_a.earthdata_login(earthdata_uid, email)

region_a.avail_granules()

region_a.order_granules(session,subset=False)


path = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/ICESAT2/Outputs'

region_a.download_granules(session,path)


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
        




        
        
