#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 09:38:17 2020

@author: arran
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
from numpy import linalg as LA
from functools import reduce
import os
import sys
import time
import datetime as dt
import pandas as pd
import geopandas as gpd
import scipy as sp
from scipy import signal

from geopandas import GeoDataFrame
from shapely.geometry import Point
import fiona
from tqdm import tqdm

path = '/Volumes/arc_02/whitefar/DATA/REMOTE_SENSING/ICESAT2/Outputs'


# glob to list of files (run block of code creating wd and path variables if starting processing here)
ATL06_list = sorted(glob.glob(path+'/*.h5'))



filename = ATL06_list[5]
with h5py.File(filename, 'r') as f:
    # List all groups
    pairs=[1, 2, 3]
    beams=['l','r']
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

# /gt1l/land_ice_segments/h_li 
# gt2l
# gt2r
# gt3l
# gt3r
# dict containing data entries to retrive
dataset_dict = {'land_ice_segments':['delta_time','longitude','h_li']}

## the data can be converted to geopandas dataframe, see ATL08_2_gdf function in topolib gda_lib
temp_gdf = gda_lib.ATL06_2_gdf(ATL06_list[0],dataset_dict)
#enter 'r'






# `arrandf = df.query(expr="x > @minx & x < @maxx & y > @miny & y < @maxy")`
# x_atc is the 'x' along track coordinate;
#  x and y are EPSG:3031 coordinates; 
# h1, h2, h3, h4, h5 are the heights at each flyover cycle (there's quite a few NaNs);
# and hrange is just the maximum height difference within h1-h5.

# Oh, and the rough dates are:
# Cycle 1 - Spring2018 - 13Oct2018 - 28Dec2018
# Cycle 2 - Summer2019 - 28Dec2018 - 29Mar2019
# Cycle 3 - Autumn2019 - 29Mar2019 - 28Jun2019
# Cycle 4 - Winter2019 - 09Jul2019 - 26Sep2019
# Cycle 5 - Spring2019 - 26Sep2019 - 26Dec2019

df = pd.read_csv(filepath_or_buffer="/Volumes/arc_02/REMOTE_SENSING/ICESAT2/kamb_atl11_subset.csv")