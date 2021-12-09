#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:47:47 2020

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

from geopandas import GeoDataFrame
from shapely.geometry import Point
import fiona


path = "/Users/home/whitefar/DATA/MINE/Holdsworth_to_Thorndon_26Jan20.gpx"
ht = gpd.read_file(path,layer='track_points')
ht.track_points['datetime'] = np.array([dt.datetime.strptime(t,"%Y-%m-%dT%H:%M:%S") for t in ht.track_points.time])
ht.track_points['timestamp'] = ht.track_points.datetime.apply(ts_func)