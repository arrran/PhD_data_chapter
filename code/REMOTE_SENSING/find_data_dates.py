#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 09:41:23 2020

@author: arran
"""
import numpy as np
import os
import sys

import matplotlib.pyplot as plt

import pandas as pd


wd = '/Users/home/whitefar/DATA/REMA_2m_strips/KAMB_CHANNEL/'

sys.path.append(os.path.abspath(wd))
sys.path.append(os.path.abspath('/Users/home/whitefar/DATA/code/'))
os.chdir(wd)

# attribute_table_stripes_over_channel
df = pd.read_csv('/Users/home/whitefar/DATA/REMA_2m_strips/KAMB_CHANNEL/attribute_table_stripes_over_channel.txt',delimiter='\t')

#acquisition date
df.acquisition
df.acquisit_1
# both are the same

#stripe name:
df.name

