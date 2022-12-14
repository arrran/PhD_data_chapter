#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:56:42 2020

@author: arran
"""

import rasterio as rio
import rasterio.mask
import fiona
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, LineString
from shapely.affinity import scale
import numpy as np
import os
import sys
import glob
import matplotlib.pyplot as plt
from scipy import interpolate
import xarray as xr
from shapely.geometry import Point
from shapely.geometry import LineString
import hvplot.xarray  # noqa
import cartopy.crs as crs

path = "/Volumes/arc_02/REMOTE_SENSING/ICESAT2/ds_subset_kamb_20200404.nc"

class icesat_dataset:
    """
    
    """
    
    def __init__(self,path):
        path = "/Volumes/arc_02/REMOTE_SENSING/ICESAT2/ds_subset_kamb_20200404.nc"
        ds = xr.open_dataset(path, engine="h5netcdf")
        
        df = ds[["h_corr", "utc_time","h_corr_sigma"]].to_dataframe().reset_index()
        
        #dataframe over channel area
        da = df.query("x > -382064.5722209641 & x < -374781.1654740692 & y > -734075.0820404041 & y < -722764.4514729496")
        
        #convert to geodataframe
        points = [Point(xy) for xy in zip(da.x,da.y)]
        gda = gpd.GeoDataFrame(da,geometry=points,crs=3031)
        
        ts_func = lambda t : pd.Timestamp.timestamp(t)
        
        # gda['timestamp'] = gda_line.utc_time.apply(ts_func)
               
        gda_lines = {}
        
        #icesat lines
        is_lines = {'is0': LineString([(-377164.7979570718, -734049.2067410639),(-374791.4999125675, -731478.0037744837)]),
                         'is1': LineString([(-381896.7461368341, -733380.6368849602),(-377024.528257302, -728195.0649845804)]),
                         'is3': LineString([(-376096.6933692924, -732027.5543375007), (-375737.272743314, -731644.9631040762)]),
                         'is4': LineString([(-381896.7461368341, -733380.6368849602), ( -377024.528257302, -728195.0649845804)]),
                         'is5': LineString([(-380587.3777984626, -722783.2054143124), ( -382058.8039400483, -729564.8165617182)]),
                         'is7': LineString([(-376691.5694304057, -722796.118348052), ( -378990.8284925717, -734052.5693451025)]),
                         'is8': LineString([(-376974.9258506759, -722791.6817553517), ( -379329.5927148979, -734035.883147971)]),
                         'is9': LineString([(-379830.9172270952, -726358.755056211), ( -377754.301185295, -724148.9706686384)]),
                         'is11': LineString([(-376096.6933692924, -732027.5543375007), ( -375737.272743314, -731644.9631040762)]),
                         'is12': LineString([(-376691.5694304057, -722796.118348052), ( -378990.8284925717, -734052.5693451025)]),
                         'is13': LineString([(-380331.1736181656, -722803.0102771204), ( -382053.63229378, -731022.4492230047)]),
                         'is14': LineString([(-380587.3777984626, -722783.2054143124), ( -382058.8039400483, -729564.8165617182)])}
        #extrapolate short lines
        is_lines['is5'] = scale(is_lines['is5'],20,20)
        is_lines['is4'] = scale(is_lines['is4'],2,2)
        is_lines['is9'] = scale(is_lines['is9'],4,4)
        
    def line(self,icesat_line_number):
        """

        Parameters
        ----------
        icesat_line_number : TYPE
            as defined in 'ICESAT2_fast_analysis_new_dataset.ipynb'
        
        Sets and Returns
        -------
        line: a geodataframe of icesat points over that icesat line

        """
        
        icesat_line_number = int(icesat_line_number)
        
        dict_entry = f'is{icesat_line_number}'
        
        line = is_lines[dict_entry]
        
        #dataset points over an icesat line
        gda_line = gda[gda.geometry.intersects(line.buffer(2.5))]
        
        gda_lines[dict_entry] = gda[gda.geometry.intersects(line.buffer(2.5))]
        
        list_of_cycle_numbers = gda_line.cycle_number.value_counts().index.to_list()
        
        print(f'icesat line {icesat_line_number}')
        print('cycle number, corresponding number of points')
        print(gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts())
        
        return line
        
    def all_lines(self):
        """
        do line() over all lines
        """
        for line in list(is_lines.values()):
            gda_lines[dict_entry] = gda[gda.geometry.intersects(line.buffer(2.5))]
            
    def plot_line_crosssection(self,icesat_line_number):
        """
        Parameters
        ----------
        icesat_line_number : TYPE
            DESCRIPTION.


        Plots
        -------
        cross sections of surface for all cycle numbers

        """
        try:
            gda_line = gda_line
        except AttributeError:
            dict_entry = f'is{icesat_line_number}'
            gda_line = gda_lines[dict_entry]
                        
        
        list_of_cycle_numbers = gda_line.cycle_number.value_counts().index.to_list()
        
        
        
        for cycle_number in list_of_cycle_numbers:
            plt.plot(gda_line[gda_line.cycle_number==cycle_number].geometry.x,gda_line[gda_line.cycle_number==cycle_number].h_corr,'1')
        plt.legend( list(map(str,list_of_cycle_numbers)) )
        plt.xlabel('x, (m)')
        plt.ylabel('elevation, (m)')
        plt.grid()
        plt.show()
        
    def plot_dhdt_crosssection(self,cycle_number_from,cycle_number_till,icesat_line_number):
        """

        Parameters
        ----------
        cycle_number_from : TYPE
            DESCRIPTION.
        cycle_number_till : TYPE
            DESCRIPTION.
        icesat_line : TYPE
            DESCRIPTION.
            
        Plots
        -------
        cross sections of surface and change in elevation between two cycles.

        """
        

        
        try:
            gda_line = gda_line
        except AttributeError:
            dict_entry = f'is{icesat_line_number}'
            gda_line = gda_lines[dict_entry]
            
        for cycle in [cycle_number_from,cycle_number_till]:
            try:
                gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts().loc[cycle]
            except KeyError:
                print(f'There is no data for line {icesat_line_number}, cycle {cycle}')
            
       
        gda_line_diff = gpd.GeoDataFrame( gda_line[gda_line.cycle_number==cycle_number_from],geometry=gda_line.geometry )
        # gda_line_diff.rename(columns{"h_corr": f"h_corr_cycle_{cycle_number_from}"}, inplace=True)
        # gda_line_diff[f"h_corr_cycle_{cycle_number_till}"] = gda_line[gda_line.cycle_number==cycle_number_till].h_corr
       

        #DH MUST ONLY SUBTRACT THE SAME POINTS
                                  
        gda_line_diff['dh'] =(gda_line[gda_line.cycle_number==cycle_number_from].h_corr.to_numpy() -
                              gda_line[gda_line.cycle_number==cycle_number_till].h_corr.to_numpy())
        
        # get the time in years between data points
        gda_line_diff['dt'] = (gda_line[gda_line.cycle_number==cycle_number_till].utc_time.to_numpy()  - 
                          gda_line[gda_line.cycle_number==cycle_number_from].utc_time.to_numpy()  ) 
        gda_line_diff.dt = gda_line_diff.dt /  np.timedelta64(1, 'Y')
        
        gda_line_diff['dhdt'] = gda_line_diff['dh'].to_numpy()/gda_line_diff['dt'].to_numpy()
        
        
        fig, ax1 = plt.subplots(figsize=(12,7),dpi=300)
        color = 'tab:blue'
        ax1.set_xlabel('x',fontsize=18)
        ax1.set_ylabel('elevation (m)', color=color,fontsize=18)
        
        ax1.plot(gda_line_diff[gda_line_diff.cycle_number==cycle_number_from].x,
                 gda_line_diff[gda_line_diff.cycle_number==cycle_number_from].h_corr,'b-')
        
        ax1.plot(gda_line_diff[gda_line_diff.cycle_number==cycle_number_till].x,
                 gda_line_diff[gda_line_diff.cycle_number==cycle_number_till].h_corr,'g-')
        
        ax1.set_label([f"{gda_line[gda_line.cycle_number==cycle_number_from].utc_time.mean().date()}",
                   f"{gda_line[gda_line.cycle_number==cycle_number_till].utc_time.mean().date()}"
                   ])
        
        ax1.tick_params(axis='x', labelcolor=color)
        ax1.grid()
        
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        
        color = 'tab:red'
        ax2.set_ylabel('rate of change in surface , (m/a)', color=color,fontsize=18)  # we already handled the x-label with ax1
        ax2.plot(gda_line_diff.x,gda_line_diff.dhdt,'ro')
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.grid()
        
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.title(f"Rate of change in surface from {gda_line[gda_line.cycle_number==cycle_number_from].utc_time.mean().date()}" 
                  f" to {gda_line[gda_line.cycle_number==cycle_number_till].utc_time.mean().date()}")
        plt.show()
        
    def plot_dhdt_map(self,cycle_number_from,cycle_number_till):
        """
        

        Parameters
        ----------
        cycle_number_from : TYPE
            DESCRIPTION.
        cycle_number_till : TYPE
            DESCRIPTION.

        Plots
        -------
        Plan view scatter plot of 

        """
        
        dadh = da
        
        dadh['dh'] =  ( dadh[cycle_number==cycle_number_from].h_corr -
                          dadh(cycle_number==cycle_number_till).h_corr )
        
        dadh['dt'] = (dadh[cycle_number==cycle_number_from].utc_time 
                          - dadh(cycle_number==cycle_number_till).utc_time )
        
        dadh['dhdt'] = dadh['dh'].to_numpy()/dadh['dt'].to_numpy()
        
        plt.figure(figsize=(12,12))
        plt.scatter(dadh.x,dadh.y,c=dadh.dhdt,cmap='winter')
        plt.legend(['melt_rate'])
        plt.colorbar()
        plt.grid()
        plt.show()
        
        
        ts_func = lambda t : pd.Timestamp.timestamp(t)
        
        gda_line.utc_time.apply(ts_func)
        
    
    