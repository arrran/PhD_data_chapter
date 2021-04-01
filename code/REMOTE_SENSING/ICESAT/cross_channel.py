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




class icesat_dataset:
    """
    This class will load an iceset database (provided by wei ji at https://github.com/weiji14/deepicedrain/blob/master/atl11_play.py with kamb as placeset subset) and will produce subsets of points over specific lines
    in the kamb channel area to compare over the years.
    
    It is supposed to streamline and speed up analysis so that it is easy to analyse new data which comes out
    
    """
    
    def __init__(self,path):
        """
        
        Parameters
        ----------
        path : string of location of dataset,
        dataset is a .nc file

        Sets
        ----------
        self.gda a geodataframe of the icesat points over the kamb channel area
        self.is_lines a dictionary of icesat repeat lines (no data)

        """
        
        ds = xr.open_dataset(path, engine="h5netcdf")
        
        df = ds[["h_corr", "utc_time","h_corr_sigma"]].to_dataframe().reset_index()
        
        #dataframe over channel area
        da = df.query("x > -382064.5722209641 & x < -374781.1654740692 & y > -734075.0820404041 & y < -722764.4514729496")
        
        #convert to geodataframe
        points = [Point(xy) for xy in zip(da.x,da.y)]
        self.gda = gpd.GeoDataFrame(da,geometry=points,crs=3031)
        
        self.gda_lines = {}
        
        #icesat lines
        self.is_lines = {'is0': LineString([(-377164.7979570718, -734049.2067410639),(-374791.4999125675, -731478.0037744837)]),
                         'is2': LineString([(-376096.6933692924, -732027.5543375007), (-375737.272743314, -731644.9631040762)]),
                         'is3': LineString([(-381706.7103145397, -733633.1081130945),(-375189.2166156959, -726633.5887804579)]),
                         'is4': LineString([(-381896.7461368341, -733380.6368849602), ( -377024.528257302, -728195.0649845804)]),
                         'is5': LineString([(-378187.3930402951, -729042.2059732673), (-377786.1458335032, -728618.909909106)]),
                         'is6': LineString([(-376438, -723116), (-378537,-733580 )]),
                         'is7': LineString([(-376691.5694304057, -722796.118348052), ( -378990.8284925717, -734052.5693451025)]),
                         'is8': LineString([(-376974.9258506759, -722791.6817553517), ( -379329.5927148979, -734035.883147971)]),
                         'is9': LineString([(-379830.9172270952, -726358.755056211), ( -377754.301185295, -724148.9706686384)]),
                         'is11': LineString([(-382062.0207341266, -727907.9477563774), (-377549.1232601753, -723185.2876206436)]),
                         'is12': LineString([(-380065.9441702637, -722818.7788760758), (-382061.9004570942, -732590.3699018088)]),
                         'is13': LineString([(-380331.1736181656, -722803.0102771204), ( -382053.63229378, -731022.4492230047)]),
                         'is14': LineString([(-380587.3777984626, -722783.2054143124), ( -382058.8039400483, -729564.8165617182)])}
        #extrapolate short lines
        self.is_lines['is2'] = scale(self.is_lines['is2'],15,15)
        self.is_lines['is3'] = scale(self.is_lines['is3'],1.5,1.5)
        self.is_lines['is5'] = scale(self.is_lines['is5'],20,20)
        self.is_lines['is6'] = scale(self.is_lines['is6'],1.5,1.5)
        self.is_lines['is4'] = scale(self.is_lines['is4'],2,2)
        self.is_lines['is9'] = scale(self.is_lines['is9'],4,4)
        self.is_lines['is7'] = scale(self.is_lines['is7'],1.2,1.2)
        self.is_lines['is11'] = scale(self.is_lines['is11'],1.5,1.5)
        
        self.meta = pd.DataFrame({})
        
        self.channel_mid = gpd.read_file("/Volumes/arc_04/FIELD_DATA/K8621920/RES/SHARE_SHAPEFILES/interpolated_channel_mid.shp")
        
        print(f'For this dataset')
        print('cycle number, corresponding number of points')
        print(self.gda.cycle_number.value_counts())
        print('cycle number, corresponding number of non NaN points')
        print(self.gda.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts())
        
    def plot_icesat_lines_map(self):
        """
        Plots
        -------
        a map of the all the icesat lines, to see whats where

        """
        
        gd_chan = gpd.read_file("/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RES/PROCESSED_LINES_GISFILE/linedownchan.shp")
        
        plt.figure(figsize=(12,10),dpi=500)
        for key,value in zip(self.is_lines.keys(),self.is_lines.values()):
            plt.plot(value.xy[0],value.xy[1],label=key)
            plt.legend()
        plt.plot(self.channel_mid.geometry.x,self.channel_mid.geometry.y,'g:',label='channel centre')
        plt.xlim([-382064.5722209641,-374781.1654740692])
        plt.ylim([-734075.0820404041, -722764.4514729496])
        plt.grid()
        
    def getdata_line(self,icesat_line_number,return_line=False,verb=True,buff=2.5):
        """

        Parameters
        ----------
        icesat_line_number : TYPE
            as defined in 'ICESAT2_fast_analysis_new_dataset.ipynb'
        
        return_line : if set to true, the function will return a geodataframe of all icesat data close to the line
        
        verb : if set to true, spits out info about how much data on the line
        
        Sets and Returns
        -------
        self.gda_lines[line_name]
        gda_line: a geodataframe of icesat points over that icesat line

        """
        
        self.icesat_line_number = int(icesat_line_number)
        
        dict_entry = f'is{icesat_line_number}'
        
        try:
            self.line = self.is_lines[dict_entry]
        except KeyError:
            print(f'There is no icesat line {self.icesat_line_number}')
            return
        
        #dataset points over an icesat line
        gda_line = self.gda[self.gda.geometry.intersects(self.line.buffer(buff))]
        
        self.gda_lines[dict_entry] =self.gda[self.gda.geometry.intersects(self.line.buffer(buff))]
        
        if verb==True:    
            print(f'icesat line {self.icesat_line_number}')
            print('cycle number, corresponding number of points')
            print(gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts())
            
        self.meta = self.meta.append( pd.DataFrame({'cycle_number' : gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts().index.to_list(),
                                           'num_of_data': gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts().to_list(),
                                           'icesat_line_number': [self.icesat_line_number]*len(gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts().to_list())
                                           })
                         , ignore_index=True)
        
        self.meta_table = pd.pivot_table(self.meta, values='num_of_data', index=['cycle_number'],
                                         columns=['icesat_line_number'], aggfunc=np.sum)
        
        if return_line==True:
            return gda_line
        
    def getdata_alllines(self,verb=True,buff=2.5):
        """
        does getdata_line() over all lines
        
        i.e. gets all iceset data over all the icesat lines in study area 
        
        
        verb : if set to true, spits out info about how much data on each line
        optional buff=2.5, buffer size
        """
        for line in list(self.is_lines.keys()):
            self.getdata_line(line[2:],verb=False,buff=buff)

        print(self.meta_table)
        
    def getdata_dhdt(self,cycle_number_from,cycle_number_till):
        
        self.dadh_cycles = [cycle_number_from,cycle_number_till]
        
               # This dataframe has h_corr from cycle_number_from, and dhdt. Times are from cycle_number_from
        self.dadh = gpd.GeoDataFrame( self.gda[self.gda.cycle_number==cycle_number_from],geometry=self.gda.geometry )
                    
        self.dadh['dh'] = (self.gda[self.gda.cycle_number==cycle_number_till].h_corr.to_numpy() - 
                              self.gda[self.gda.cycle_number==cycle_number_from].h_corr.to_numpy())
        
        # get the time in years between data points
        self.dadh['dt'] = (self.gda[self.gda.cycle_number==cycle_number_till].utc_time.to_numpy()  - 
                          self.gda[self.gda.cycle_number==cycle_number_from].utc_time.to_numpy()  ) 
        self.dadh.dt = self.dadh.dt /  np.timedelta64(1, 'Y')
        
        self.dadh['dhdt'] = self.dadh['dh'].to_numpy() / self.dadh['dt'].to_numpy()  
        
        
    def plot_line_crosssection(self,icesat_line_number):
        """
        Parameters
        ----------
        icesat_line_number 

        Plots
        -------
        cross sections of surface for all cycle numbers

        """
        
        # if data over that line hasnt been found using getdataline or all_lines, find it.
        dict_entry = f'is{icesat_line_number}'
        try:
            gda_line = self.gda_lines[dict_entry]
        except KeyError:
            self.getdata_line(icesat_line_number)
            gda_line = self.gda_lines[dict_entry]
        
        list_of_cycle_numbers = gda_line.cycle_number.value_counts().index.to_list()
        legend = ["cycle number {c_n}" for c_n in list_of_cycle_numbers]
        
        plt.figure(figsize=(12,7),dpi=500)
        for cycle_number in list_of_cycle_numbers:
            plt.plot(gda_line[gda_line.cycle_number==cycle_number].sort_values(by ='x' ).geometry.x,gda_line[gda_line.cycle_number==cycle_number].sort_values(by ='x' ).h_corr,'1')
        plt.legend(legend )
        plt.xlabel('x, (m)')
        plt.ylabel('elevation, (m)')
        plt.grid()
        plt.show()
        
    def plot_dhdt_crosssection(self,cycle_number_from,cycle_number_till,icesat_line_number):
        """

        Parameters
        ----------
        cycle_number_from : 
        cycle_number_till : 
        icesat_line_number : 
            
        Plots
        -------
        cross sections of surface and change in elevation between two cycles.

        """        
        
        # if data over that line hasnt been found using getdataline or all_lines, find it.
        dict_entry = f'is{icesat_line_number}'
        try:
            gda_line = self.gda_lines[dict_entry]
        except KeyError:
            self.getdata_line(icesat_line_number)
            gda_line = self.gda_lines[dict_entry]
            
        for cycle in [cycle_number_from,cycle_number_till]:
            try:
                gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts().loc[cycle]
            except KeyError:
                print(f'There is no data for line {icesat_line_number}, cycle {cycle}')
                return
            
        point_midchannel = np.array( LineString(gda_line.geometry.tolist()).intersection(LineString(self.channel_mid.geometry.tolist())) )
        
        if len(point_midchannel)>1:
            point_midchannel = np.mean(point_midchannel,axis=0)
        

        # This dataframe has h_corr from cycle_number_from, and dhdt. Times are from cycle_number_from
        gda_line_diff = gpd.GeoDataFrame( gda_line[gda_line.cycle_number==cycle_number_from],geometry=gda_line.geometry ).reset_index(drop=True)
        gda_line_diff.rename(columns={"h_corr": f"h_corr_cycle_{cycle_number_from}"}, inplace=True)
        gda_line_diff[f"h_corr_cycle_{cycle_number_till}"] = gda_line[gda_line.cycle_number==cycle_number_till].h_corr.reset_index(drop=True)
       
                    
        gda_line_diff['dh'] =- (gda_line[gda_line.cycle_number==cycle_number_from].h_corr.to_numpy() -
                              gda_line[gda_line.cycle_number==cycle_number_till].h_corr.to_numpy())
        
        # get the time in years between data points
        gda_line_diff['dt'] = (gda_line[gda_line.cycle_number==cycle_number_till].utc_time.to_numpy()  - 
                          gda_line[gda_line.cycle_number==cycle_number_from].utc_time.to_numpy()  ) 
        gda_line_diff.dt = gda_line_diff.dt /  np.timedelta64(1, 'Y')
        
        gda_line_diff['dhdt'] = gda_line_diff['dh'].to_numpy()/gda_line_diff['dt'].to_numpy()
                

        fig, ax1 = plt.subplots(figsize=(12,7),dpi=500)
        color = 'tab:blue'
        ax1.set_xlabel('x',fontsize=18)
        ax1.set_ylabel('elevation (m)', color=color,fontsize=18)
        
        
        
        ax1.plot(gda_line_diff.sort_values(by ='x' ).x,
                 gda_line_diff.sort_values(by ='x' )[f"h_corr_cycle_{cycle_number_from}"],'b-')
        
        ax1.plot(gda_line_diff.sort_values(by ='x' ).x,
                 gda_line_diff.sort_values(by ='x' )[f"h_corr_cycle_{cycle_number_till}"],'g-')
        
        ax1.set_label([f"{gda_line[gda_line.cycle_number==cycle_number_from].utc_time.mean().date()}",
                   f"{gda_line[gda_line.cycle_number==cycle_number_till].utc_time.mean().date()}"
                   ])
        
        ax1.tick_params(axis='x', labelcolor=color)
        ax1.legend([f"{gda_line[gda_line.cycle_number==cycle_number_from].utc_time.mean().date()}",
                   f"{gda_line[gda_line.cycle_number==cycle_number_till].utc_time.mean().date()}"
                   ])
        ax1.grid()
        
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        
        color = 'tab:red'
        ax2.set_ylabel('rate of change in surface , (m/a)', color=color,fontsize=18)  # we already handled the x-label with ax1
        ax2.plot(gda_line_diff.x,gda_line_diff.dhdt,'ro')
        ax2.vlines(point_midchannel[0],gda_line_diff.dhdt.min(),gda_line_diff.dhdt.max(),label='Channel midpoint')
        
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.grid()
        
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.title(f"Rate of change in surface from {gda_line[gda_line.cycle_number==cycle_number_from].utc_time.mean().date()}" 
                  f" to {gda_line[gda_line.cycle_number==cycle_number_till].utc_time.mean().date()}")
        plt.show()
        
    def plot_dhdt_map(self,cycle_number_from,cycle_number_till,**kwargs):
        """
        
        Parameters
        ----------
        cycle_number_from : 
        cycle_number_till : 

        Plots
        -------
        Plan view scatter plot of dhdt

        """
        
            # check to see if dhdt has been calculated with the came cycle numbers
        try:
            if not self.dadh_cycles == [cycle_number_from,cycle_number_till]:
                self.getdata_dhdt(cycle_number_from,cycle_number_till)         
        except AttributeError:
            self.getdata_dhdt(cycle_number_from,cycle_number_till)
        
        plt.figure(figsize=(15,15))
        plt.scatter(self.dadh.x,self.dadh.y,c=self.dadh.dhdt,cmap='bwr_r',**kwargs)
        plt.plot(self.channel_mid.geometry.x,self.channel_mid.geometry.y,'g:',label='channel centre')
        plt.legend([f'surface channel low',f'melt_rate from cycle {cycle_number_from} to {cycle_number_till} '])
        cb = plt.colorbar()
        cb.set_label('rate of elevation change, m/a')
        plt.grid()
        # plt.show()
        
    def plot_dhdt_crosssection_map(self,cycle_number_from,cycle_number_till,icesat_line_number,buff=2.5):
        """

        This is just for troubleshooting, see what points the buffer region picks up

        """
        
        
        # if data over that line hasnt been found using getdataline or all_lines, find it.
        dict_entry = f'is{icesat_line_number}'
        try:
            gda_line = self.gda_lines[dict_entry]
        except KeyError:
            self.getdata_line(icesat_line_number)
            
            gda_line = self.gda_lines[dict_entry]
            
        for cycle in [cycle_number_from,cycle_number_till]:
            try:
                gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts().loc[cycle]
            except KeyError:
                print(f'There is no data for line {icesat_line_number}, cycle {cycle}')
                return
            
       
            # check to see if dhdt has been calculated with the came cycle numbers
        try:
            if not self.dadh_cycles == [cycle_number_from,cycle_number_till]:
                self.getdata_dhdt(cycle_number_from,cycle_number_till)         
        except AttributeError:
            self.getdata_dhdt(cycle_number_from,cycle_number_till)
        
        polygon = self.is_lines[dict_entry].buffer(buff)
        
        plt.figure(figsize=(15,15))
        plt.scatter(gda_line_diff.x,gda_line_diff.y,c=gda_line_diff.dhdt,cmap='bwr_r',vmin=-0.8, vmax=1,label ='melt_rate')
        plt.plot(polygon.exterior.xy[0],polygon.exterior.xy[1],label='is polygon')
        plt.plot(self.channel_mid.geometry.x,self.channel_mid.geometry.y,'g:',label='channel middle')
        plt.legend()
        cb = plt.colorbar()
        cb.set_label('rate of elevation change, m/a')
        plt.grid()
        plt.show()
    

    def plot_multiple_dhdt_crosssection_map(self,cycle_number_from,cycle_number_till,icesat_lines,buff=300,
                                            **kwargs):
        """
        

        Parameters
        ----------
        cycle_number_from : TYPE integer
            DESCRIPTION. icesat cycle number
        cycle_number_till : TYPE integer
            DESCRIPTION. icesat cycle number
        icesat_lines : TYPE list of lists of number and string 'up' or 'down'
            DESCRIPTION e.g. icesat_lines = [[12,'up'],[11,'down'],[7,'up'],[3,'up']]

        Returns
        -------
        None.
        
        Plot of the map of dhdt but with labelled lines

        """
            # check to see if dhdt has been calculated with the came cycle numbers
        try:
            if not self.dadh_cycles == [cycle_number_from,cycle_number_till]:
                self.getdata_dhdt(cycle_number_from,cycle_number_till)         
        except AttributeError:
            self.getdata_dhdt(cycle_number_from,cycle_number_till)
        
        dict_entrys = [[f'is{icesat_line_number[0]}',icesat_line_number[1]] for icesat_line_number in icesat_lines]
        
        polygons = [[self.is_lines[dict_entry[0]].buffer(buff),dict_entry[1]] for dict_entry in dict_entrys]

        plt.figure(figsize=(10,10))
        plt.scatter(self.dadh.x,self.dadh.y,c=self.dadh.dhdt,cmap='bwr_r',**kwargs)
        plt.plot(self.channel_mid.geometry.x,self.channel_mid.geometry.y,'g-')
        i=65
        for polygon, text_loc in polygons:   
            
            plt.plot(polygon.exterior.xy[0],polygon.exterior.xy[1],':',label='is polygon')
            if text_loc == 'up':
                plt.text(np.max(polygon.exterior.xy[0]),np.max(polygon.exterior.xy[1])+ buff/2,chr(i),size=30)
            elif text_loc == 'down':
                    plt.text(np.min(polygon.exterior.xy[0])- buff,np.min(polygon.exterior.xy[1])- buff,chr(i),size=30)
            i += 1
        plt.legend(['surface channel low','melt_rate'],fontsize=12)
        cb = plt.colorbar()
        cb.set_label('rate of elevation change, m/a',fontsize=12)
        plt.grid()
        plt.show()
                
        # icesat_line_numbers = [b[0] for b in icesat_lines]


    def plot_multiple_dhdt_crosssection(self,cycle_number_from,cycle_number_till,icesat_lines):
        """
        MUST DO  plot_multiple_dhdt_map() first.

        Parameters
        ----------
        cycle_number_from : TYPE integer
            DESCRIPTION. icesat cycle number
        cycle_number_till : TYPE integer
            DESCRIPTION. icesat cycle number
        icesat_lines : TYPE list of lists of number and string 'up' or 'down'
            DESCRIPTION e.g. icesat_lines = [[12,'up'],[11,'down'],[7,'up'],[3,'up']]

        Returns
        -------
        None.
        
        Plots of cross sections which correspond to the map made with  plot_multiple_dhdt_crosssection_map() 
        """

        # if data over that line hasnt been found using getdataline or all_lines, find it.
        dict_entrys = [f'is{icesat_line_number[0]}' for icesat_line_number in icesat_lines]
        icesat_line_numbers = [b[0] for b in icesat_lines]
        
        gda_line_diffs = []
        point_midchannels = []
        
        # for each line do the cross section diff
        for icesat_line_number in icesat_line_numbers:
            # if data over that line hasnt been found using getdataline or all_lines, find it.
            dict_entry = f'is{icesat_line_number}'
            try:
                gda_line = self.gda_lines[dict_entry]
            except KeyError:
                self.getdata_line(icesat_line_number)
                gda_line = self.gda_lines[dict_entry]
        
            for cycle in [cycle_number_from,cycle_number_till]:
                try:
                    gda_line.dropna(axis='index',subset=['h_corr']).cycle_number.value_counts().loc[cycle]
                except KeyError:
                    print(f'There is no data for line {icesat_line_number}, cycle {cycle}')
                    continue
        
            point_midchannel = np.array( LineString(gda_line.geometry.tolist()).intersection(LineString(self.channel_mid.geometry.tolist())) ) 
        
            point_midchannels.append( np.mean(point_midchannel,axis=0) )

           # This dataframe has h_corr from cycle_number_from, and dhdt. Times are from cycle_number_from
            gda_line_diff = gpd.GeoDataFrame( gda_line[gda_line.cycle_number==cycle_number_from],geometry=gda_line.geometry ).reset_index(drop=True)
            gda_line_diff.rename(columns={"h_corr": f"h_corr_cycle_{cycle_number_from}"}, inplace=True)
            gda_line_diff[f"h_corr_cycle_{cycle_number_till}"] = gda_line[gda_line.cycle_number==cycle_number_till].h_corr.reset_index(drop=True)
        
        
            gda_line_diff['dh'] =- (gda_line[gda_line.cycle_number==cycle_number_from].h_corr.to_numpy() -
                                  gda_line[gda_line.cycle_number==cycle_number_till].h_corr.to_numpy())
        
            # get the time in years between data points
            gda_line_diff['dt'] = (gda_line[gda_line.cycle_number==cycle_number_till].utc_time.to_numpy()  - 
                              gda_line[gda_line.cycle_number==cycle_number_from].utc_time.to_numpy()  ) 
            gda_line_diff.dt = gda_line_diff.dt /  np.timedelta64(1, 'Y')
        
            gda_line_diff['dhdt'] = gda_line_diff['dh'].to_numpy()/gda_line_diff['dt'].to_numpy()
            
            gda_line_diffs.append(gda_line_diff)
            
            del gda_line_diff
            
        n_plots = len(gda_line_diffs)
        
        fig, ax = plt.subplots(n_plots,1,figsize=(12,3*n_plots),dpi=500)

        i= 65
        for ax1,gda_line_diff,point_midchannel in zip(ax,gda_line_diffs,point_midchannels):
        
            color = 'tab:blue'
        
            ax1.set_ylabel(f'{chr(i)} \n elevation (m)', color=color,fontsize=10)
            i += 1
        
            ax1.plot(gda_line_diff.sort_values(by ='x' ).x,
                     gda_line_diff.sort_values(by ='x' )[f"h_corr_cycle_{cycle_number_from}"],'b-',
                    label = f"{gda_line[gda_line.cycle_number==cycle_number_from].utc_time.mean().date()}")
        
            ax1.plot(gda_line_diff.sort_values(by ='x' ).x,
                     gda_line_diff.sort_values(by ='x' )[f"h_corr_cycle_{cycle_number_till}"],'g-',
                    label = f"{gda_line[gda_line.cycle_number==cycle_number_till].utc_time.mean().date()}")
        
            ax1.tick_params(axis='x', labelcolor=color)
            
            ax1.grid()
        #     ax1.legend()
        
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        
            color = 'tab:red'
             # we already handled the x-label with ax1
            ax2.plot(gda_line_diff.x,gda_line_diff.dhdt,'ro', label='dh')
            ax2.vlines(point_midchannel[0],gda_line_diff.dhdt.min(),gda_line_diff.dhdt.max(),label='Channel midpoint')
        
            ax2.tick_params(axis='y', labelcolor=color)
            ax2.grid()
            ax2.set_ylabel('surface change, (m/a)', color=color,fontsize=10)
            
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc=0)
            
            
        ax1.set_xlabel('UTM, x',fontsize=18)
        plt.show()

                

        
        
