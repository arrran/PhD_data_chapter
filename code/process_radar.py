#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:22:58 2020

@author: whitefar
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

from load_ppp import load_ppp_date

sys.path.append('/Users/home/whitefar/DATA/code')
# =============================================================================
# #check all radar files are in folder
# 
# folder = "/Volumes/arc_04/FIELD_DATA/K8621920/RES"
# 
# #do this later
# #converters = {"date_nzst": lambda t : int(round(float(s)))\
# #              "time_nzst": lambda t : } 
# 
# dtype = [('line_name',"U50"), ("date_nzst","U12"), ("started_file_nzdt","U12"), ("stopped_file_nzdt","U12"),\
#          ("waveforms","int32"), ("filecode","U12")]
# 
# meta = np.genfromtxt("/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata",skip_header=1,dtype=dtype)
# 
# for filecode in meta['filecode'][:20]:
#     if len(glob.glob(os.path.join(folder,"**",str(filecode)+"*"),recursive=True)) ==0:
#         print("files are missing for {}".format(filecode))
#     else:
#         print("file exists for {}".format(filecode))
 
# 
# =============================================================================

# =============================================================================
# Problems:

#the timesyncing on 31st is wrong
#in densprof the a, a= argwhere, im not entirely convinced it should be argwhere? maybe where?
#in auto_depthimage_t the rolling mean movmean is along spatial axes for signal but along time axis for time

# I dont think it works getting closest gps point, see line5.geodata.datetime.iloc[0] is diff to line5.radata.time.iloc[0]
# =============================================================================


# 
def metadata_func(fc):
    """
    This function reads the metadata for all radar lines for use with the radarsurvey data class
    
    input: a filecode (fc) 
    output: a callable dataframe with entries [line_name date_nzdt started_file_nzdt stopped_file_nzdt waveforms filecode]
    
    """
        
    metadata_path = "/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata"
   
    #not sure the convertors work, it seems to use pandas time format, which seem compatible with datetime
    converters1 = {"date_nzdt" :lambda t : dt.datetime.strptime(t,"%Y-%m-%d"),\
                "started_file_nzdt": lambda t : dt.datetime.strptime(t,"%H:%M"),\
                "stopped_file_nzdt" :lambda t : dt.datetime.strptime(t,"%H:%M")}
   
   
    metadata = pd.read_csv(metadata_path,delimiter=' ',converters=converters1)
   
    #date was in a separate column to time, add the date to the time series to make datetimes
    metadata.started_file_nzdt = [dt.datetime.combine(metadata.date_nzdt[i].date(),metadata.started_file_nzdt[i].time()) for i,_ in enumerate(metadata.started_file_nzdt)]
    metadata.stopped_file_nzdt = [dt.datetime.combine(metadata.date_nzdt[i].date(),metadata.stopped_file_nzdt[i].time()) for i,_ in enumerate(metadata.stopped_file_nzdt)]
   
    #list of filecodes
    filecodes = np.loadtxt("/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata",dtype=str,skiprows=1,usecols=5)
   
    #dictionary which returns the row in metadata given filecode
    filecode2metarow = {filecode:row for filecode,row in zip(filecodes,np.arange(0,len(filecodes)))}
   
    return metadata.loc[filecode2metarow[fc]]

def get_metadata(fc):
    """
    This function reads the metadata for all radar lines. outputs just a pandas dataframe with metadata
    
    input: a filecode (fc) 
    output: a callable dataframe with entries [line_name date_nzdt started_file_nzdt stopped_file_nzdt waveforms filecode]
    
    """
        
    metadata_path = "/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata"
   
    #not sure the convertors work, it seems to use pandas time format, which seem compatible with datetime
    converters1 = {"date_nzdt" :lambda t : dt.datetime.strptime(t,"%Y-%m-%d"),\
                "started_file_nzdt": lambda t : dt.datetime.strptime(t,"%H:%M"),\
                "stopped_file_nzdt" :lambda t : dt.datetime.strptime(t,"%H:%M")}
   
   
    metadata = pd.read_csv(metadata_path,delimiter=' ',converters=converters1)
   
    #date was in a separate column to time, add the date to the time series to make datetimes
    metadata.started_file_nzdt = [dt.datetime.combine(metadata.date_nzdt[i].date(),metadata.started_file_nzdt[i].time()) for i,_ in enumerate(metadata.started_file_nzdt)]
    metadata.stopped_file_nzdt = [dt.datetime.combine(metadata.date_nzdt[i].date(),metadata.stopped_file_nzdt[i].time()) for i,_ in enumerate(metadata.stopped_file_nzdt)]
   
    #list of filecodes
    filecodes = np.loadtxt("/Volumes/arc_04/FIELD_DATA/K8621920/RES/radar_metadata",dtype=str,skiprows=1,usecols=5)
   
    #dictionary which returns the row in metadata given filecode
    filecode2metarow = {filecode:row for filecode,row in zip(filecodes,np.arange(0,len(filecodes)))}
   
    return metadata

def set_timesync(date_in):
    """
    This reads the timesync data
    
    input: a date string, format "%Y-%m-%d"
    output: a time delta, the time difference between pixie_time and utc_time written on the file "time_sync"
    """
    converters = {"exact_nzdt" :lambda t : dt.datetime.strptime(str(t),"%d-%m-%YT%H:%M:%S"),\
                "pixie_time": lambda t : dt.datetime.strptime(str(t),"%d-%m-%YT%H:%M:%S"),\
                "exact_utc_time" :lambda t : dt.datetime.strptime(str(t),"%d-%m-%YT%H:%M:%S")}
    
        
    timesync = pd.read_csv("/Volumes/arc_04/FIELD_DATA/K8621920/RES/time_sync",delimiter=' ',converters=converters)
    #get a time delta
    timesync["dt"] = timesync.pixie_time.array - timesync.exact_utc_time.array
    #make a dictionary which returns the time delta given a date
        
    timesync_dict = {date:dt for date,dt in zip([D.date().strftime("%Y-%m-%d") for D in timesync.exact_nzdt], timesync.dt)}

    
    return timesync_dict[date_in]


def density_profile(separation_distance = 58.37):
            """
            """
          
            rho_o=292  # surface snow density
            rho_f=917  # ice density
            s = separation_distance
            z=np.arange(0,8000+1)
            
            # ridge BC r^2(fit)=0.8380 
            #rhoRBC=rho_o+(rho_f-rho_o)*(1-exp(-z/26.73));
            rhoRBC=rho_o+(rho_f-rho_o)*(1-np.exp(-0.0386*z))
            
            volfracRBC=rhoRBC/rho_f  # volfrac is the density relative to ice
            
            # determine effective dielectric constant for a mixture of two materials
            # using diel_m.m where the two materials are air (dielectic constant of 1)
            
            diel_m = lambda VolFrac2, e_1, e_2 : (VolFrac2*( e_2**(1/3) - e_1**(1/3) )+e_1**(1/3))**3
            
            # diel_mix_r.m
            # yields the effective dielectric constant for a mixture of e_1 and e_2
            # where VolFrac2 is the volume fraction of e_2
            #  According to: Looyenga's Equation
            
            # Rock: e_1=8
            # Water: e_2=88
            # Typical: VolFrac2=0.25
            # and ice (dielectric constant of 3.12).
    
            e_effRBC = diel_m(volfracRBC, 1,3.12)
            
            velRBC=300/np.sqrt(e_effRBC) #velocity of radio waves in air is 300 m/mus 
            
            szz=len(z)
            int_avg_velRBC_c=np.cumsum(velRBC)
            int_avg_velRBC=int_avg_velRBC_c/np.arange(1,szz+1)
            
            #int_avg_vel(i) is the integrated average velocity to depth i
            # the tt vs depth curve is calculated by multiplying z./int_avg_vel
            # this is the travel time curve as a function of z
    
            ttimeRBC=z/int_avg_velRBC
            
                        
            # to use this to find an actual depth, calculate the two-way travel time t to
            # an object.  The depth of an object is then related to that and the
            # separation distance s.  e.g. if t=5uSsecond and s=90m
            # a=max(find(ttime<=(t/2+s/300)));  depth=sqrt(a^2-(s/2)^2);
            #
            # for t=5.56uS, s=100m, depth ~= 506.3m
            # not accounting for geometry or s, one would otherwise interpret a 5.56uS 
            # reflection as being 5.56*86.5=480m;
                        
            depthRBC = []
            for i,ttime in enumerate(ttimeRBC):
                a = np.max( np.argwhere( ttimeRBC <= ttime/2 ) ).astype(complex)
                depthRBC.append( np.sqrt( a**2-((s+0j)/2)**2 ) )
                
            # note that the first few values of depth1 will be complex.  This
            # is due to the fact that there are travel times for the reflected
            # wave that cannot exist.  This is because the reflected wave travels 
            # slower than the direct wave.  Even for an infinitely shallow reflection,
            # the travel time for the reflected wave must at least slightly lag the 
            # direct wave.  This lag
            # depends on the separation distance and velocity difference, and is
            # the duration of time for which travel times cannot occur.
            
            # to convert to measured travel time with T=0 being the arrival of the direct
            # wave (not the actual travel time which is what is used so far...  
            # subtract the time it takes for the direct wave to travel.
                        
            ttimeRBC_2=ttimeRBC-s/300
            depthRBC_r=np.real(depthRBC)
            
            Xinc = 10e-9
            time = np.arange(-Xinc*125,Xinc*(2500-125),Xinc)
            
            if np.log10(np.amax(time))<-3:
                time1=time*1e6
            else:
                time1=time
              
            # eliminate negative numbers;
            #  i=find(time1<0);
            #  time1(i)=0;
            
            #this makes variable 'depth', which reduces depthRBC to 2500 elements, where every element where the travel time is less than the time1 is defaulted to depthRBC_r[0]
            depth_tmp = []
            for i, t in enumerate(time1):
                if np.argwhere(ttimeRBC_2<t).size == 0:
                    j = 0
                else:
                    j = np.max( np.argwhere(ttimeRBC_2<t) )                
                depth_tmp.append( depthRBC_r[j] )
                
            #    if(real(depth(i)) <= 0)
            #	  depth(i)=time1(i)*86;
            #	end
            
            #this wee loop replaces anything with depth0
            depth2 = depth_tmp
            if np.where(depth_tmp==0)[0].size != 0:
                k = np.amax( np.where(depth_tmp==0) )
                
                for i  in np.arange(1,k+1):
                    depth2[i]=-(k-i)*Xinc[0]*84e6
            
            global depth
            
            depth = np.real(depth2)
                
            # # if max_bottom value exists,then calculate ice_thickness
            # 
            # if(exist('max_bottom'))
            #   sz=length(max_bottom);
            #   for i=1:sz
            # 	ice_thickness(i)=depth(max(find(time<=max_bottom(i))));
            #   end
            # end

density_profile()            
         

    

        

class radarsurvey:
    """
    One period of continuously recording radar
    """
    
    def __init__(self,filecode):
            """
            filecode is the code given to the radar line eg 06348013011
            """
            self.filecode = filecode
            self.ch1_filename = filecode + "ch1"
            self.ch0_filename = filecode + "ch0"
            self.info_filename = filecode + "info.txt"
            self.time_filename = filecode + "time.txt"
            self.filenames = [self.ch0_filename,self.ch1_filename,self.info_filename, self.time_filename]
            self.metadata = metadata_func(filecode)
            self.timesync = set_timesync(self.metadata.date_nzdt.strftime("%Y-%m-%d"))
            
            
            
                
    def load_radar_data(self,path = "/Volumes/arc_04/FIELD_DATA/K8621920/RES/"):
            """
            load radar data from path and all directories beneath it
            
            ...could recode this to import with pandas...
            
            This also stacks temporaly so that there is only one trace per timestamp
            """
            
            if len(glob.glob(os.path.join(path,"**",self.filenames[0]),recursive=True)) == 0:
                raise ValueError("Can't find files with that filecode")
            
            self.files_paths = [glob.glob(os.path.join(path,"**",filename),recursive=True)[0] for filename in  self.filenames]
            
            print("loading files from:")
            for file_path in self.files_paths:
                print(file_path)
            
    
            
            #two binary files in big endian
            #ch0 = np.fromfile( files_paths[0],dtype=">f8", count=-1).reshape(-1,2500)
            #ch1 = np.fromfile( files_paths[1],dtype=">f8", count=-1).reshape(-1,2500)
            self.info = np.genfromtxt( self.files_paths[2] ,delimiter=',') # two text files
            
            #for consitancy, all time data in pd timstamps
            pixietimes_raw = [pd.Timestamp( dt.datetime(2019,1,1) - dt.timedelta(days=1) + dt.timedelta(t)) for t in np.genfromtxt( self.files_paths[3] )]
            
            datetime_raw = [pixietime - self.timesync for pixietime in pixietimes_raw]
            self.pixietimes_raw = pixietimes_raw
            
            
            datetime_df_raw = pd.DataFrame({'datetime':datetime_raw})
            ts_func = lambda t : pd.Timestamp.timestamp(t)
            timestamp_raw = datetime_df_raw.datetime.apply(ts_func)
            
            #self.radata['ch0'] = list( np.fromfile( self.files_paths[0],dtype=">f8", count=-1).reshape(-1,2500) )
            #self.radata['ch1'] = list( np.fromfile( self.files_paths[1],dtype=">f8", count=-1).reshape(-1,2500) )
            
            ch0_raw =  np.fromfile( self.files_paths[0],dtype=">f8", count=-1).reshape(-1,2500) 
            ch1_raw =  np.fromfile( self.files_paths[1],dtype=">f8", count=-1).reshape(-1,2500)
            
            #deal with the change in radar signals recorded per unit time
            splittimes_index = np.hstack([np.argwhere(timestamp_raw.diff().to_numpy() !=0 ).flatten(),-1])
            
            stacked_ch0 = np.zeros([len(splittimes_index)-1,ch0_raw.shape[1]])
            for i, index in enumerate(splittimes_index[:-1]):
                stacked_ch0[i,:] = np.mean( ch0_raw[ index:splittimes_index[i+1],: ],axis=0 )    
            self.ch0 = stacked_ch0
            
            stacked_ch1 = np.zeros([len(splittimes_index)-1,ch1_raw.shape[1]])
            for i, index in enumerate(splittimes_index[:-1]):
                stacked_ch1[i,:] = np.mean( ch1_raw[ index:splittimes_index[i+1],: ],axis=0 )    
            self.ch1 = stacked_ch1
            
            self.radata = pd.DataFrame({'timestamp':timestamp_raw.to_numpy()[splittimes_index]})
            
            dt_func = lambda t : pd.Timestamp.utcfromtimestamp(t)
            self.radata['datetime'] =  self.radata.timestamp.apply(dt_func)
            
            self.radata = self.radata.reset_index()
            
#            self.time_offset_start =  self.metadata.started_file_nzdt - self.radata.datetime.iloc[0]
#            self.time_offset_stopped = self.metadata.stopped_file_nzdt - self.radata.datetime.iloc[-1]
       
       
        

            
    def reset_data(self,channel=0):
            if channel==0:
                self.ch0 =  self.ch0_raw
            elif channel == 1:
                self.ch1 =  self.ch1_raw
    
    
    def load_garmin_data(self,gps_path = "/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RES_GPS/2019-12-30 181325.gpx"):
            """
            output:
                - radarsurvey.track_points, the gps file as a geodataframe
                #- radarsurvey.geodata, has (x,y) points which are recorded closest in time to the timestamps on radar instances.
            """
            ts_func = lambda t : t.timestamp()
            #load the track as geodataframe
            self.track_points = gpd.read_file(gps_path,layer='track_points')
            
            #convert to datetime object
            self.track_points['datetime'] = np.array([dt.datetime.strptime(t,"%Y-%m-%dT%H:%M:%S") for t in self.track_points.time])
            self.track_points['timestamp'] = self.track_points.datetime.apply(ts_func)
            
            #self.radar_to_gps_index = [np.argwhere( abs(self.track_points.datetime - t)==abs(self.track_points.datetime - t).min() )[0][0] for t in self.datetime]
            
            #geometry = line5.track_points.geometry.iloc[self.radar_to_gps_index]
            #self.geodata = GeoDataFrame(geometry, geometry=geometry)
            #self.geodata['datetime'] = self.track_points.datetime.iloc[self.radar_to_gps_index]
            
            #self.geodata['timestamp'] = self.geodata.datetime.apply(ts_func)
            #self.geodata.reset_index(drop=True,inplace=True)
            
    def extra_position_info(self):
            """
            adds a smoothed velocity and acceleration
            
            """
            window = 5
            self.track_points["dt"] = self.track_points["timestamp"].diff().rolling(window=window).mean()
            if np.argwhere(self.track_points.dt==0).shape[0] != 0:
                raise ValueError("problem with zero valued dt, need bigger window")
            
            self.track_points["geometry_m"] = self.track_points.geometry.to_crs(epsg=3031)
            tmp_dfp = [Point.distance(self.track_points.geometry_m.iloc[i]) for i,Point in enumerate(self.track_points.geometry_m.iloc[1:])]
            tmp_dfp[:0] = [0]
            self.track_points['distance_from_prev'] = pd.Series(tmp_dfp)
            self.track_points['distance_along_line'] = self.track_points.distance_from_prev.cumsum()
            
            
            self.track_points["velocity"] = pd.Series([d/self.track_points.dt[i] for i,d in enumerate(self.track_points.distance_from_prev.rolling(window=window).mean())])
            self.track_points["acc"] = pd.Series([d/self.track_points.dt[i] for i,d in enumerate(self.track_points.velocity)]).rolling(window=window).mean()
     
    def load_gnss_data(self):
            """
            loads ppp gnss data over one or two days
            """
            
            start_date_utm = self.radata.datetime.iloc[0].strftime('%Y-%m-%d')
            end_date_utm = self.radata.datetime.iloc[-1].strftime('%Y-%m-%d')
            
            if start_date_utm == end_date_utm:
                self.track_points = load_ppp_date([start_date_utm])
            if start_date_utm != end_date_utm:
                self.track_points = load_ppp_date([start_date_utm,end_date_utm])
           
    
    def interpolate_gnss(self):
            """
            INPUT:  self.timestamp - timestamp taken from the radar
                    self.track_points.geometry - position data from the gnss
                    
            OUTPUT: self.radata.geometry, self.radata.geometry, -  position from gnss interpolated and picked at points where the radar pipped
                    self.distance_from_origin, distance_from_prev (dx), distance_m
            
            Must first run radarsurvey.load_gnss_data
            
            NB if you try run this twice the to_crs method does not work!!!
            """
            
            #get locations for each radar pulse
            
            x_interp_fn = sp.interpolate.interp1d(self.track_points.timestamp, self.track_points.geometry.x,kind='linear')
            x_locations = x_interp_fn(self.radata.timestamp)
            
            y_interp_fn = sp.interpolate.interp1d(self.track_points.timestamp, self.track_points.geometry.y,kind='linear')
            y_locations = y_interp_fn(self.radata.timestamp)
            
            geometry = [Point(xy) for xy in zip(x_locations, y_locations)]
            self.radata = GeoDataFrame(self.radata, geometry=geometry, crs = {'init': 'epsg:3031'} )
            self.radata["geometry_m"] = self.radata.geometry.to_crs(epsg=3031)
            
            #self.radata['distance_m'] is not quite right, its distance from origin not cumalative distance over track... not sure how to...
            #its hard to do cumulative distance, as most points are the same as ones after...
            self.radata['distance_from_origin'] = self.radata.to_crs(epsg=3031).distance(self.radata.to_crs(epsg=3031).geometry.iloc[0])
            
            #this is actual distance
            tmp_dfp = [Point.distance(self.radata.geometry_m.iloc[i]) for i,Point in enumerate(self.radata.geometry_m.iloc[1:])]
            tmp_dfp[:0] = [0]
            self.radata['distance_from_prev'] = pd.Series(tmp_dfp) #note the 1:, equivalent to i+1
            self.radata['distance_along_line'] = self.radata.distance_from_prev.cumsum()
            
            self.radata['distance_m'] = self.radata.distance_from_prev.cumsum()
            
            
            
            
            #        line5.radata.iloc[40:80].distance(line5.radata.geometry.iloc[30:80])
    #        
    #        plt.plot(line5.radata.iloc[65:75].geometry.x)
    #        
    #        line5.radata.iloc[69].geometry.distance(line5.radata.geometry.iloc[71])
    #        
    #        line5.radata.iloc[69:70].geometry.distance(line5.radata.geometry.iloc[70:71])
    
            
    
            
            #line5.radata['distance_m'] = line5.radata.to_crs({'epsg:3031'}).distance(line5.radata.to_crs({'init': 'epsg:3031'}).geometry.iloc[0])
            
            
            
            
    #        line5.radata['distance_m'] = line5.radata.to_crs('epsg:3031').iloc[1:].distance(line5.radata.to_crs('epsg:3031').geometry.iloc[:-1])
    #        
    #    
    #        line5.radata['distance_m'] = line5.radata.to_crs('epsg:3031').geometry.iloc[1:].distance(line5.radata.to_crs('epsg:3031').geometry.iloc[:-1]).to_numpy().cumsum()
    #        
    #        line5.radata.to_crs('epsg:3031').geometry.iloc[66].distance(line5.radata.to_crs('epsg:3031').geometry.iloc[67])
    #        
    #        distance_from_prev = [Point.distance(line5.radata.to_crs('epsg:3031').geometry[i]) for i,Point in enumerate(line5.radata.to_crs('epsg:3031').geometry[1:])] #note the 1:, equivalent to i+1
    
    
    def detrend_data(self,channel=0):
            """
            """
            
            if channel==0:
                self.ch0 =  signal.detrend(self.ch0, axis=1, type='constant', bp=0) #centres the signal about zero
            elif channel == 1:
                self.ch1 =  signal.detrend(self.ch1, axis=1, type='constant', bp=0) #centres the signal about zero
            
    
    def radargram(self,channel=0,bound=0.008,title='radargram',x_axis='time'):
        """
        """
        if channel==0:
            data = self.ch0
        elif channel == 1:
            data = self.ch1
                    
        
        
        
        if x_axis == "time":
                
            ttim   = self.radata.timestamp
                    
            extent = [ttim.to_numpy()[0],ttim.to_numpy()[-1],depth[-1]/2,depth[0]]
            
            x_label = 'time, s'
            

        elif x_axis == "space":
            
            #this is a lie, putting start and end as two sides of radargram, as assumes constant speed                    
            extent = [self.radata.distance_m.iloc[0],self.radata.distance_m.iloc[-1],depth[-1]/2,depth[0]]
            
            x_label = 'distance, m'
            
        
            
        fig, ax = plt.subplots(figsize=(12,12),dpi=180)
        ax.imshow(data[:,:1250].T,vmin=-bound, vmax=bound,extent=extent,aspect='auto'  )
        ax.set_title(title)
        ax.xaxis.set_tick_params(rotation=90)
        ax.set_xlabel(x_label)
            
               
            
    def split_lines_choose(self,threshold_type = 'acc',moving_threshold=1,window = 2,plot_radargram = True):
        """
        """
        
        self.radata["dt"] = self.radata["timestamp"].diff().rolling(window=window).mean()
        if np.argwhere(self.radata.dt==0).shape[0] != 0:
            raise ValueError("problem with zero valued dt, need bigger window")
    
        self.radata["velocity"] = pd.Series([d/self.radata.dt[i] for i,d in enumerate(self.radata.distance_from_prev.rolling(window=window).mean())])
        self.radata["acc"] = pd.Series([d/self.radata.dt[i] for i,d in enumerate(self.radata.velocity)]).rolling(window=window).mean()
                    
        self.radata.acc.fillna(0,inplace=True)
        
        plt.figure()
        plt.plot(self.radata.velocity)
        plt.title('velocity profile')
        plt.hlines(moving_threshold,0,len(self.radata.velocity))
        
        plt.figure()
        plt.plot(self.radata.acc)
        plt.title('acc profile')
        plt.hlines(moving_threshold,0,len(self.radata.velocity))
        
        if threshold_type == 'acc':
            self.index_moving = np.argwhere(self.radata.acc>moving_threshold).flatten()
        elif threshold_type == 'velocity':
            self.index_moving = np.argwhere(self.radata.velocity>moving_threshold).flatten()
            
                
        
        self.dindex_moving = np.diff(self.index_moving)
        self.segment_splits = np.argwhere(self.dindex_moving>1).flatten()+1
        self.index_moving_segments = np.split(self.index_moving,self.segment_splits)
        self.number_of_segments = len(self.index_moving_segments)
        print(f"line has {self.number_of_segments} segments of moving, where "+threshold_type+" < " + str(moving_threshold) )
        
        
        if plot_radargram == True:
            
            
            data = signal.detrend(self.ch0, axis=1, type='constant', bp=0)
        
            #this is a lie, putting start and end as two sides of radargram, as assumes constant speed                    
            #extent = [self.radata.distance_m.iloc[0],self.radata.distance_m.iloc[-1],depth[-1]/2,depth[0]]
            
            #x_label = 'distance, m'
                           
            bound=0.008
            
            fig, ax = plt.subplots(figsize=(12,12),dpi=180)
            ax.imshow(data[:,:1250].T,vmin=-bound, vmax=bound,aspect='auto'  )
            ax.xaxis.set_tick_params(rotation=90)
            #ax.set_xlabel(x_label)
            for segment in self.index_moving_segments:
                ax.axvline(segment[0])
                ax.axvline(segment[-1])
                
            
#    def split_lines_show(self,names):
#        plt.plot(self.radata.velocity,'x')
#        
#        for i,segment_indicies in enumerate(self.index_moving_segments):
#            plt.text(5,segment_indicies[0],names[i])
                
            
    def split_lines_plot(self,names):
            
            for i,segment_indicies in enumerate(self.index_moving_segments):
                self.radata.iloc[segment_indicies].plot()
                plt.title(names[i])
                
            
                
    def split_lines_output(self):
            """
            
            e.g. self,turn,line3,turn2,line4 = radarsurvey.split_lines_outout()
            """
            
            sections = []
            for segment_indicies in self.index_moving_segments:
                
                section = {'radata': self.radata.iloc[segment_indicies],
                           'ch0':  self.ch0[segment_indicies,:], 
                           'ch1': self.ch1[segment_indicies,:],
                           'info': self.info,
                           }
                
                sections.append(section)
                
#            print(f'returning {len(sections)} sections as list. ' )
                
            return sections
        
    def refine_timesync(self,Dt):
        
        print('positive dt moves lines left')
        self.timesync =  set_timesync(self.metadata.date_nzdt.strftime("%Y-%m-%d")) + pd.Timedelta(Dt)
        
        refine_timesync_func = lambda t : t + pd.Timedelta(Dt)
        
        self.radata.datetime =  self.radata.datetime.apply(refine_timesync_func)
        
        ts_func = lambda t : pd.Timestamp.timestamp(t)
        
        self.radata.timestamp = self.radata.datetime.apply(ts_func)

        self.interpolate_gnss()
        
        
    

# =============================================================================
            
class radarline:
    
    def __init__(self,input_dictionary):
        """
        The input dictionary from split_lines_output(self)
        """
        
        self.radata = input_dictionary["radata"]
        self.ch0_raw = input_dictionary["ch0"]
        self.ch1_raw = input_dictionary["ch1"]
        self.ch0 = input_dictionary["ch0"]
        self.ch1 = input_dictionary["ch1"]
        self.info =  input_dictionary["info"]
            
    def stack_spatially(self,stack_distance=5):
        """
        stack distance in m
        """
        
        bins = np.arange(-5,self.radata.distance_along_line.iloc[-1]+10,stack_distance)
        
        self.radata['distance_bins'] = pd.cut(self.radata.distance_along_line, bins ,labels=(bins[:-1]+5)   )     #then average each bin see https://stackoverflow.com/questions/45273731/binning-column-with-python-pandas#45273750
        
                #now stack over the distance bins
        splitdistance_index = np.hstack([np.argwhere(self.radata.distance_bins.diff().to_numpy() !=0 ).flatten(),-1])
        
        stacked_ch0 = np.zeros([len(splitdistance_index)-1,self.ch0.shape[1]])
        for i, index in enumerate(splitdistance_index[:-1]):
            stacked_ch0[i,:] = np.mean( self.ch0[ index:splitdistance_index[i+1],: ],axis=0 )    
        self.ch0 = stacked_ch0
        
        stacked_ch1 = np.zeros([len(splitdistance_index)-1,self.ch1.shape[1]])
        for i, index in enumerate(splitdistance_index[:-1]):
            stacked_ch1[i,:] = np.mean( self.ch1[ index:splitdistance_index[i+1],: ],axis=0 )    
        self.ch1 = stacked_ch1
        
        self.radata = self.radata.iloc[splitdistance_index]
        self.radata = self.radata.reset_index()
        
        
        
            
    def stack_data(self,channel=0,stack=30):   
        """
        """
        
        if channel==0:
            data = self.ch0
        elif channel == 1:
            data = self.ch1
        else:
            print('Channel 0 or 1 not chosen')
        
        #this is stacking over space, not time
        filtdata_stacked = []
        for sig in data:
            sig_stacked = pd.Series(sig.astype("<f8")).rolling(stack,center=True,min_periods=1).mean().to_numpy()
            filtdata_stacked.append(sig_stacked)
            
        if channel==0:
            self.ch0 = np.array(filtdata_stacked)
        elif channel == 1:
            self.ch1 = np.array(filtdata_stacked)
         
            
    def detrend_data(self,channel=0):
            """
            """
            
            if channel==0:
                self.ch0 =  signal.detrend(self.ch0, axis=1, type='constant', bp=0) #centres the signal about zero
            elif channel == 1:
                self.ch1 =  signal.detrend(self.ch1, axis=1, type='constant', bp=0) #centres the signal about zero
            
    
    
    
    def filter_data(self,channel=0,High_Corner_Freq = 3e6):
            """
            """
            
            
            #        
            Xinc = self.info[0] # Xinc is in units: seconds/sample
            Yoffset = self.info[1]
            Yinc = self.info[2]
            Number_of_data = self.info[3]
            
            # low pass filter - 4th order butterworth design
            # first determine the cut-off frequency - expressed as a 
            #	fraction of the Nyquist frequency (half the sample freq).
            # 	all of this is in Hz
            
            #High_Corner_Freq = 3e6 #0.5e6       # high cut-off freq in hz
            print(f'Lowpassing below {High_Corner_Freq/1e6} MHz')
            Sample_Freq = int(1/Xinc)
            Nyquist_Freq = int(Sample_Freq/2)
            
            Corner_Freq = High_Corner_Freq/Nyquist_Freq
            
            
            
            # calculate the filter polynomials for 4th order butterworth lowpass
            b, a = signal.butter(4, Corner_Freq, btype='low', analog=False, output='ba')       
            
            # now, sweep through the data and filter each waveform using filtfilt
            
            if channel==0:
                self.ch0 = signal.filtfilt(b, a, self.ch0, axis=1,
                                                       padtype='odd', padlen=None,
                                                       method='pad', irlen=None)
            elif channel == 1:
                self.ch1 = signal.filtfilt(b, a, self.ch1, axis=1,
                                                       padtype='odd', padlen=None,
                                                       method='pad', irlen=None)
            
            
    def density_profile(self,separation_distance = 58.37):
            """
            """
          
            rho_o=292  # surface snow density
            rho_f=917  # ice density
            s = separation_distance
            z=np.arange(0,8000+1)
            
            # ridge BC r^2(fit)=0.8380 
            #rhoRBC=rho_o+(rho_f-rho_o)*(1-exp(-z/26.73));
            rhoRBC=rho_o+(rho_f-rho_o)*(1-np.exp(-0.0386*z))
            
            volfracRBC=rhoRBC/rho_f  # volfrac is the density relative to ice
            
            # determine effective dielectric constant for a mixture of two materials
            # using diel_m.m where the two materials are air (dielectic constant of 1)
            
            diel_m = lambda VolFrac2, e_1, e_2 : (VolFrac2*( e_2**(1/3) - e_1**(1/3) )+e_1**(1/3))**3
            
            # diel_mix_r.m
            # yields the effective dielectric constant for a mixture of e_1 and e_2
            # where VolFrac2 is the volume fraction of e_2
            #  According to: Looyenga's Equation
            
            # Rock: e_1=8
            # Water: e_2=88
            # Typical: VolFrac2=0.25
            # and ice (dielectric constant of 3.12).
    
            e_effRBC = diel_m(volfracRBC, 1,3.12)
            
            velRBC=300/np.sqrt(e_effRBC) #velocity of radio waves in air is 300 m/mus 
            
            szz=len(z)
            int_avg_velRBC_c=np.cumsum(velRBC)
            int_avg_velRBC=int_avg_velRBC_c/np.arange(1,szz+1)
            
            #int_avg_vel(i) is the integrated average velocity to depth i
            # the tt vs depth curve is calculated by multiplying z./int_avg_vel
            # this is the travel time curve as a function of z
    
            ttimeRBC=z/int_avg_velRBC
            
                        
            # to use this to find an actual depth, calculate the two-way travel time t to
            # an object.  The depth of an object is then related to that and the
            # separation distance s.  e.g. if t=5uSsecond and s=90m
            # a=max(find(ttime<=(t/2+s/300)));  depth=sqrt(a^2-(s/2)^2);
            #
            # for t=5.56uS, s=100m, depth ~= 506.3m
            # not accounting for geometry or s, one would otherwise interpret a 5.56uS 
            # reflection as being 5.56*86.5=480m;
                        
            depthRBC = []
            for i,ttime in enumerate(ttimeRBC):
                a = np.max( np.argwhere( ttimeRBC <= ttime/2 ) ).astype(complex)
                depthRBC.append( np.sqrt( a**2-((s+0j)/2)**2 ) )
                
            # note that the first few values of depth1 will be complex.  This
            # is due to the fact that there are travel times for the reflected
            # wave that cannot exist.  This is because the reflected wave travels 
            # slower than the direct wave.  Even for an infinitely shallow reflection,
            # the travel time for the reflected wave must at least slightly lag the 
            # direct wave.  This lag
            # depends on the separation distance and velocity difference, and is
            # the duration of time for which travel times cannot occur.
            
            # to convert to measured travel time with T=0 being the arrival of the direct
            # wave (not the actual travel time which is what is used so far...  
            # subtract the time it takes for the direct wave to travel.
                        
            ttimeRBC_2=ttimeRBC-s/300
            depthRBC_r=np.real(depthRBC)
            
            Xinc = self.info[0]
            time = np.arange(-Xinc*125,Xinc*(2500-125),Xinc)
            
            if np.log10(np.amax(time))<-3:
                time1=time*1e6
            else:
                time1=time
              
            # eliminate negative numbers;
            #  i=find(time1<0);
            #  time1(i)=0;
            
            #this makes variable 'depth', which reduces depthRBC to 2500 elements, where every element where the travel time is less than the time1 is defaulted to depthRBC_r[0]
            depth = []
            for i, t in enumerate(time1):
                if np.argwhere(ttimeRBC_2<t).size == 0:
                    j = 0
                else:
                    j = np.max( np.argwhere(ttimeRBC_2<t) )                
                depth.append( depthRBC_r[j] )
                
            #    if(real(depth(i)) <= 0)
            #	  depth(i)=time1(i)*86;
            #	end
            
            #this wee loop replaces anything with depth0
            depth2 = depth
            if np.where(depth==0)[0].size != 0:
                k = np.amax( np.where(depth==0) )
                
                for i  in np.arange(1,k+1):
                    depth2[i]=-(k-i)*Xinc[0]*84e6
    
            depth = np.real(depth2)
                
            # # if max_bottom value exists,then calculate ice_thickness
            # 
            # if(exist('max_bottom'))
            #   sz=length(max_bottom);
            #   for i=1:sz
            # 	ice_thickness(i)=depth(max(find(time<=max_bottom(i))));
            #   end
            # end


         
            
 
        
    def radargram(self,channel=0,bound=0.008,title='radargram',x_axis='time'):
        """
        """
        if channel==0:
            data = self.ch0
        elif channel == 1:
            data = self.ch1
                    
        
        
        
        if x_axis == "time":
                
            ttim   = self.radata.timestamp
                    
            extent = [ttim.to_numpy()[0],ttim.to_numpy()[-1],depth[-1]/2,depth[0]]
            
            x_label = 'time, s'
            

        elif x_axis == "space":
            
            #this is a lie, putting start and end as two sides of radargram, as assumes constant speed                    
            extent = [self.radata.distance_m.iloc[0],self.radata.distance_m.iloc[-1],depth[-1]/2,depth[0]]
            
            x_label = 'distance, m'
            
        
            
        fig, ax = plt.subplots(figsize=(12,12),dpi=180)
        ax.imshow(data[:,:1250].T,vmin=-bound, vmax=bound,extent=extent,aspect='auto'  )
        ax.set_title(title)
        ax.xaxis.set_tick_params(rotation=90)
        ax.set_xlabel(x_label)
