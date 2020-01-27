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
import time
import datetime as dt
import pandas as pd
import geopandas as gpd
import scipy as sp
from scipy import signal

from geopandas import GeoDataFrame
from shapely.geometry import Point
import fiona

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
    This function reads the metadata for all radar lines.
    
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


    

        

class radarline:
    
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
    
    
                  
    def set_filecode(self,filecode):
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
            
    def load_radar_data(self,path = "/Volumes/arc_04/FIELD_DATA/K8621920/RES/"):
            """
            load radar data from path and all directories beneath it
            
            ...could recode this to import with pandas...
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
            self.pixietimes = [pd.Timestamp( dt.datetime(2019,1,1) - dt.timedelta(days=1) + dt.timedelta(t)) for t in np.genfromtxt( self.files_paths[3] )]
            
            self.timesync = set_timesync(self.metadata.date_nzdt.strftime("%Y-%m-%d"))
            
            self.datetime = [pixietime - self.timesync for pixietime in self.pixietimes]
            self.time_str = np.array([t.strftime("%H:%M:%S %d%b%y") for t in self.datetime])
            self.pixie_time_str = np.array([t.strftime("%H:%M:%S %d%b%y") for t in self.pixietimes])
            
            self.radata = pd.DataFrame({'time':self.datetime})
            ts_func = lambda t : t.timestamp()
            self.radata['timestamp'] = self.radata.time.apply(ts_func)
            #self.radata['ch0'] = list( np.fromfile( self.files_paths[0],dtype=">f8", count=-1).reshape(-1,2500) )
            #self.radata['ch1'] = list( np.fromfile( self.files_paths[1],dtype=">f8", count=-1).reshape(-1,2500) )
            self.ch0_raw =  np.fromfile( self.files_paths[0],dtype=">f8", count=-1).reshape(-1,2500) 
            self.ch1_raw =  np.fromfile( self.files_paths[1],dtype=">f8", count=-1).reshape(-1,2500)
            self.ch0 =  self.ch0_raw
            self.ch1 =  self.ch1_raw
            
    def reset_data(self,channel=0):
            if channel==0:
                self.ch0 =  self.ch0_raw
            elif channel == 1:
                self.ch1 =  self.ch1_raw
    
    
    def load_gps_data(self,gps_path = "/Users/home/whitefar/DATA/ANT_DATA_1920/RES_GPS/2019-12-30 181325.gpx"):
            """
            output:
                - radarline.track_points, the gps file as a geodataframe
                - radarline.geodata, has (x,y) points which are recorded closest in time to the timestamps on radar instances.
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
    
            self.depth = np.real(depth2)
                
            # # if max_bottom value exists,then calculate ice_thickness
            # 
            # if(exist('max_bottom'))
            #   sz=length(max_bottom);
            #   for i=1:sz
            # 	ice_thickness(i)=depth(max(find(time<=max_bottom(i))));
            #   end
            # end

    def interpolate_gps(self):
            """
            """
            #get locations for each radar pulse
            
            x_interp_fn = sp.interpolate.interp1d(self.track_points.timestamp, self.track_points.geometry.x,kind='linear')
            x_locations = x_interp_fn(self.radata.timestamp)
            
            y_interp_fn = sp.interpolate.interp1d(self.track_points.timestamp, self.track_points.geometry.y,kind='linear')
            y_locations = y_interp_fn(self.radata.timestamp)
            
            geometry = [Point(xy) for xy in zip(x_locations, y_locations)]
            self.radata = GeoDataFrame(self.radata, geometry=geometry, crs = {'init': 'epsg:4326'} )
            line5.radata["geometry_m"] = line5.radata.geometry.to_crs(epsg=3031)
            
            #self.radata['distance_m'] is not quite right, its distance from origin not cumalative distance over track... not sure how to...
            #its hard to do cumulative distance, as most points are the same as ones after...
            self.radata['distance_from_origin'] = self.radata.to_crs(epsg=3031).distance(self.radata.to_crs(epsg=3031).geometry.iloc[0])
            
            #this is actual distance
            tmp_dfp = [Point.distance(line5.radata.geometry_m.iloc[i]) for i,Point in enumerate(line5.radata.geometry_m.iloc[1:])]
            tmp_dfp[:0] = [0]
            line5.radata['distance_from_prev'] = pd.Series(tmp_dfp) #note the 1:, equivalent to i+1
            
            line5.radata['distance_m'] = line5.radata.distance_from_prev.cumsum()
            
            
            
            
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
    

         
            
            
            
    def split_lines_a(self):
            
        
            window = 10
            line5.radata["dt"] = line5.radata["timestamp"].diff().rolling(window=window).mean()
        
            line5.radata["velocity"] = pd.Series([d/line5.radata.dt[i] for i,d in enumerate(line5.radata.distance_from_prev.rolling(window=window).mean())])
            line5.radata.velocity.fillna(0,inplace=True)
            plt.plot(line5.radata.velocity,'x')
        
            line5.index_moving = np.argwhere(line5.radata.velocity>0.01).flatten()
            
            line5.dindex_moving = np.diff(line5.index_moving)
            line5.segment_start_finish = line5.dindex_moving[line5.dindex_moving!=1]
            print(f"line has {line5.segment_start_finish.shape[0]/2} segments of moving")
        
                
    def split_lines_b(self,names=[]):
        
        if len(names) != int(line5.segment_start_finish.shape[0]/2):
            raise ValueError(f'number of names is {len(names)} but detected {line5.segment_start_finish.shape[0]/2} segments of moving')
        
        line5.section = {}
        for i,name in enumerate(names):
            start_index = line5.segment_start_finish[2*i]
            end_index = line5.segment_start_finish[2*i+1]
            
            line5.section = radarline()
            line5.section[name]['radata'] = line5.radata.iloc[range(start_index,end_index)]
            line5.section[name]['ch0'] = line5.ch0[range(start_index,end_index),:]
            line5.section[name]['ch1'] = line5.ch1[range(start_index,end_index),:]
        
    def radargram(self,channel=0,bound=0.008,title='radargram',x_axis='time'):
        """
        """
        if channel==0:
            data = self.ch0
        elif channel == 1:
            data = self.ch1
                    
        
        
        
        if x_axis == "time":
                
            ttim   = self.radata.timestamp
                    
            extent = [ttim.to_numpy()[0],ttim.to_numpy()[-1],self.depth[-1]/2,self.depth[0]]
            
            x_label = 'time, s'
            

        elif x_axis == "space":
            
            #this is a lie, putting start and end as two sides of radargram, as assumes constant speed                    
            extent = [self.radata.distance_m.iloc[0],self.radata.distance_m.iloc[-1],self.depth[-1]/2,self.depth[0]]
            
            x_label = 'distance, m'
            
        
            
        fig, ax = plt.subplots(figsize=(12,12),dpi=180)
        ax.imshow(data[:,:1250].T,vmin=-bound, vmax=bound,extent=extent,aspect='auto'  )
        ax.set_title(title)
        ax.xaxis.set_tick_params(rotation=90)
        ax.set_xlabel(x_label)
            
            
        
            
            
   

               
        
       
        
        
       
        
        
        
            
#    
#           

##1-1-2020
#linescamp = radarline("06001001502")
#linescamp.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#
#linescamp.detrend_data()
#linescamp.density_profile()
#linescamp.filter_data(High_Corner_Freq = 2.5e7)
#linescamp.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz')
#            
##31-1-2020
#lat_apres = radarline("06001000411")
#lat_apres.load_radar_data()
#lat_apres.load_gps_data()
#lat_apres.time_str



#up_chan = radarline()
#up_chan.set_filecode("06001000235")
#up_chan.load_radar_data()
#up_chan.load_gps_data()
#
##30-12-2019


#2460 11272
        
        
line5 = radarline("06364035101")
line5.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
line5.detrend_data()
line5.density_profile()
line5.filter_data(High_Corner_Freq = 2.5e7)
line5.load_gps_data()

line5.track_points.to_crs(epsg=3031)

line5.interpolate_gps()

line5.cut_stopped()

line5.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
#

#line5.datetime[0]
#
#line5.radata.time.iloc[0]
#pd.Timestamp.utcfromtimestamp(line5.radata.time.iloc[0].timestamp())
#
#
#line5.radata.timestamp.iloc[0]
#pd.Timestamp.utcfromtimestamp(line5.radata.timestamp.iloc[0])
#
#plt.plot(line5.geodata.datetime,line5.geodata.geometry.y)
#plt.xticks(rotation=90)
#
#plt.plot(line5.track_points.datetime,line5.track_points.geometry.x,'x')
#plt.xticks(rotation=90)
#plt.grid()
##
#
#df = line5.track_points
#
#
#geometry[line5.radar_to_gps_index]
#
#gpd.GeoDataFrame()
#            line5.geodata['geometry_datetime'] = line5.track_points.datetime[line5.radar_to_gps_index].to_numpy
#

#line5.density_profile()
#
#line5.reset_data()
#line5.detrend_data()
#line5.radargram(channel=0,bound=0.008,title=f'nice radargram')
#line5.reset_data()
#
#
#line5.reset_data()
#line5.detrend_data()
#line5.filter_data(High_Corner_Freq = 2.5e7)
#line5.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz')
    

#best bound is 0.01

#line5.filter_data()
#
#line5.load_gps_data()
    
    
##
###29-12-2019
#line14 = radarline("06363041031")
#line14.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#line14.detrend_data()
#line14.density_profile()
#line14.filter_data(High_Corner_Freq = 2.5e7)
#line14.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz')
#
##28-12-2019
#line11 = radarline()
#line11.set_filecode("06362023503")
#line11.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#line11.load_gps_data()


        
        #stack
#        ts_func = lambda t : t.timestamp()
#        
#        if filtdata.shape[1]>2000:
#            
#            filtdata_stacked = []
#            for sig in filtdata:
#                sig_stacked = pd.Series(sig.astype("<f8")).rolling(4,center=True,min_periods=1).mean().to_numpy()
#                filtdata_stacked.append(sig_stacked)
#            ftdata = np.array(filtdata_stacked)
#            
#            #rolling mean on the POSIX timestamp of panads series of timestamps
#            ttim = pd.Series(line5.radata.time.apply(ts_func)).rolling(4,center=True,min_periods=1).mean() 
#        else:
#            ftdata = filtdata
#            ttim   = pd.Series(line5.radata.time.apply(ts_func))
        
        #not sure i need these
        #ttim = ttim-ttim[1]
        #ttim = ttim*24*60;
        
         #x_locations = sp.interpolate.spline(line5.track_points.timestamp, line5.track_points.geometry.x, line5.radata.timestamp, order=3, kind='smoothest', conds=None)