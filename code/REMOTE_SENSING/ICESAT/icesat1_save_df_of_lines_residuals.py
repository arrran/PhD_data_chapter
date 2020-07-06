#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:00:37 2020

@author: whitefar


Save a df of all attributes x y z zp etc so that i can load and plot them

The coordinates for zp, are the average x and y of the bins. Draw a line through this I think.

NB extrapolation is done by x for 0211 and y for 0099, v important, must redo REMA projection if changing this
"""

# You can calculate the residuals from the track*_all_dzdt.mat files in /Users/home/horganhu/ICESAT_LINK/TAMATA_ICESAT/GLA12_633_DZDT/

# track0099_all_dzdt.mat
# dzdt.track0099.bin2963.grad,
# dzdt.track0099.bin2963.dzdt,
# dzdt.track0099.bin2963.x,
# dzdt.track0099.bin2963.y,
# dzdt.track0099.bin2963.z,
# dzdt.track0099.bin2963.retide,
# dzdt.track0099.bin2963.sdatedays,
# dzdt.track0099.bin2963.elev_wgs84_tpxtide,
# dzdt.track0099.bin2963.lat,
# dzdt.track0099.bin2963.lon,
# dzdt.track0099.bin2963.ret

# zp = [x-mean(x); y-mean(y)]'*grad + (t-t0)'*(bindzdt) + mean(z); % Note here we use z here (not elev)

# zp=zp';

# res = z - zp;     % transient changes (footprint elevation changes corrected for
#                 % slope and secular change. Note here we use z here (not elev)

# dz=res + ((t - t0).*bindzdt);   % Elevation change corrected for each segment. This is
#                 % corrected for gradient. Note, this is per day.

import h5py
import xarray as xr
import glob
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point, LineString
import geopandas as gpd
from shapely.ops import nearest_points
import numpy as np
from tqdm import tqdm
from scipy import interpolate
import datetime as dt
from scipy.signal import savgol_filter
from shapely.affinity import scale

## I ran this on ipython terminal
def icesat1_alldzdt_todataframe(track):
    """
    this function will import the .mat and arrange into a dataframe
    track = eg track
    
    da is a Dataframe restricted to the channel Area
    """
    
    
    
    path = f'/Users/home/horganhu/ICESAT_LINK/TAMATA_ICESAT/GLA12_633_DZDT/{track}_all_dzdt.mat'
    da = pd.DataFrame({'x':[],
                       'pass_num':[],
                      })
    psys = []
    h = []
    UTCtime = []
    dh = []

    with h5py.File(path, 'r') as f:
        
        bins = list(f['dzdt'][track].keys())
        
        #make dataframe
        for b in tqdm(bins):
            
            x_arr = f['dzdt'][track][b]['x'][0]
            x_len = x_arr.shape[0]
            x_list = list(x_arr)
        
            grad = (list(f['dzdt'][track][b]['grad'][0])[0],list(f['dzdt'][track][b]['grad'][1])[0])
            df_temp =  pd.DataFrame({'x': x_list,
                                    'y': list(f['dzdt'][track][b]['y'][0]),
                                    'bin_number': [b]*x_len,
                                    'dzdt': list(f['dzdt'][track][b]['dzdt'][0])*x_len,
                                    'elev_wgs84_tpxtide': list(f['dzdt'][track][b]['elev_wgs84_tpxtide'][0]),
                                    'grad':[grad]*x_len,
                                    'ret': list(f['dzdt'][track][b]['ret'][0])*x_len,
                                    'retide': list(f['dzdt'][track][b]['retide'][0]),
                                    'sdatedays': list(f['dzdt'][track][b]['sdatedays'][0]),
                                    'z': list(f['dzdt'][track][b]['z'][0])
                                   })
 
            da_temp = df_temp.query("x > -382064.5722209641 & x < -374781.1654740692 & y > -734075.0820404041 & y < -722764.4514729496")
            
            da = da.append(da_temp, ignore_index=True )
        
            del df_temp, da_temp
            
    da.reset_index(drop=True,inplace=True)
    
    da.to_pickle(f'/Users/home/whitefar/DATA/REMOTE_SENSING/ICESAT1/dataframe_{track}_alldzdt.pkl')
    
# #all the paths of relevant .mat icesat1 files
# paths = glob.glob('/Users/home/horganhu/ICESAT_LINK/TAMATA_ICESAT/GLA12_633_DZDT/*_all_dzdt.mat')
# tracks = [path[-22:-13] for path in paths]
# for track in tqdm(tracks):
#     icesat1_alldzdt_todataframe(track)
    

# track = 'track0099'
# b = da.bin_number.unique()[0]

def icesat1_load_alldzdt_pickle(track):
    """
    It was too slow to do this at the same time as reading files
    """
        
    path = f'/Users/home/whitefar/DATA/REMOTE_SENSING/ICESAT1/dataframe_{track}_alldzdt.pkl'
    
    da = pd.read_pickle(path)
    
    #convert sdatedays UTCtime to a pd.Timestamp object. # see work_out_utc_time_to_a_timestamp.ipynb for more details
    da['UTCtime'] = da.sdatedays*86400 - (730486+0.5)*86400 
    da['timestamp'] = [pd.Timestamp.utcfromtimestamp(t) + 
                       (pd.Timestamp(2000,1,1,12) - pd.Timestamp(1970,1,1))
                       for t in da.UTCtime.to_list()
                       ]
    da['zp'] = np.nan
    da['residual'] = np.nan
    da['dz'] = np.nan
    da['x_bin'] = np.nan
    da['y_bin'] = np.nan
    
  
    for b in da.bin_number.unique():
        
        df_temp = da[da.bin_number==b].copy()
        
        df_len = df_temp.shape[0]
            
        df_temp['zp'] = ( np.matmul( np.vstack([df_temp.x.to_numpy()-df_temp.x.mean(),df_temp.y.to_numpy()-df_temp.y.mean()]).T,np.array(df_temp.grad.iloc[0]))    # rough size array([ 0.29046047,  0.12173804, -0.20730361, -0.30691026,  0.31946077,       -0.4548414 ,  0.237396  ])
                         + (( (df_temp.timestamp-df_temp.timestamp.iloc[0]) /  np.timedelta64(1, 'Y'))*df_temp.dzdt.iloc[0] # size (0   -0.0000001   -0.000569 2   -0.001384 3   -0.001777 4   -0.001999 5   -0.002551 6   -0.003171 )
                            + df_temp.z.mean() ).to_numpy() 
                         )
            
        df_temp['residual'] = df_temp.z - df_temp.zp;     # transient changes (footprint elevation changes corrected for
                                            # slope and secular change. Note here we use z here (not elev)

        df_temp['dz'] = df_temp.residual + ( (df_temp.timestamp-df_temp.timestamp.iloc[0]) /  np.timedelta64(1, 'Y'))*df_temp.dzdt.iloc[0];   # Elevation change corrected for each segment. This is
                # corrected for gradient. Note, this is per year.
        
        #add one coordinate per bin, and the rest nans
        df_temp['x_bin'] = [df_temp.x.mean()]+[np.nan]*(df_len-1) 
        df_temp['y_bin'] = [df_temp.y.mean()]+[np.nan]*(df_len-1)

        da[da.bin_number==b] = df_temp.copy()
    
    points = [Point(xy) for xy in zip(da.x,da.y)]
    gda = gpd.GeoDataFrame(da,geometry=points,crs=3031)
    
    # distance_cum = []
    # for pass_date in gda.timestamp.dt.date.unique():
    #     da_date = gda[gda.timestamp.dt.date==pass_date].copy()
    #     #make new column of along track distance
    #     #temp distance from previous
    #     tmp_dx = [Point.distance(da_date.geometry.iloc[i]) for i,Point in enumerate(da_date.geometry.iloc[1:])] #note the 1:, equivalent to i+1
    #     tmp_dx[:0] = [0]
    #     da_date['dx'] = tmp_dx
    #     distance_cum.extend( da_date.dx.cumsum().tolist() )
    # gda['distance_cum'] = distance_cum
            
    return gda


# =============================================================================
# 



def calculate_zp(da):
    """
    this is the same as above but for a dataframe.
    For use after smooth_grad
    """
 
    for b in da.bin_number.unique():
        
        df_temp = da[da.bin_number==b].copy()
    
        grads = np.array([[grad[0],grad[1]] for grad in df_temp.grad.tolist()]).T
        df_temp['zp'] = ( np.matmul( np.vstack([df_temp.x.to_numpy()-df_temp.x.mean(),
                                                             df_temp.y.to_numpy()-df_temp.y.mean()]).T,
                                                 grads)[:,0] #the first column of matrix product
                                       +
                                       (( (df_temp.timestamp-df_temp.timestamp.iloc[0]) /  np.timedelta64(1, 'Y'))*df_temp.dzdt.iloc[0] 
                                        + df_temp.z.mean() ).to_numpy() 
                                       )
            
        df_temp['residual'] = df_temp.z - df_temp.zp;     # transient changes (footprint elevation changes corrected for
                                            # slope and secular change. Note here we use z here (not elev)

        df_temp['dz'] = df_temp.residual + ( (df_temp.timestamp-df_temp.timestamp.iloc[0]) /  np.timedelta64(1, 'Y'))*df_temp.dzdt.iloc[0];   # Elevation change corrected for each segment. This is
                # corrected for gradient. Note, this is per year.
        
                
        
        da[da.bin_number==b] = df_temp.copy()
    
    points = [Point(xy) for xy in zip(da.x,da.y)]
    gda = gpd.GeoDataFrame(da,geometry=points,crs=3031)
                
    return gda
# =============================================================================



def smooth_grad(da,window_length  = 5,degree=2):
    """
    

    Parameters
    ----------
    da : pandas dataframe

    Returns
    -------
    pandas dataframe, where ['grad'] has been smoothed over each day

    """  
    
    for pass_date in da.timestamp.dt.date.unique():

        da_date = da[da.timestamp.dt.date==pass_date].copy()
        grads = np.array([ [grad[0],grad[1]] for grad in da_date.grad.tolist() ])
        if da_date.shape[0]<window_length:
            window_length = da_date.shape[0]
        smooth_grads = [(a,b) for a,b in zip(savgol_filter(grads[:,0],window_length,degree), savgol_filter(grads[:,1],window_length,degree)) ] 
        da_date['grad'] = smooth_grads
        da[da.timestamp.dt.date==pass_date] = da_date.copy()
    
    return da
        


def save_zps(track,x_or_y):
    """
    

    Parameters
    ----------
    track : track name e.g. 'track0211'
    x_or_y : whether to interpolate x or y coordinate (choose the one the line traverses most)

    Saves
    -------
    zps_{track}.pkl a pickle of the dataframe with zp and dzpdt_rate  for different passes

    """
    
    #smoothing
    window = 3
    degree = 1
    
    if x_or_y == 'x':
        x_or_y_other = 'y'
    elif x_or_y == 'y':
        x_or_y_other = 'x'
    
    da = icesat1_load_alldzdt_pickle(track)
    
    #GET COORDS
    #the actual coords of all zps run through the bin centres. We'll make a line through the bin centres then
    #use spacing from the first pass (zp_t0), smoothed
    zp_t0 = da[da.timestamp.dt.date==min(da.timestamp.dt.date.unique())].copy()
    coord_index = zp_t0.x_bin.dropna().index
    coord_line = LineString(list( zip(zp_t0.x_bin.loc[coord_index].tolist(),zp_t0.y_bin.loc[coord_index].tolist()) ) )
    coord_line = scale(coord_line,1.5,1.5)
    coord_bins = {'x': coord_line.coords.xy[0],'y':coord_line.coords.xy[1]} #dictionary so that x_or_y calls correct coords.
    ###
    zp_t0[x_or_y_other+'_smooth'] = savgol_filter(zp_t0[x_or_y_other],window,degree)
    zp_t0[x_or_y+'_smooth'] = savgol_filter(zp_t0[x_or_y],window,degree)
    
    f_bincoords = interpolate.interp1d(coord_bins['x'], coord_bins['y'],fill_value="extrapolate")
    zp_t0[x_or_y_other] = f_bincoords(zp_t0[x_or_y+'_smooth'])
    
    #Set zero line for differences
    z_0 = zp_t0.zp.to_numpy() 
    z_00 = z_0.copy() #00 for cumulatulative difference
    date_0 = zp_t0.timestamp.dt.date.iloc[0]
    date_00 = date_0
    
    for pass_date in da.sort_values(['timestamp'],axis=0).timestamp.dt.date.unique()[1:]:
        
        #get dataframe restricted to the area, restricted to a cetrain date
        da_date = da[da.timestamp.dt.date==pass_date].copy() 
        if da_date.shape[0] < len(zp_t0[x_or_y]) - 5:            #discard dates with too few points
            print(f'not enough data for {pass_date}, only {da_date.shape[0]} points')
            continue
        
        #interpolate zp so that we can colocate points and suptract
        da_date['x_smooth'] = savgol_filter(da_date.x,window,degree)
        da_date['y_smooth'] = savgol_filter(da_date.y,window,degree)
        da_date['zp_smooth'] = savgol_filter(da_date.zp,window,degree)
        
        f = interpolate.interp1d(da_date[x_or_y + '_smooth'], da_date.zp_smooth,fill_value="extrapolate")
        z_1 = f(zp_t0[x_or_y + '_smooth'])
  
        zp_t0[f"zp_{pass_date}"] = z_1
        
        date_1 = pass_date
    
        zp_t0[f"dzdt_{pass_date}"]  =( (z_1 - z_0 ) / (date_1 - date_0).days /365.2422 )
        
        zp_t0[f"dzdt_cum_{pass_date}"] = ( (z_1 - z_00 ) / (date_1 - date_00).days /365.2422 )
            
        z_0 = z_1.copy()
        date_0 = date_1
    
        del da_date
    
    zp_t0.to_pickle(f'/Users/home/whitefar/DATA/REMOTE_SENSING/ICESAT1/zps_{track}smoothpoints.pkl')
    print(f'Pickle of {track} saved')

    #Also save to shapefile for use with icesat1 over REMA 
    zp_t0.drop(['timestamp','grad'],axis=1,inplace=True)
    points = [Point(xy) for xy in zip(zp_t0.x_smooth,zp_t0.y_smooth)]
    gda = gpd.GeoDataFrame(zp_t0,geometry=points,crs="EPSG:3031")
    gda.to_file(f'/Users/home/whitefar/DATA/REMOTE_SENSING/ICESAT1/shapefiles_of_icesat1_over_channel/{track}smoothpoints.shp')

# =============================================================================
save_zps('track0099','y')
save_zps('track0211','x')
