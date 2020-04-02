#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 09:37:06 2020

@author: arran
"""

#this is just a draft script, fiddling with the process


# An option could be to load your SEGY file using obspy 
# (https://github.com/obspy/obspy/wiki)[https://github.com/obspy/obspy/wiki], 
# here is a description of how to do that 
# (https://docs.obspy.org/master/packages/obspy.io.segy.html)[https://docs.obspy.org/master/packages/obspy.io.segy.html], 
# and then export them as a depth - time binary file .dt1, you can then hand-write a header
#  in the style of the example data .h. This should allow you to load the data.







runfile('/home/arran/PHD/DATA/code/process_radar.py', wdir='/home/arran/PHD/DATA/code')

# start_date_utm = survey79.radata.datetime.iloc[0].strftime('%Y-%m-%d')
# end_date_utm = survey79.radata.datetime.iloc[-1].strftime('%Y-%m-%d')

# dates = ['2019-12-27']
# files_paths = []
    
# for date in dates:
#     input_files = f"/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/PPP/**/"+date+".pos"
#     files_paths = files_paths + glob.glob(input_files)

# #files_paths = '/Volumes/arc_04/FIELD_DATA/K8621920/GNSS/PROCESSED/PPP/full_output_2019-12-27/2019-12-27.pos'

# survey79.track_points = load_ppp_date([start_date_utm])




# =============================================================================
# ON CO522PC01
# #R7_L7_L9_R9 2019-12-28 10:49 12:36 17850 06361214828 survey79

survey79 = radarsurvey("06361214828")
survey79.load_radar_data()
survey79.load_gnss_data()
survey79.interpolate_gnss()
survey79.refine_timesync('20 seconds')
survey79.split_lines_choose(moving_threshold=0.5)
#survey79.split_lines_plot(["dud","go", "line7","loop1","left79",'loop2',"line9"])
_,_, line7dict,_,left79dict,_,line9dict= survey79.split_lines_output()

line7 = radarline(line7dict,'line7')
line7.stack_spatially()
line7.detrend_data()
line7.export_segy()

line9 = radarline(line9dict,'line9')
line9.stack_spatially()
line9.detrend_data()
line9.export_segy()


#https://docs.obspy.org/tutorial/code_snippets/anything_to_miniseed.html


import obspy as ob


# =============================================================================
# a single trace
# data = line7.ch0[0,:]
# 
# # Fill header attributes
# stats = {'station': 'PX2', 
#          'location': (str(line7.radata.geometry.x[0])+', '+str(line7.radata.geometry.y[0])+', '+str(line7.radata.height[0])),
#          'starttime': line7.radata.datetime.iloc[0].strftime('%Y-%m-%dT%H:%M:%SZ'),
#          'channel': 'ch0',
#          'npts': len(data),}
# 
# 
# test = Stream([Trace(data=data, header=stats)])
# 
# =============================================================================

traces = []

for i,data in enumerate(line7.ch0):

    # Fill header attributes
    stats = {'station': 'PX2', 
             'location': (str(line7.radata.geometry.x.iloc[i])+', '+str(line7.radata.geometry.y.iloc[i])+', '+str(line7.radata.height.iloc[i])),
             'starttime': line7.radata.datetime.iloc[i].strftime('%Y-%m-%dT%H:%M:%SZ'),
             'channel': 'ch0',
             'npts': len(data),}
    
    traces.append( Trace(data=data, header=stats) )

line7 = Stream(traces)
    


# =============================================================================









#TO LOAD LOCALLY
# #R7_L7_L9_R9 2019-12-28 10:49 12:36 17850 06361214828 survey79

survey79 = radarsurvey("06361214828",metadata_path = "/home/arran/PHD/DATA/RADAR/EXAMPLE_RADARLINE/radar_metadata",timesync_path = "/home/arran/PHD/DATA/RADAR/EXAMPLE_RADARLINE/time_sync")
survey79.load_radar_data(path = "/home/arran/PHD/DATA/RADAR/EXAMPLE_RADARLINE")
survey79.load_gnss_data(ppp_path = "/home/arran/PHD/DATA/RADAR")
survey79.interpolate_gnss()
survey79.refine_timesync('20 seconds',timesync_path = "/home/arran/PHD/DATA/RADAR/EXAMPLE_RADARLINE/time_sync")
survey79.split_lines_choose(moving_threshold=0.5)
#survey79.split_lines_plot(["dud","go", "line7","loop1","left79",'loop2',"line9"])
_,_, line7dict,_,left79dict,_,line9dict= survey79.split_lines_output()

line7 = radarline(line7dict,'line7')
line7.stack_spatially()
line7.detrend_data()
















#line7.density_profile()
#line7.filter_data(High_Corner_Freq = 2.5e7)
#line7.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
line7.export()

left79 = radarline(left79dict,'left79')
left79.stack_spatially()
left79.detrend_data()
#left79.density_profile()
#left79.filter_data(High_Corner_Freq = 2.5e7)
#left79.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
left79.export()

line9 = radarline(line9dict,'line9')
line9.stack_spatially()
line9.detrend_data()
#line9.density_profile()
#line9.filter_data(High_Corner_Freq = 2.5e7)
#line9.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
line9.export()