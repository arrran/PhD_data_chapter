#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:36:30 2020

@author: whitefar
"""
# =============================================================================
#           

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

line7.radargram()
plt.show()

line7.dewow(window=1000)
line7.radargram()
plt.show()

line7.reset()
line7.dewow(window=100)
line7.radargram()
plt.show()

line7.reset()
line7.dewow(window=10)
line7.radargram()
plt.show()

line7.reset()
line7.dewow(window=1)
line7.radargram()
plt.show()


 def dewow(line7,window=100):
        '''
        Subtracts from each sample along each trace an 
        along-time moving average.
        Can be used as a low-cut filter.
        INPUT:
        data       data matrix whose columns contain the traces 
        window     length of moving average window 
                   [in "number of samples"]
        OUTPUT:
        newdata    data matrix after dewow
        '''
        data = np.asmatrix(line7.ch0.T)
        totsamps = data.shape[0]
        
        # If the window is larger or equal to the number of samples,
        # then we can do a much faster dewow
        if (window >= totsamps):
            newdata = data-np.matrix.mean(data,0)            
        else:
            newdata = np.asmatrix(np.zeros(data.shape))
            halfwid = int(np.ceil(window/2.0))
            
            # For the first few samples, it will always be the same
            avgsmp=np.matrix.mean(data[0:halfwid+1,:],0)
            newdata[0:halfwid+1,:] = data[0:halfwid+1,:]-avgsmp
    
            # for each sample in the middle
            for smp in range(halfwid,totsamps-halfwid+1):
                winstart = int(smp - halfwid)
                winend = int(smp + halfwid)
                avgsmp = np.matrix.mean(data[winstart:winend+1,:],0)
                newdata[smp,:] = data[smp,:]-avgsmp
    
            # For the last few samples, it will always be the same
            avgsmp = np.matrix.mean(data[totsamps-halfwid:totsamps+1,:],0)
            newdata[totsamps-halfwid:totsamps+1,:] = data[totsamps-halfwid:totsamps+1,:]-avgsmp
            
        print('done with dewow')
        line7.ch0 = newdata.T

line7.dewow()
line7.radargram









# # line7.detrend_data()
# #line7.density_profile()
# #line7.filter_data(High_Corner_Freq = 2.5e7)
# #line7.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
# line7.export_gis()

# left79 = radarline(left79dict,'left79')
# left79.stack_spatially()
# # left79.detrend_data()
# #left79.density_profile()
# #left79.filter_data(High_Corner_Freq = 2.5e7)
# #left79.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
# left79.export_gis()

# line9 = radarline(line9dict,'line9')
# line9.stack_spatially()
# # line9.detrend_data()
# #line9.density_profile()
# #line9.filter_data(High_Corner_Freq = 2.5e7)
# #line9.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
# line9.export_gis()
