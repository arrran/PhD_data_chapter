#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 09:47:09 2020

@author: arran
"""
# https://github.com/emanuelhuber/RGPR/blob/12578623e9e70744b250abd3ac2a51f1f81880a4/tests/testthat/test_arg_checking.R
# test_that("dewow > OK", {
#   expect_true(class(dewow(x, type = "runmed", w = 20)) == "GPR")
#   expect_true(class(dewow(x, type = "runmean", w = 20)) == "GPR")
#   expect_true(class(dewow(x, type = "gaussian", w = 20)) == "GPR")
#   expect_true(class(dewow(x, type = "Gaussian", w = 20)) == "GPR")
#   expect_warning(dewow(x, type = "mad", w = 20))
#   expect_true(class(dewow(x)) == "GPR")
#   expect_warning(dewow(x, type = "MAD", w = 20))
#   expect_length(dewow(x)@proc, 1)
#   expect_length(dewow(x, track = FALSE)@proc, 0)
# })

# Dewow

# Remove the low-frequency components (the so-called “wow”) of the GPR record using:

#     a running median filter (type = "runnmed")
#     a running mean filter (type = "runmean")
#     a Gaussian filter (type = "Gaussian")

# For the two first cases, the argument w is the length of the filter in time units. For the Gaussian filter, w is the standard deviation.

# x3 <- dewow(x2, type = "runmed", w = 50)     # dewowing:
# plot(x3)                                     # plot the result

# https://github.com/NSGeophysics/GPRPy/blob/master/gprpy/gprpy.py

#     def dewow(self,window):
#         '''
#         Subtracts from each sample along each trace an 
#         along-time moving average.
#         Can be used as a low-cut filter.
#         INPUT:
#         window     length of moving average window 
#                    [in "number of samples"]
#         '''
#         # Store previous state for undo
#         self.storePrevious()
#         self.data = tools.dewow(self.data,window)
#         # Put in history
#         histstr = "mygpr.dewow(%d)" %(window)
#         self.history.append(histstr)

# def dewow(data,window):
#     '''
#     Subtracts from each sample along each trace an 
#     along-time moving average.
#     Can be used as a low-cut filter.
#     INPUT:
#     data       data matrix whose columns contain the traces 
#     window     length of moving average window 
#                [in "number of samples"]
#     OUTPUT:
#     newdata    data matrix after dewow
#     '''
#     totsamps = data.shape[0]
#     # If the window is larger or equal to the number of samples,
#     # then we can do a much faster dewow
#     if (window >= totsamps):
#         newdata = data-np.matrix.mean(data,0)            
#     else:
#         newdata = np.asmatrix(np.zeros(data.shape))
#         halfwid = int(np.ceil(window/2.0))
        
#         # For the first few samples, it will always be the same
#         avgsmp=np.matrix.mean(data[0:halfwid+1,:],0)
#         newdata[0:halfwid+1,:] = data[0:halfwid+1,:]-avgsmp

#         # for each sample in the middle
#         for smp in tqdm(range(halfwid,totsamps-halfwid+1)):
#             winstart = int(smp - halfwid)
#             winend = int(smp + halfwid)
#             avgsmp = np.matrix.mean(data[winstart:winend+1,:],0)
#             newdata[smp,:] = data[smp,:]-avgsmp

#         # For the last few samples, it will always be the same
#         avgsmp = np.matrix.mean(data[totsamps-halfwid:totsamps+1,:],0)
#         newdata[totsamps-halfwid:totsamps+1,:] = data[totsamps-halfwid:totsamps+1,:]-avgsmp
        
#     print('done with dewow')
#     return newdata

# Dewowing and standard bandpass filtering

# Many GPR data show a significantly very low 
# frequency component either due to inductivephenomena or possible instrumentation restrictions. 
# This low frequency range must be removedbefore applying any other digital filter algorithms. 
# There exist many different ways. A simple dewow filter acts within the time domain.
#  A running mean value is calculated for eachvalue of each trace. This running mean is subtracted 
#  from the central point. As filter parameter thetime range for the calculation of the running mean 
#  value must be entered which should be set toabout one or two principal periods. A possible static 
#  shift will also be removed using this filter. Alternatives to the dewow filter may be a high pass 
# bandpass filter working either within thefrequency or time domain or a simple subtract DC-shift filter 
# if only a constant value shall beremoved.



import numpy as np

def dewow(radata):
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
    data = np.asmatrix(radata.trData)
    totsamps = data.shape[0]
    window = radata.userArg
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
    return newdata










# data = np.loadtxt('/home/arran/PHD/DATA/RADAR/EXAMPLE_RADARLINE/radardata-line7-.csv').T

# window = 50

# from tqdm import tqdm

# for a in tqdm(range(10000)):
#     a**a
    
# import matplotlib.pyplot as plt

# newdata = dewow(data,20)

# bound=0.008
# fig, ax = plt.subplots(figsize=(12,12),dpi=180)
# ax.imshow(data[:1250,:],vmin=-bound, vmax=bound,aspect='auto'  )
# ax.xaxis.set_tick_params(rotation=90)
# ax.set_title('old')


# bound=0.008
# fig, ax = plt.subplots(figsize=(12,12),dpi=180)
# ax.imshow(newdata[:1250,:],vmin=-bound, vmax=bound,aspect='auto'  )
# ax.xaxis.set_tick_params(rotation=90)
