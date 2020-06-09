#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:36:30 2020

@author: whitefar
"""
# =============================================================================

#start0m_kis1 2019-12-14 13:40 13:56 6568 06348013919 line0a
#done


survey0a = radarsurvey("06348013919")
survey0a.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
survey0a.load_gnss_data()
survey0a.interpolate_gnss()

survey0a.refine_timesync('4  seconds')
survey0a.split_lines_choose(moving_threshold=1.2,window = 5,plots = False)
# survey0a.plots = Falseot(["line0a","dud","dud"])
line0adict = survey0a.split_lines_output()[0]

line0a = radarline(line0adict,'line0a')
line0a.offset()
# line0a timestamp is monotonic

# line0a.stack_spatially()
# line0a.offset()

# line0a.export()

#kis1_end0 2019-12-14 10:45 11:49 13480 06347224428 line0b
#done
survey0b = radarsurvey("06347224428")
survey0b.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
survey0b.load_gnss_data()
survey0b.interpolate_gnss()

survey0b.refine_timesync('4 seconds')
survey0b.split_lines_choose(moving_threshold=1,window = 3,plots = False)
survey0b.split_lines_plot(["dud","0","1","2","3"])
_,dict0,dict1,dict2,dict3 = survey0b.split_lines_output()



line0bdict = {'radata': pd.concat([dict0['radata'],dict1['radata'],dict2['radata'],dict3['radata']],0),
              'ch0': np.concatenate([dict0['ch0'],dict1['ch0'],dict2['ch0'],dict3['ch0']],0),
              'ch1': np.concatenate([dict0['ch1'],dict1['ch1'],dict2['ch1'],dict3['ch1']],0),
              'info': dict0['info']  }



line0b = radarline(line0bdict,'line0bKIS1')


# =============================================================================
# 


from scipy.signal import savgol_filter
from shapely.affinity import translate


#for each point, find rolling "heading_angle"

gradients=( (line0b.radata.geometry.y.to_numpy()[1:] - line0b.radata.geometry.y.to_numpy()[:-1])
                             / (line0b.radata.geometry.x.to_numpy()[1:]- line0b.radata.geometry.x.to_numpy()[:-1]) ) 
gradients = np.hstack([gradients[0],gradients])

line0b.radata['raw_gradient'] = gradients
line0b.radata.raw_gradient.iloc[ 399] =   -0.58

offset_by = 27.185 #in metres

# mean_gradient = line0b.radata['raw_gradient'].rolling(window=window,center=True).mean().to_list()

# smoothed_gradient = [mean_gradient[8]]*int(window/2) + mean_gradient + [mean_gradient[-8]]*int(window/2)
window_gradient = 31
line0b.radata['smoothed_gradient'] = savgol_filter(line0b.radata['raw_gradient'],window_gradient , 3)
line0b.radata['theta'] = np.arctan(line0b.radata.smoothed_gradient.to_numpy())

offset_location = []
for i,row in line0b.radata.iterrows():
    offset_location.append(  translate(row.geometry ,
                                       xoff= offset_by*np.cos(row.theta) ,
                                       yoff= offset_by*np.sin(row.theta) ) )
line0b.radata['geometry'] = offset_location


#height
window_height=15
line0b.radata['height'] = savgol_filter(line0b.radata.height,window_height , 3)

# =============================================================================




# line0b.stack_spatially()
line0b.offset()
# line0b.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
# line0b.export()

#line 0
#done
line0dict = {'radata': pd.concat([line0adict['radata'],line0bdict['radata']],0),
              'ch0': np.concatenate([line0adict['ch0'],line0bdict['ch0']],0),
              'ch1': np.concatenate([line0adict['ch1'],line0bdict['ch1']],0),
              'info': line0adict['info']  }



line0KIS1 = radarline(line0dict,'line0KIS1_dc')

# line0.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line0KIS1.export()