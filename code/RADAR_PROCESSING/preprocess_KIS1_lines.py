#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:36:30 2020

@author: whitefar
"""
# =============================================================================
#           

# THIS SCRIPT PROCESSES THE RADARLINES AT KIS1 and at the seismic lines

# #RADAR METADATA        
##         
#seiswp5_seiswp6 2019-12-18 18:34 19:41 14240 06352053256 seis12
#seiswp3_seiswp4 2019-12-18 15:15 17:36 20692 06352021458 seis34
#seiswp1_seiswp2 2019-12-18 11:28 14:14 22359 06351222831 seis56
#startn250m_endn250m 2019-12-17 15:34 16:38 13715 06351043316 linen250
#end500m_end250m_start250m_startn250m 2019-12-17 12:45 14:21 16825 06351014411 line250
#start500m_end500m 2019-12-17 10:52 12:09 14871 06350235133 line500
#start0m_kis1 2019-12-14 13:40 13:56 6568 06348013919 line0a
#startn500m_start0m 2019-12-14 13:31 13:38 4261 06348013011 linestart
#endn500m_startn500m 2019-12-14 12:12 13:13 13163 06348001055 linen500
#end0_endn500m 2019-12-14 12:00 12:07 4373 06347235903 lineend
#kis1_end0 2019-12-14 10:45 11:49 13480 06347224428 line0b
# =============================================================================


# ###############    2019-12-18     

#++++++++++++++++++++++++++++++++++++++++++++        
#seiswp5_seiswp6 2019-12-18 18:34 19:41 14240 06352053256 seis56


surveyseis56 = radarsurvey("06352053256")
surveyseis56.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyseis56.load_gnss_data()
surveyseis56.interpolate_gnss()

surveyseis56.refine_timesync('-182 seconds')
surveyseis56.split_lines_choose(moving_threshold=1,window = 3,threshold_type='velocity',plots = False)
# surveyseis56.split_lines_plot(["lineseis56"])
lineseis56dict = surveyseis56.split_lines_output()

lineseis56 = radarline(lineseis56dict,'lineseis56')
lineseis56.stack_spatially()
lineseis56.detrend_data()
# lineseis56.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
lineseis56.export_segy()

#seiswp3_seiswp4 2019-12-18 15:15 17:36 20692 06352021458 seis34

surveyseis34 = radarsurvey("06352021458")
surveyseis34.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyseis34.load_gnss_data()
surveyseis34.interpolate_gnss()

surveyseis34.refine_timesync('-179 seconds')
surveyseis34.split_lines_choose(moving_threshold=1,window = 3,plots = False)
surveyseis34.plots = Falseot(["lineseis34a","lineseis34b",'dud'])
lineseis34adict,lineseis34bdict,_ = surveyseis34.split_lines_output()

lineseis34dict = {'radata': pd.concat([lineseis34adict['radata'],lineseis34adict['radata']],0),
              'ch0': np.concatenate([lineseis34adict['ch0'],lineseis34bdict['ch0']],0),
              'ch1': np.concatenate([lineseis34adict['ch1'],lineseis34bdict['ch1']],0),
              'info': lineseis34adict['info']  }

lineseis34 = radarline(lineseis34dict,'lineseis34')
lineseis34.stack_spatially()
lineseis34.detrend_data()
# lineseis34.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
lineseis34.export_segy()

#seiswp1_seiswp2 2019-12-18 11:28 14:14 22359 06351222831 seis12


surveyseis12 = radarsurvey("06351222831")
surveyseis12.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyseis12.load_gnss_data()
surveyseis12.interpolate_gnss()

surveyseis12.refine_timesync('-180 seconds')
surveyseis12.split_lines_choose(moving_threshold=1,window = 3,plots = False)
# surveyseis12.plots = Falseot(["dud","lineseis12a","lineseis12b"])
_,lineseis12adict,lineseis12bdict = surveyseis12.split_lines_output()

lineseis12dict = {'radata': pd.concat([lineseis12adict['radata'],lineseis12bdict['radata']],0),
              'ch0': np.concatenate([lineseis12adict['ch0'],lineseis12bdict['ch0']],0),
              'ch1': np.concatenate([lineseis12adict['ch1'],lineseis12bdict['ch1']],0),
              'info': lineseis12adict['info']  }

lineseis12 = radarline(lineseis12dict,'lineseis12')
lineseis12.stack_spatially()
lineseis12.detrend_data()
# lineseis12.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

lineseis12.export_segy()


# ###############  
# ###############    2019-12-17    

#startn250m_endn250m 2019-12-17 15:34 16:38 13715 06351043316 linen250

surveyn250 = radarsurvey("06351043316")
surveyn250.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyn250.load_gnss_data()
surveyn250.interpolate_gnss()

surveyn250.refine_timesync('-12 seconds')
surveyn250.split_lines_choose(moving_threshold=1,window = 3,plots = False)
# surveyn250.plots = Falseot(["linen250"])
linen250dict = surveyn250.split_lines_output()

linen250 = radarline(linen250dict,'linen250')
linen250.stack_spatially()
linen250.detrend_data()
# linen250.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
linen250.export_segy()

#end500m_end250m_start250m_startn250m 2019-12-17 12:45 14:21 16825 06351014411 line250

survey250 = radarsurvey("06351014411")  #NOT SURE THIS IS THE RIGHT FILECODE
survey250.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
survey250.load_gnss_data()
survey250.interpolate_gnss()

survey250.refine_timesync('-12 seconds')
survey250.split_lines_choose(moving_threshold=1.8,window = 8,plots = False)
# survey250.plots = Falseot(["dud","dud","line250","dud","dud","dud","dud","dud","dud"])
_,_,line250dict,_,_,_,_,_,_ = survey250.split_lines_output()

line250 = radarline(line250dict,'line250')
line250.stack_spatially()
line250.detrend_data()
# line250.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
line250.export_segy()

#start500m_end500m 2019-12-17 10:52 12:09 14871 0635025133 line500


survey500 = radarsurvey("06350235133")
survey500.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
survey500.load_gnss_data()
survey500.interpolate_gnss()

survey500.refine_timesync('-12 seconds')
survey500.split_lines_choose(moving_threshold=1.9,window = 5,plots = False)
# survey500.plots = Falseot(["dud","line500","dud","dud","dud","dud","dud","dud"])
line500dict = survey500.split_lines_output()[1]

line500 = radarline(line500dict,'line500')
line500.stack_spatially()
line500.detrend_data()
# line500.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
line500.export_segy()



# ###############  
# ###############    2019-12-14   


#startn500m_start0m 2019-12-14 13:31 13:38 4261 06348013011 linestart
#done
surveystart = radarsurvey("06348013011")
surveystart.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveystart.load_gnss_data()
surveystart.interpolate_gnss()

surveystart.refine_timesync('4 seconds')
surveystart.split_lines_choose(moving_threshold=1,window = 3,plots = False)
# surveystart.plots = Falseot(["1"])
linestartdict = surveystart.split_lines_output()

linestart = radarline(linestartdict,'linestart')
linestart.clip_line_choose(clip_start_by=100,clip_end_by=65)
linestart.clip_line()

linestart.stack_spatially()
linestart.detrend_data()
# linestart.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
linestart.export_segy()


#endn500m_startn500m 2019-12-14 12:12 13:13 13163 06348001055 linen500
#done
surveyn500 = radarsurvey("06348001055")
surveyn500.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyn500.load_gnss_data()
surveyn500.interpolate_gnss()

surveyn500.refine_timesync('4 seconds')
surveyn500.split_lines_choose(moving_threshold=0.25,window = 5,plots = False)
# surveyn500.plots = Falseot(["1"])
linen500dict = surveyn500.split_lines_output()

linen500 = radarline(linen500dict,'linen500')
linen500.stack_spatially()
linen500.detrend_data()
# linen500.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
linen500.export_segy()

#end0_endn500m 2019-12-14 12:00 12:07 4373 06347235903 lineend
#done
surveyend = radarsurvey("06347235903")
surveyend.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyend.load_gnss_data()
surveyend.interpolate_gnss()

surveyend.refine_timesync('4 seconds')
surveyend.split_lines_choose(moving_threshold=1.61,window = 5,plots = False)
# surveyend.plots = Falseot(["dud","yep","dud"])
_,lineenddict,_ = surveyend.split_lines_output()

lineend = radarline(lineenddict,'lineend')
lineend.clip_line_choose(clip_start_by=10,clip_end_by=1)
lineend.clip_line()
lineend.stack_spatially()
lineend.detrend_data()
# lineend.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
lineend.export_segy()



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


#kis1_end0 2019-12-14 10:45 11:49 13480 06347224428 line0b
#done
survey0b = radarsurvey("06347224428")
survey0b.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
survey0b.load_gnss_data()
survey0b.interpolate_gnss()

survey0b.refine_timesync('4 seconds')
survey0b.split_lines_choose(moving_threshold=1,window = 3,plots = False)
# survey0b.plots = Falseot(["dud","0","1","2","3"])
_,dict0,dict1,dict2,dict3 = survey0b.split_lines_output()

line0bdict = {'radata': pd.concat([dict0['radata'],dict1['radata'],dict2['radata'],dict3['radata']],0),
              'ch0': np.concatenate([dict0['ch0'],dict1['ch0'],dict2['ch0'],dict3['ch0']],0),
              'ch1': np.concatenate([dict0['ch1'],dict1['ch1'],dict2['ch1'],dict3['ch1']],0),
              'info': dict0['info']  }

line0b = radarline(line0bdict,'line0b')
line0b.stack_spatially()
line0b.detrend_data()
# line0b.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

#line 0
#done
line0dict = {'radata': pd.concat([line0adict['radata'],line0bdict['radata']],0),
              'ch0': np.concatenate([line0adict['ch0'],line0bdict['ch0']],0),
              'ch1': np.concatenate([line0adict['ch1'],line0bdict['ch1']],0),
              'info': line0adict['info']  }



line0 = radarline(line0dict,'line0')
line0.stack_spatially()
line0.detrend_data()
# line0.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


line0.export_segy()