#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:36:30 2020

@author: whitefar
"""
# =============================================================================
#           
# #RADAR METADATA        
#         
# #camp_L7p5_R7p5_R7p25_L7p25_L7p75_R7p75_camp 2020-01-01 10:07 11:39 15474 06001001502 surveycamp
# #camp_G0_G1_G2_G3 2019-12-31 20:53 22:28 15498 06001000235 surveyupchan
# #Cp01_Cp02_ddd_Cp11 2019-12-31 14:57 15:38 8374 06001000411 surveyAPREScross
# #L5_R5 2019-12-30 16:51 17:36 11538 06364035101 survey5
# #R3_L3_L5 2019-12-30 15:05 16:29 15543 06364020457 survey3
# #Cp25_Cp24_ddd_Cp16_ddd_L1_R1_R3 2019-12-30 11:14 13:52 21704 06363221309 surveyAPRESdown
# #R14_L14_L15 2019-12-29 17:10 18:15 13752 06363041031 survey14
# #R11_R12_L12_L13_R13_R14 2019-12-29 13:49 16:13 20700 06363004826 survey1213
# #L11_R11 2019-12-28 15:35 16:24 12052 06362023503 survey11
# #R9_R10_L10_L11 2019-12-28 13:53 15:25 16662 06362005244 survey10
# #R7_L7_L9_R9 2019-12-28 10:49 12:36 17850 06361214828 survey79
# #L6_R6_R8_L8_L10 2019-12-27 14:30 17:02 21166 06361013051 survey68
# #R4_L4_L6 2019-12-24 17:22 18:30 14205 06358042135 survey4
# #L0_L2_R2_R4 2019-12-24 15:00 16:38 17215 06358015929 survey2
# #C0_R0_L0 2019-12-24 12:33 13:40 14067 06357233238 survey0
# #camp_C7_C6_ddd_C0 2019-12-24 10:52 12:29 16930 06357215137 surveydownchan        
#             
# =============================================================================

        


# ###############    2019-12-31      

#++++++++++++++++++++++++++++++++++++++++++++        

# #camp_G0_G1_G2_G3 2019-12-31 20:53 22:28 15498 06001000235 surveyupchan
surveyupchan = radarsurvey("06001000235")
surveyupchan.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyupchan.load_gnss_data()
surveyupchan.interpolate_gnss()
surveyupchan.refine_timesync('-21 seconds')
surveyupchan.split_lines_choose(moving_threshold=0.5)
surveyupchan.split_lines_plot(["line5"])
line5 = radarline(surveyupchan.split_lines_output()[0])

# #Cp01_Cp02_ddd_Cp11 2019-12-31 14:57 15:38 8374 06001000411 surveyAPREScross
        
survey5 = radarsurvey("06364035101")
survey5.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
survey5.load_gnss_data()
survey5.interpolate_gnss()
survey5.refine_timesync('-21 seconds')
survey5.split_lines_choose(moving_threshold=0.5)
survey5.split_lines_plot(["line5"])
line5 = radarline(survey5.split_lines_output()[0])
#

# =============================================================================
###############    2019-12-30      

#++++++++++++++++++++++++++++++++++++++++++++        

#    survey 5     one segment only
#        L5_R5 2019-12-30 16:51 17:36 11538 06364035101 survey5
        
survey5 = radarsurvey("06364035101")
survey5.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
survey5.load_gnss_data()
survey5.interpolate_gnss()
survey5.refine_timesync('-21 seconds')
survey5.split_lines_choose(moving_threshold=0.5)
survey5.split_lines_plot(["line5"])
line5 = radarline(survey5.split_lines_output()[0])
  
line5.detrend_data()
line5.density_profile()
line5.filter_data(High_Corner_Freq = 2.5e7)
line5.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

#++++++++++++++++++++++++++++++++++++++++++++ 
# #R3_L3_L5 2019-12-30 15:05 16:29 15543 06364020457 survey3

#survey 3

survey3 = radarsurvey("06364020457")
survey3.load_radar_data()
survey3.load_gnss_data()
survey3.interpolate_gnss()
survey3.split_lines_choose(moving_threshold=0.5)
survey3.refine_timesync('-21 seconds')
survey3.split_lines_choose(moving_threshold=0.5)



survey3.split_lines_plot(["line3","loop","L35","loop"])
line3dict, _, left35dict, _ = survey3.split_lines_output()

line3 = radarline(line3dict)
line3.detrend_data()
line3.density_profile()
line3.filter_data(High_Corner_Freq = 2.5e7)
line3.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


left35 = radarline(left35dict)
left35.detrend_data()
left35.density_profile()
left35.filter_data(High_Corner_Freq = 2.5e7)
left35.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


#++++++++++++++++++++++++++++++++++++++++++++ 
# #Cp25_Cp24_ddd_Cp16_ddd_L1_R1_R3 2019-12-30 11:14 13:52 21704 06363221309 surveyAPRESdown

#surveyAPRESdown

surveyAPRESdown = radarsurvey("06363221309")
surveyAPRESdown.load_radar_data()
surveyAPRESdown.load_gnss_data()
surveyAPRESdown.interpolate_gnss()
surveyAPRESdown.refine_timesync('-21 seconds')
surveyAPRESdown.split_lines_choose(moving_threshold=0.5)
surveyAPRESdown.split_lines_plot(['lineAPRESdown','dud1','shittyhalfline1','loop1','dud4','line1','dud2','loop3','right13','dud3'])

lineAPRESdowndict,_,_,_,_,line1dict,_,_,right13dict,_ = surveyAPRESdown.split_lines_output()

lineAPRESdown = radarline(lineAPRESdowndict)
lineAPRESdown.detrend_data()
lineAPRESdown.density_profile()
lineAPRESdown.filter_data(High_Corner_Freq = 2.5e7)
lineAPRESdown.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line1 = radarline(line1dict)
line1.detrend_data()
line1.density_profile()
line1.filter_data(High_Corner_Freq = 2.5e7)
line1.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

right13 = radarline(right13dict)
right13.detrend_data()
right13.density_profile()
right13.filter_data(High_Corner_Freq = 2.5e7)
right13.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


# =============================================================================
###############    2019-12-29     


#++++++++++++++++++++++++++++++++++++++++++++ 
# #R14_L14_L15 2019-12-29 17:10 18:15 13752 06363041031 survey14

survey14 = radarsurvey("06363041031")
survey14.load_radar_data()
survey14.load_gnss_data()
survey14.interpolate_gnss()
survey14.refine_timesync('-21 seconds')
survey14.split_lines_choose(moving_threshold=0.5)
survey14.split_lines_plot(['loop1','line14','loop2','left1415'])

_,line14dict,_,left1415dict = survey14.split_lines_output()

line14 = radarline(line14dict)
line14.detrend_data()
line14.density_profile()
line14.filter_data(High_Corner_Freq = 2.5e7)
line14.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

left1415 = radarline(left1415dict)
left1415.detrend_data()
left1415.density_profile()
left1415.filter_data(High_Corner_Freq = 2.5e7)
left1415.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')



#++++++++++++++++++++++++++++++++++++++++++++ 
# #R11_R12_L12_L13_R13_R14 2019-12-29 13:49 16:13 20700 06363004826 survey1213


survey1213 = radarsurvey("06363004826")
survey1213.load_radar_data()
survey1213.load_gnss_data()
survey1213.interpolate_gnss()
survey1213.refine_timesync('-20 seconds')
survey1213.split_lines_choose(moving_threshold=0.5)
survey1213.split_lines_plot(['dud','right1112','loop2','line12','loop3','left1213','loop4','dud','line13','loop','right1314'])

_,right1112dict,_,line12dict,_,left1213dict,_,_,line13dict,_,right1314dict = survey1213.split_lines_output()

right1112 = radarline(right1112dict)
right1112.detrend_data()
right1112.density_profile()
right1112.filter_data(High_Corner_Freq = 2.5e7)
right1112.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


line12 = radarline(line12dict)
line12.detrend_data()
line12.density_profile()
line12.filter_data(High_Corner_Freq = 2.5e7)
line12.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

left1213 = radarline(left1213dict)
left1213.detrend_data()
left1213.density_profile()
left1213.filter_data(High_Corner_Freq = 2.5e7)
left1213.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line13 = radarline(line13dict)
line13.detrend_data()
line13.density_profile()
line13.filter_data(High_Corner_Freq = 2.5e7)
line13.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

right1314 = radarline(right1314dict)
right1314.detrend_data()
right1314.density_profile()
right1314.filter_data(High_Corner_Freq = 2.5e7)
right1314.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


# #L11_R11 2019-12-28 15:35 16:24 12052 06362023503 survey11
# #R9_R10_L10_L11 2019-12-28 13:53 15:25 16662 06362005244 survey10
# #R7_L7_L9_R9 2019-12-28 10:49 12:36 17850 06361214828 survey79
# #L6_R6_R8_L8_L10 2019-12-27 14:30 17:02 21166 06361013051 survey68
# #R4_L4_L6 2019-12-24 17:22 18:30 14205 06358042135 survey4
# #L0_L2_R2_R4 2019-12-24 15:00 16:38 17215 06358015929 survey2
# #C0_R0_L0 2019-12-24 12:33 13:40 14067 06357233238 survey0
# #camp_C7_C6_ddd_C0 2019-12-24 10:52 12:29 16930 06357215137 surveydownchan        