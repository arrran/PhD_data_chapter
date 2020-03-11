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


# ###############    2020-12-01      

#++++++++++++++++++++++++++++++++++++++++++++        

# #camp_L7p5_R7p5_R7p25_L7p25_L7p75_R7p75_camp 2020-01-01 10:07 11:39 15474 06001001502 surveycamp

surveycamp = radarsurvey("06001001502")
surveycamp.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveycamp.load_gnss_data()
surveycamp.interpolate_gnss()
surveycamp.refine_timesync('21 seconds')
surveycamp.split_lines_choose(moving_threshold=1.5,window = 3,threshold_type='velocity')
surveycamp.split_lines_plot(["dud","dud1","dud2","dud3","dud4","loop","line7p5","dud5","dudloop2","line7p25","loop2","dud6","line7p75","loop3","back2camp"])
_,_,_,_,_,_,line7p5dict,_,_,line7p25dict,_,_,line7p75dict,_,lineback2campdict = surveycamp.split_lines_output()
        
line7p5 = radarline(line7p5dict)
line7p5.detrend_data()
line7p5.density_profile()
line7p5.filter_data(High_Corner_Freq = 2.5e7)
line7p5.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line7p25 = radarline(line7p25dict)
line7p25.detrend_data()
line7p25.density_profile()
line7p25.filter_data(High_Corner_Freq = 2.5e7)
line7p25.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line7p75 = radarline(line7p75dict)
line7p75.detrend_data()
line7p75.density_profile()
line7p75.filter_data(High_Corner_Freq = 2.5e7)
line7p75.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

lineback2camp = radarline(lineback2campdict)
lineback2camp.detrend_data()
lineback2camp.filter_data(High_Corner_Freq = 2.5e7)
lineback2camp.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

# ###############    2019-12-31      

#++++++++++++++++++++++++++++++++++++++++++++        

# #camp_G0_G1_G2_G3 2019-12-31 20:53 22:28 15498 06001000235 surveyupchan

surveyupchan = radarsurvey("06001000235")
surveyupchan.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyupchan.load_gnss_data()
surveyupchan.interpolate_gnss()
surveyupchan.refine_timesync('-55 seconds')
surveyupchan.split_lines_choose(moving_threshold=0.5)
#surveyupchan.split_lines_plot(["lineupchan"])
lineupchan = radarline(surveyupchan.split_lines_output()[0],'lineupchan')

lineupchan.detrend_data()
lineupchan.stack_spatially()
lineupchan.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

lineupchan.export()

# #Cp01_Cp02_ddd_Cp11 2019-12-31 14:57 15:38 8374 06001000411 surveyAPREScross

#should start at utc 31 01:57 end at 31 02:35
        
surveyAPREScross = radarsurvey("06001000411")
surveyAPREScross.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
surveyAPREScross.load_gnss_data()
surveyAPREScross.interpolate_gnss()
surveyAPREScross.refine_timesync('-5 hours 59 minutes 17 seconds')
surveyAPREScross.split_lines_choose(moving_threshold=0.5)
surveyAPREScross.split_lines_plot(["lineAPREScross"])
lineAPREScross = radarline(surveyAPREScross.split_lines_output()[0])

lineAPREScross.detrend_data()
lineAPREScross.density_profile()
lineAPREScross.filter_data(High_Corner_Freq = 2.5e7)
lineAPREScross.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
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
survey5.refine_timesync('21 seconds')
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
survey3.refine_timesync('21 seconds')
survey3.split_lines_choose(moving_threshold=0.5)



survey3.split_lines_plot(["line3","loop","L35","loop"])
line3dict, _, left35dict, _ = survey3.split_lines_output()

line3 = radarline(line3dict)
line3.stack_spatially(stack_distance=4)
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
surveyAPRESdown.refine_timesync('21 seconds')
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
survey14.refine_timesync('21 seconds')
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

# =============================================================================
###############    2019-12-28   

#++++++++++++++++++++++++++++++++++++++++++++ 
# #L11_R11 2019-12-28 15:35 16:24 12052 06362023503 survey11

survey11 = radarsurvey("06362023503")
survey11.load_radar_data()
survey11.load_gnss_data()
survey11.interpolate_gnss()
survey11.refine_timesync('20 seconds')
survey11.split_lines_choose(moving_threshold=0.5)
survey11.split_lines_plot(["loop", "line11"])
 _, line11dict = survey11.split_lines_output()

line11 = radarline(line11dict)
line11.detrend_data()
line11.density_profile()
line11.filter_data(High_Corner_Freq = 2.5e7)
line11.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

#++++++++++++++++++++++++++++++++++++++++++++ 
# #R9_R10_L10_L11 2019-12-28 13:53 15:25 16662 06362005244 survey10


survey10 = radarsurvey("06362005244")
survey10.load_radar_data()
survey10.load_gnss_data()
survey10.interpolate_gnss()
survey10.refine_timesync('20 seconds')
survey10.split_lines_choose(moving_threshold=0.5)
survey10.split_lines_plot(["loop","dud", "right910","loop","line10",'loop',"left1011"])
_,_, right910dict,_,line10dict,_,left1011dict = survey10.split_lines_output()

right910 = radarline(right910dict)
right910.detrend_data()
right910.density_profile()
right910.filter_data(High_Corner_Freq = 2.5e7)
right910.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line10 = radarline(line10dict)
line10.detrend_data()
line10.density_profile()
line10.filter_data(High_Corner_Freq = 2.5e7)
line10.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

left1011 = radarline(left1011dict)
left1011.detrend_data()
left1011.density_profile()
left1011.filter_data(High_Corner_Freq = 2.5e7)
left1011.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

#++++++++++++++++++++++++++++++++++++++++++++ 
# #R7_L7_L9_R9 2019-12-28 10:49 12:36 17850 06361214828 survey79

survey79 = radarsurvey("06361214828")
survey79.load_radar_data()
survey79.load_gnss_data()
survey79.interpolate_gnss()
survey79.refine_timesync('20 seconds')
survey79.split_lines_choose(moving_threshold=0.5)
survey79.split_lines_plot(["dud","go", "line7","loop1","left79",'loop2',"line9"])
_,_, line7dict,_,left79dict,_,line9dict= survey79.split_lines_output()

line7 = radarline(line7dict)
line7.detrend_data()
line7.density_profile()
line7.filter_data(High_Corner_Freq = 2.5e7)
line7.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

left79 = radarline(left79dict)
left79.detrend_data()
left79.density_profile()
left79.filter_data(High_Corner_Freq = 2.5e7)
left79.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line9 = radarline(line9dict)
line9.detrend_data()
line9.density_profile()
line9.filter_data(High_Corner_Freq = 2.5e7)
line9.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

# =============================================================================
###############    2019-12-27

#++++++++++++++++++++++++++++++++++++++++++++ 
# #L6_R6_R8_L8_L10 2019-12-27 14:30 17:02 21166 06361013051 survey68


survey68 = radarsurvey("06361013051")
survey68.load_radar_data()
survey68.load_gnss_data()
survey68.interpolate_gnss()
survey68.refine_timesync('60 seconds')
survey68.split_lines_choose(moving_threshold=1.23)
survey68.split_lines_plot(["1","2", "line6","4","5","6","right68","8","line8","10","left810"])
_,_, line6dict,_,_,_,right68dict,_,line8dict,_,left810dict = survey68.split_lines_output()

line6 = radarline(line6dict)
line6.detrend_data()
line6.density_profile()
line6.filter_data(High_Corner_Freq = 2.5e7)
line6.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

right68 = radarline(right68dict)
right68.detrend_data()
right68.density_profile()
right68.filter_data(High_Corner_Freq = 2.5e7)
right68.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line8 = radarline(line8dict)
line8.detrend_data()
line8.density_profile()
line8.filter_data(High_Corner_Freq = 2.5e7)
line8.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

left810 = radarline(left810dict)
left810.detrend_data()
left810.density_profile()
left810.filter_data(High_Corner_Freq = 2.5e7)
left810.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


# =============================================================================
###############    2019-12-24

#++++++++++++++++++++++++++++++++++++++++++++ 
# #R4_L4_L6 2019-12-24 17:22 18:30 14205 06358042135 survey4

survey4 = radarsurvey("06358042135")
survey4.load_radar_data()
survey4.load_gnss_data()
survey4.interpolate_gnss()
survey4.refine_timesync('-10 seconds')
survey4.split_lines_choose(moving_threshold=0.5)
survey4.split_lines_plot(["line4","loop","left46"])
line4dict,_,left46dict= survey4.split_lines_output()

line4 = radarline(line4dict)
line4.detrend_data()
line4.density_profile()
line4.filter_data(High_Corner_Freq = 2.5e7)
line4.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

left46 = radarline(left46dict)
left46.detrend_data()
left46.density_profile()
left46.filter_data(High_Corner_Freq = 2.5e7)
left46.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')



#++++++++++++++++++++++++++++++++++++++++++++ 
# #L0_L2_R2_R4 2019-12-24 15:00 16:38 17215 06358015929 survey2


survey2 = radarsurvey("06358015929")
survey2.load_radar_data()
survey2.load_gnss_data()
survey2.interpolate_gnss()
survey2.refine_timesync('-10 seconds')
survey2.split_lines_choose(moving_threshold=0.5)
survey2.split_lines_plot(["0","left02","loop","line2","loop2","right24","loop3"])
_,left02dict,_,line2dict,_,right24dict,_ = survey2.split_lines_output()

left02 = radarline(left02dict)
left02.detrend_data()
left02.density_profile()
left02.filter_data(High_Corner_Freq = 2.5e7)
left02.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

line2 = radarline(line2dict)
line2.detrend_data()
line2.density_profile()
line2.filter_data(High_Corner_Freq = 2.5e7)
line2.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

right24 = radarline(right24dict)
right24.detrend_data()
right24.density_profile()
right24.filter_data(High_Corner_Freq = 2.5e7)
right24.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')



#++++++++++++++++++++++++++++++++++++++++++++ 
# #C0_R0_L0 2019-12-24 12:33 13:40 14067 06357233238 survey0

survey0 = radarsurvey("06357233238")
survey0.load_radar_data()
survey0.load_gnss_data()
survey0.interpolate_gnss()
survey0.refine_timesync('-12 seconds')
survey0.split_lines_choose(moving_threshold=0.5)
survey0.split_lines_plot(["dud","centretoright","loop","line0"])
_,_,_,line0dict = survey0.split_lines_output()

line0 = radarline(line0dict)
line0.detrend_data()
line0.density_profile()
line0.filter_data(High_Corner_Freq = 2.5e7)
line0.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


#++++++++++++++++++++++++++++++++++++++++++++ 
# #camp_C7_C6_ddd_C0 2019-12-24 10:52 12:29 16930 06357215137 surveydownchan   


surveydownchan = radarsurvey("06357215137")
surveydownchan.load_radar_data()
surveydownchan.load_gnss_data()
surveydownchan.interpolate_gnss()
surveydownchan.refine_timesync('-12 seconds')
surveydownchan.split_lines_choose(moving_threshold=0.5)
surveydownchan.split_lines_plot(["dud","linedownchan","loo"])
_,linedownchandict,_ = surveydownchan.split_lines_output()

linedownchan = radarline(linedownchandict)
linedownchan.detrend_data()
linedownchan.density_profile()
linedownchan.filter_data(High_Corner_Freq = 2.5e7)
linedownchan.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
