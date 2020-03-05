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
########  # 24-12-2019
        
#camp_C7_C6_ddd_C0 2019-12-24 10:52 12:29 16930 06357215137 surveydownchan 

surveydownchan = radarsurvey("06357215137")
surveydownchan.load_radar_data()
surveydownchan.load_gnss_data()
surveydownchan.interpolate_gnss()
surveydownchan.split_lines_choose(moving_threshold=1,window=100)
surveydownchan.split_lines_plot(["linedownchan","?","???"])
line3dict, _, L35dict, _ = surveydownchan.split_lines_output()




        
        
# 1st      
# camp_L7p5_R7p5_R7p25_L7p25_L7p75_R7p75_camp 2020-01-01 10:07 11:39 15474 06001001502
surveycamp = radarsurvey("06001001502") 
surveycamp.load_radar_data()
#print(f"offset varies by {surveycamp.time_offset_variation}")
#
## offset varies by -1 days +23:59:51.024000
## ie less than a minute
# 
#31st        
##camp_G0_G1_G2_G3 2019-12-31 00:00 22:28 15498 06001000235
#surveyAPRESdownchan= radarsurvey("06001000235")
#surveyAPRESdownchan.load_radar_data()
#print(f"offset varies by {surveyAPRESdownchan.time_offset_variation}")
#
## offset varies by -1 days +03:06:26.976000
## less than a minute
#
##Cp01_Cp02_ddd_Cp11 2019-12-31 14:57 15:38 8374 06001000411
#surveyAPREScrosschan= radarsurvey("06001000411")
#surveyAPREScrosschan.load_radar_data()
#print(f"offset varies by {surveyAPREScrosschan.time_offset_variation}")
#
## offset varies by 0 days 00:00:24.950400
##less than a minute
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
survey5.split_lines_choose()
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
survey3.split_lines_choose(moving_threshold=3)
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








#COMPARE INTERPOLATED POSITION TO POSITION
survey3 = radarsurvey("06364020457")
survey3.load_radar_data()
survey3.load_gnss_data()
survey3.extra_position_info()


plt.plot(survey3.track_points.datetime,survey3.track_points.velocity)






survey3.metadata

survey3.radata.datetime.iloc[0] 
survey3.radata.datetime.iloc[-1]

#the radar timestamp is recording 02:04:11 till 03:27:48
#the metadata says line from 02:05:00 till 03:29:00, so pretty accurate

#plot the line period on the raw GNSS data
plt.plot(np.hstack([-379400*np.ones(1470),survey3.track_points.geometry.x.to_numpy()[7745:12300]]))  #- index 7750 till 12295
#plot the line period on intepd radar timestamp
#plt.figure()
plt.plot(survey3.radata.geometry.x.to_numpy()[3635:15320:2])


#is the RADAR recording more and more peeps with time?

data = signal.detrend(survey3.ch0, axis=1, type='constant', bp=0)
fig, ax = plt.subplots(figsize=(12,12),dpi=180)
ax.imshow(data[:,:1250].T,vmin=-0.008, vmax=0.008,aspect='auto'  )
ax.xaxis.set_tick_params(rotation=90)

plt.plot(survey3.track_points.velocity)

#First period in raw radargram is from index 3598 to 11760
#second large period is from index 12670 to 15221
len1_rad = abs(3598 - 11760)
len2_rad = abs(12670 - 15221)

len1_rad/len2_rad
#First period in raw gnss is from index 7749 to 10370
#second large period is from index 10810 to 12289
len1_gnss = abs(7749 - 10370)
len2_gnss = abs(10810 - 12289)
len1_gnss/len2_gnss

#++++++++++++++++++++++++++++++++++++++++++++ 
# #Cp25_Cp24_ddd_Cp16_ddd_L1_R1_R3 2019-12-30 11:14 13:52 21704 06363221309 surveyAPRESdown

#surveyAPRESdown

surveyAPRESdown = radarsurvey("06363221309")
surveyAPRESdown.load_radar_data()
surveyAPRESdown.load_gnss_data()
surveyAPRESdown.interpolate_gnss()
surveyAPRESdown.split_lines_choose(moving_threshold=2,threshold_type = 'acc')
surveyAPRESdown.split_lines_plot(list(range(0,4))+['lineAPRESdown']+list(range(5,22)))
#here i have to use different thresholds to separate different sections so just doing APRESdown first
lineAPRESdowndict = surveyAPRESdown.split_lines_output()[0]

lineAPRESdown = radarline(lineAPRESdowndict)
lineAPRESdown.detrend_data()
lineAPRESdown.density_profile()
lineAPRESdown.filter_data(High_Corner_Freq = 2.5e7)
lineAPRESdown.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')

surveyAPRESdown.split_lines_choose(moving_threshold=2)
surveyAPRESdown.split_lines_plot(['dud1','dud2','line1','loop2','right13','loop3'])
dud1,dud2,line1dict,loop2,right13dict,loop3 = surveyAPRESdown.split_lines_output()


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

# #R14_L14_L15 2019-12-29 17:10 18:15 13752 06363041031 survey14
# #R11_R12_L12_L13_R13_R14 2019-12-29 13:49 16:13 20700 06363004826 survey1213



# =============================================================================

#surveyAPRESdown = radarsurvey("06363221309")
#surveyAPRESdown.load_radar_data()



     
# LINE 14         R14_L14_L15 
#line has 5 segments of moving - blip,turn,line14,turn,left14
#
#        
#survey14 = radarsurvey("06363041031")
#survey14.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#survey14.load_gps_data()
#survey14.interpolate_gps()
#survey14.split_lines_choose()
#survey14.split_lines_plot(names = ['blip','turn','line14','turn','left14'])
#blip,turn,line14dict,turn,left14dict = survey14.split_lines_output()
#
#line14 = radarline(line14dict)
#line14.detrend_data()
#line14.density_profile()
#line14.detrend_data()
#line14.filter_data(High_Corner_Freq = 2.5e7)
#line14.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
#
#left14 = radarline(left14dict)
#left14.detrend_data()
#left14.density_profile()
#left14.detrend_data()
#left14.filter_data(High_Corner_Freq = 2.5e7)
#left14.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')


# =============================================================================
#Cp25_Cp24_ddd_Cp16_ddd_L1_R1_R3  
#line has 6 segments of moving - blip, downapres,halfline1&turn,line1,turn,right13
#
#surveydownapres = radarsurvey("06363221309")
#surveydownapres.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#surveydownapres.load_gps_data()
#surveydownapres.interpolate_gps()
#surveydownapres.split_lines_choose()
#surveydownapres.split_lines_plot(names = ['blip','downapres','halfline1&turn','line1','turn','right3'])
#blip, downapresdict, halfline1turn,line1dict,turn,right3dict = surveydownapres.split_lines_output()
#
#line1= radarline(line1dict)
#line1.detrend_data()
#line1.density_profile()
#line1.detrend_data()
#line1.filter_data(High_Corner_Freq = 2.5e7)
#line1.radargram(channel=0,bound=0.008,title='line1 filtered to 2.5e7 Hz',x_axis='space')
#
#downapres= radarline(downapresdict)
#downapres.detrend_data()
#downapres.density_profile()
#downapres.detrend_data()
#downapres.filter_data(High_Corner_Freq = 2.5e7)
#downapres.radargram(channel=0,bound=0.008,title='downapres filtered to 2.5e7 Hz',x_axis='space')
#
#right3= radarline(right13dict)
#right3.detrend_data()
#right3.density_profile()
#right3.detrend_data()
#right3.filter_data(High_Corner_Freq = 2.5e7)
#right3.radargram(channel=0,bound=0.008,title='right3 filtered to 2.5e7 Hz',x_axis='space')


# =============================================================================

# =============================================================================

# =============================================================================
# R7_L7_L9_R9 2019-12-28 10:49 12:36 17850 06361214828
# line has 7 segments of moving, where acc < 2.8 - blip, blip, line7, loop,left79,loop,line9
#
        
#survey79 = radarsurvey("06361214828")
#survey79.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#survey79.load_gps_data()
#survey79.interpolate_gps()
#survey79.radata.plot()
#survey79.split_lines_choose(moving_threshold=2.8)
#survey79.detrend_data()
#survey79.radargram()
#survey79.split_lines_plot(names = ['blip', 'blip', 'line7', 'loop','left79','loop','line9'])
#blip, blip, line7dict, loop,left79dict,loop,line9dict  = survey79.split_lines_output()
#
#line7 = radarline(line7dict)
#line7.detrend_data()
#line7.density_profile()
#line7.detrend_data()
#line7.filter_data(High_Corner_Freq = 2.5e7)
#line7.radargram(channel=0,bound=0.008,title='line7 filtered to 2.5e7 Hz',x_axis='space')
#
#left79= radarline(left79dict)
#left79.detrend_data()
#left79.density_profile()
#left79.detrend_data()
#left79.filter_data(High_Corner_Freq = 2.5e7)
#left79.radargram(channel=0,bound=0.008,title='left79 filtered to 2.5e7 Hz',x_axis='space')
#
#line9= radarline(line9dict)
#line9.detrend_data()
#line9.density_profile()
#line9.filter_data(High_Corner_Freq = 2.5e7)
#line9.radargram(channel=0,bound=0.008,title='line9 filtered to 2.5e7 Hz',x_axis='space')
        
# =============================================================================
# C0_R0_L0 2019-12-24 12:33 13:40 14067 06357233238
#       line has 9 segments of moving, where acc < 1 - b1,b2,b3,b4,b5,b6,b7,line0,loop
#survey0 = radarsurvey("06357233238")
#survey0.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#survey0.load_gps_data()    
#survey0.interpolate_gps()
#survey0.radata.plot()
#survey0.detrend_data()
#survey0.radargram(channel=0,bound=0.008) 
#survey0.split_lines_choose(moving_threshold=1,window =55)
#survey0.split_lines_plot(names = ['b1','b2','b3','line0','loop'])
#b1,b2,b3,line0dict,loop  = survey0.split_lines_output()
##
#line0= radarline(line0dict)
#line0.detrend_data()
#line0.density_profile()
#line0.detrend_data()
#line0.filter_data(High_Corner_Freq = 2.5e7)
#line0.radargram(channel=0,bound=0.008,title='line0 filtered to 2.5e7 Hz',x_axis='space')

# =============================================================================
# 
##L0_L2_R2_R4 2019-12-24 15:00 16:38 17215 06358015929
#        
#
#survey2 = radarsurvey("06358015929")
#survey2.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#survey2.load_gps_data()       
#survey2.detrend_data()
#survey2.radargram(channel=0,bound=0.008) 
#survey2.interpolate_gps()
#survey2.split_lines_choose(moving_threshold=1)
#survey2.split_lines_plot(names = ['b1','b2','b3','b4','b5','b6','b7','line0','loop'])     
# =============================================================================


###
#plt.figure()
#plt.plot(survey0.track_points.time,survey0.track_points.geometry.y,'x')
#plt.xticks(rotation=90)
#plt.grid()      

#
#
#
#survey2 = radarsurvey("06358015929")
#survey2.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#survey2.load_gps_data()  
#survey2.interpolate_gps()
#survey2.radata.plot()
#survey2.split_lines_choose(moving_threshold=1)


#
#
#pd.Timestamp('28-12-2019T10:49:00') - pd.Timestamp('27-12-2019T21:47:35') 
#
#
#
#pd.Timestamp('30-12-2019T11:14:00') - pd.Timestamp('30-12-2019T02:04:12') 
















        
        
#        First need to sort timesync for the 24th
# C0_R0_L0 2019-12-24 12:33 13:40 14067 06357233238
# get the point where stopped at L_0 - 1577140507.7
#        
#
#survey0 = radarsurvey("06357233238")
#survey0.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#survey0.load_gps_data()       
#survey0.detrend_data()
#survey0.radargram(channel=0,bound=0.008)  
#plt.axvline(1577137937.7, color='k', linestyle='solid')
#plt.axvline(1577140507.7, color='k', linestyle='solid')
# 
#survey0.extra_gps()             
##
#plt.plot(survey0.track_points.datetime,survey0.track_points.geometry.y,'x')
#plt.xticks(rotation=90)
#plt.grid()         
#
#
#timesync_aa = pd.Timestamp('24-12-2019T00:10:20') - pd.Timestamp('23-12-2019T22:35:07.7')     
#timesync_ab = pd.Timestamp('24-12-2019T02:57:02.5') - pd.Timestamp('24-12-2019T01:21:26.1')
#
##L0_L2_R2_R4 2019-12-24 15:00 16:38 17215 06358015929
## sort timesync
#        
#survey2 = radarsurvey("06358015929")
#survey2.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#survey2.load_gps_data()       
#survey2.detrend_data()
#survey2.radargram(channel=0,bound=0.008)  
##plt.axvline(1577148753.0, color='k', linestyle='solid')
# plt.axvline(1577148754.45, color='k', linestyle='solid')
#survey2.extra_gps()             
##
#plt.plot(survey2.track_points.datetime,survey2.track_points.geometry.y,'x')
#plt.xticks(rotation=90)
#plt.grid()      
##        
#
#
#timesync_b = pd.Timestamp('24-12-2019T02:13:44') - pd.Timestamp('24-12-2019T00:52:34.45')      
# 1577150486.1
#        
##camp_C7_C6_ddd_C0 2019-12-24 10:52 12:29 16930 06357215137
#
#        
#surveydownchan = radarsurvey("06357215137")
#surveydownchan.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#surveydownchan.load_gps_data()       
#surveydownchan.detrend_data()
#surveydownchan.radargram(channel=0,bound=0.008)   
#plt.axvline(1577132409.5, color='k', linestyle='solid')
#plt.axvline(1577136124.7, color='k', linestyle='solid')
#
#surveydownchan.extra_gps()             
#
#plt.plot(surveydownchan.track_points.datetime,surveydownchan.track_points.geometry.x,'x')
#plt.xticks(rotation=90)
#plt.grid()      
#        
#timesync_c_start = pd.Timestamp('23-12-2019T21:59:54') - pd.Timestamp('23-12-2019T20:20:09.50')
##timesync_a_stop =  pd.Timestamp( '24-12-2019T00:10:20') - pd.Timestamp('23-12-2019T21:22:04.7') 
# 
#       
##R4_L4_L6 2019-12-24 17:22 18:30 14205 06358042135    
#
#surveydownchan = radarsurvey("06358042135")
#surveydownchan.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#surveydownchan.load_gps_data()       
#surveydownchan.detrend_data()
#surveydownchan.radargram(channel=0,bound=0.008)   
#plt.axvline(1577132409.5, color='k', linestyle='solid')
#
#timesync_d = pd.Timestamp('24-12-2019T05:09:08.5') - pd.Timestamp('24-12-2019T03:24:53.8')       
# =============================================================================

        
        
# =============================================================================
#Cp01_Cp02_ddd_Cp11 2019-12-31 14:57 15:38 8374 06001000411


#surveycrossapres = radarsurvey("06001000411")
#surveycrossapres.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#surveycrossapres.load_gps_data()
#surveycrossapres.detrend_data()
#
#surveycrossapres.extra_gps()
##surveycrossapres.interpolate_gps()
##surveycrossapres.split_lines_choose(moving_threshold=3)
##surveycrossapres.split_lines_plot(names = ['
#
#surveycrossapres.radargram(channel=0,bound=0.008)
#
#plt.plot(surveycrossapres.track_points.datetime,surveycrossapres.track_points.velocity)
#plt.xticks(rotation=90)
#plt.grid()







































#line5.radata.keys()
#
#line5.detrend_data()
#line5.filter_data(High_Corner_Freq = 2.5e7)
#line5.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz',x_axis='space')
##

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
#line14 = radarsurvey("06363041031")
#line14.load_radar_data("/Volumes/arc_04/FIELD_DATA/K8621920/RES/")
#line14.detrend_data()
#line14.density_profile()
#line14.filter_data(High_Corner_Freq = 2.5e7)
#line14.radargram(channel=0,bound=0.008,title='filtered to 2.5e7 Hz')
#
##28-12-2019
#line11 = radarsurvey()
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