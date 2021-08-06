addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))

t1 = '2019-12-07_225600.dat'
t2 = '2019-12-07_231044.dat'
t3 = '2019-12-07_232530.dat'
t4 = '2019-12-07_233943.dat'
t5 = '2019-12-07_235231.dat'
t6 = '2019-12-08_000409.dat'
t7 = '2019-12-08_001936.dat'
t8 = '2019-12-08_003236.dat'
t9 = '2019-12-08_013943.dat'
t10 = '2019-12-08_015742.dat'
t11 = '2019-12-08_021102.dat'
t12 ='2019-12-08_022603.dat'
t13 ='2019-12-08_024034.dat'
t14 ='2019-12-08_025312.dat'
t15 ='2019-12-08_030609.dat'
g1 = "2019-12-21_215758.dat"
g2 = "2019-12-21_221136.dat"
g3 = "2019-12-21_222317.dat"
g4 = "2019-12-21_223725.dat"
g5 = "2019-12-21_224933.dat"
g6 =  "2019-12-21_230256.dat" %craig did this one
g7 = "2019-12-21_232339.dat"
g8 =  "2019-12-22_015010.dat"
g9 = "2019-12-22_021325.dat"
g10 = "2019-12-22_022327.dat"
g11 = "2019-12-22_023324.dat"
g12 = "2019-12-22_020319.dat"
g13 = "2019-12-21_235003.dat"
g14 = "2019-12-21_231333.dat"
g15 = "2019-12-21_233833.dat"


f = strcat("/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/",t6)
g =  strcat("/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/",g6)
% 
% fmcw_plot(tt1,'maxrange',800)
% fmcw_plot('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/2019-12-21_215758.dat','maxrange',800)


depth = 454
global cfg
cfg.polyorder=1
cfg = fmcw_process_config_vsr;
cfg.bedSearchRange = [depth-0.1 depth+0.1];
cfg.maxDepthConfig = depth-30;
cfg.maxRange = depth+30;
cfg.doPlotMelt = 1
cfg.doPlotAll = 1
cfg.doSaveOutput = 0
ya = fmcw_melt2(f,g)

% 
% 
% fmcw_plot(filename1)
% 
% fmcw_plot({'Survey_2019-12-26_042742.dat'},'plotop','a','burstlist','all','chirplist','all','maxrange',3500)
% 
% vdat = fmcw_load(filename1)
% 
% filename1 = '/Volumes/arc_04/FIELD_DATA/K8621920/APRES/Survey_2019-12-08_001936.dat'
% filename2 = '/Volumes/arc_04/FIELD_DATA/K8621920/APRES/Survey_2019-12-21_232339.dat'
% fmcw_melt(filename1,filename2)
