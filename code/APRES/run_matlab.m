addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))




tt =['/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/',test]


fmcw_plot(t2,'maxrange',800)

t1 = '/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/2019-12-07_225600.dat'
t2 = '/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/2019-12-21_215758.dat'
depth = 712
global cfg
cfg = fmcw_process_config_vsr;

cfg.bedSearchRange = [depth-10 depth+10];
cfg.maxDepthConfig = depth-30;
cfg.maxRange = depth+30;


melt = fmcw_melt2(t1,t2)

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
