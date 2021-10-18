addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))

addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))

addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8622021/KIS Science/ApRES data files/KIS2'))

tt =['/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/',test]


fmcw_plot(t4,'maxrange',800)

t1 = '/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/2019-12-07_225600.dat'
t2 = '/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/2019-12-21_215758.dat'
t3 = '/Volumes/arc_04/FIELD_DATA/K8622021/KIS Science/ApRES data files/KIS2/KIS2 Old WO ApRES Card2/DIR2020-11-25-0213/DATA2020-12-24-2017.DAT'
t4 = 'Survey_2020-12-23_023455.dat'
t5 = '2019-12-08_022603.dat'

depth = 629;
global cfg
cfg.polyorder=1
cfg = fmcw_process_config_vsr;
cfg.bedSearchRange = [depth-0.1 depth+0.1];
cfg.maxDepthConfig = depth-30;
cfg.maxRange = depth+30;
cfg.doPlotMelt = 1
cfg.doPlotAll = 1
cfg.doSaveOutput = 0

fmcw_melt2(t5,t4);

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
