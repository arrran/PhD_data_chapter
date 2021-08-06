addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))
addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8622021/KIS Science/ApRES data files/KIS2'))
meta = readtable('/Users/home/whitefar/DATA/ApRES/kis2_apres_metadata_by_Site_ID.csv');
data = meta;


%data.output(1) = ya;

global cfg
cfg.polyorder=1;
cfg = fmcw_process_config_vsr;

cfg.doPlotMelt = 0;
cfg.doPlotAll = 0;
cfg.doSaveOutput = 0;

for i = 1:height(meta)
    display(i)
    if string(meta.Site_ID(i))==string(meta.Site_ID(i+2));
        next_i = 2;
    elseif string(meta.Site_ID(i))==string(meta.Site_ID(i+1));
        next_i = 1;
    else 
        display('not printing melt')
        continue
    end

    depth = meta.Approx_Depth(i);
    cfg.bedSearchRange = [depth-0.1 depth+0.1];
    cfg.maxDepthConfig = depth-30;
    cfg.maxRange = depth+30;
    f = string(meta.File_Name(i));
    g = string(meta.File_Name(i+next_i));
    data_out = fmcw_melt2(f,g);
    data.meltRate(i) = data_out.meltRate;     
    %print(ya.meltRate)
    display('printed melt')
        
    

end


f = strcat("/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/",t6)
g =  strcat("/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/",g6)
% 
% fmcw_plot(tt1,'maxrange',800)
% fmcw_plot('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/2019-12-21_215758.dat','maxrange',800)





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
