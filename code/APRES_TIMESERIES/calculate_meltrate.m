%(* runstring ='matlab '+'-nodisplay '+'-nosplash '+"-r "+  '\"run(\'' +'/Users/home/whitefar/MODELLING/MITgcm/code/make_phi0surf.m'+ %"\');exit;\""
%        subprocess.call(runstring,  shell=True)
%        print('phi0surf written to '+self.model_input_folder+'phi0surf') *)


addpath(genpath('/Users/home/whitefar/DATA/code/APRES_TIMESERIES/fmcw_210118/'))
addpath(genpath('/Users/home/whitefar/DATA/code/APRES_TIMESERIES/'))
addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))
addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8622021/KIS Science/ApRES data files/KIS2'))
meta = readtable('/Users/home/whitefar/DATA/APRES_TIMESERIES/kis2_apres_timeseries_meta_15days.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = meta;

global cfg
cfg.polyorder=1;
cfg = fmcw_process_config_vsr;

cfg.doPlotMelt = 0;
cfg.doPlotAll = 0;
cfg.doSaveOutput = 0;

%bulk calculations
for i = 1:height(meta)
    
    depth = data.Approx_Depth(i);
    cfg.bedSearchRange = [depth-0.1 depth+0.1];
    cfg.maxDepthConfig = depth-30;
    cfg.maxRange = depth+30;
    f = string(data.File_Name1(i));
    g = string(data.File_Name2(i));
    data_out = fmcw_melt2(f,g);
        
    data.meltRate(i) = data_out.meltRate;
    data.meltRateError(i) = data_out.meltRateError;
    data.Approx_Depth(i) = data_out.bed.range;
    data.Approx_Depth(i+1) = data_out.bed.range+data_out.bed.dh; %set the next guess to the first shot
    data.dhStrain(i) = data_out.bed.dhStrain;
    data.dheStrain(i) =data_out.bed.dheStrain;
    data.dheStrain(i) = data_out.bed.dheStrain;
    data.dheStrain(i) = data_out.bed.dheStrain;            
end


writetable(data,'/Users/home/whitefar/DATA/APRES_TIMESERIES/meltrates_output.csv')
