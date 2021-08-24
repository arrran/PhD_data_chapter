addpath(genpath('/DATA/code/APRES/fmcw_210118/'))
addpath(genpath('/DATA/code/APRES/'))
addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))
addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8622021/KIS Science/ApRES data files/KIS2'))
meta = readtable('/Users/home/whitefar/DATA/ApRES/kis2_apres_metadata_by_Site_ID.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%WITH MANUALLY ADJUSTED DEPTHS
data = meta;

data.meltRate = nan*zeros(height(meta),1);
data.meltRateError = nan*zeros(height(meta),1);
data.bed = nan*zeros(height(meta),1);
data.dhStrain= nan*zeros(height(meta),1);
data.dheStrain= nan*zeros(height(meta),1);
data.dheStrain= nan*zeros(height(meta),1);
data.dheStrain= nan*zeros(height(meta),1);
%data.output(1) = ya;

global cfg
cfg.polyorder=1;
cfg = fmcw_process_config_vsr;

cfg.doPlotMelt = 0;
cfg.doPlotAll = 0;
cfg.doSaveOutput = 0;


data.Approx_Depth(45) = 522;
data.Approx_Depth(2) = 711;
data.Approx_Depth(11) = 717.5;
data.Approx_Depth(14) = 614.6;
data.Approx_Depth(23) = 568.6;
data.Approx_Depth(39) = 548;
data.Approx_Depth(42) = 433;

%bulk calculations
for i = 2:height(meta)-1
    disp(i)
    
    if i>2
        if string(meta.Site_ID(i-2))==string(meta.Site_ID(i))
            if abs(meta.Epoch(i-2)-meta.Epoch(i))~=0
                next_i=2;
            end
        elseif string(meta.Site_ID(i-1))==string(meta.Site_ID(i))
            if abs(meta.Epoch(i-1)-meta.Epoch(i))~=0
                next_i=1;
            end
        else 
            disp('not printing melt');
            continue
        end
    else
        if string(meta.Site_ID(i-1))==string(meta.Site_ID(i))
            if abs(meta.Epoch(i-1)-meta.Epoch(i))~=0
                next_i=1;
            end
        else 
            disp('not printing melt');
            continue
        end
    end

    depth = data.Approx_Depth(i);
    cfg.bedSearchRange = [depth-0.1 depth+0.1];
    if i==15
        cfg.bedSearchRange = [614 616];
    end
    if i==18
        cfg.bedSearchRange = [453.75 454.75];
    end
    if i==27
        cfg.bedSearchRange = [714.25 715.5];
    end
    if i==47
        cfg.bedSearchRange = [453 454.6];
    end
    if i==49
        cfg.bedSearchRange = [476.4 477.8];
    end
    if i==57
        cfg.bedSearchRange = [472 474];
    end
    if i==61
        cfg.bedSearchRange = [451.5 453.5];
    end
    if i==63
        cfg.bedSearchRange = [471 473];
    end
    if i==65
        cfg.bedSearchRange = [431.5 434];
    end
    cfg.maxDepthConfig = depth-30;
    cfg.maxRange = depth+30;
    f = string(meta.File_Name(i-next_i));
    g = string(meta.File_Name(i));
    data_out = fmcw_melt2(f,g);
    
    
    data.meltRate(i) = data_out.meltRate;
    data.meltRateError(i) = data_out.meltRateError;
    data.bed(i) = data_out.bed.range;
    data.dhStrain(i) = data_out.bed.dhStrain;
    data.dheStrain(i) =data_out.bed.dheStrain;
    data.dheStrain(i) = data_out.bed.dheStrain;
    data.dheStrain(i) = data_out.bed.dheStrain;
    %print(ya.meltRate)
    disp('printed melt')
            
end


writetable(data,'/Users/home/whitefar/DATA/ApRES/kis2_meltrates_matlab.csv')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.Approx_Depth(43) = 522;
data.Approx_Depth(11) = 717.5;


data.Approx_Depth(64) = 432.5;
data.Approx_Depth(15) = 699.223;
%individual detailed search down channel.
i=15;

cfg.doPlotMelt = 1;
disp(i)
if string(meta.Site_ID(i-2))==string(meta.Site_ID(i))
    if abs(meta.Epoch(i-2)-meta.Epoch(i))~=0
        next_i=2;
    end
    dorun = 1; 
elseif string(meta.Site_ID(i-1))==string(meta.Site_ID(i))
    if abs(meta.Epoch(i-1)-meta.Epoch(i))~=0
        next_i=1;
    end
    dorun = 1;
else 
    disp('not printing melt');
    dorun = 0;    
end

if dorun
    depth = data.Approx_Depth(i);
    cfg.bedSearchRange = [depth-0.1 depth+0.1];
    if i==45
        cfg.bedSearchRange = [453 454.6];
    endd
    if i==47
        cfg.bedSearchRange = [476.4 477.8];
    end
    if i==55
        cfg.bedSearchRange = [472 474];
    end
    if i==59
        cfg.bedSearchRange = [451.5 453.5];
    end
    if i==61
        cfg.bedSearchRange = [471 473];
    end
    if i==63
        cfg.bedSearchRange = [431.5 434];
    end
    cfg.maxDepthConfig = depth-30;
    cfg.maxRange = depth+30;
    f = string(meta.File_Name(i-next_i));
    g = string(meta.File_Name(i));
    
    
%individual detailed search cross channel epoch 2.
i=15;

cfg.doPlotMelt = 1;
disp(i)
if string(meta.Site_ID(i-1))==string(meta.Site_ID(i))
    if abs(meta.Epoch(i-1)-meta.Epoch(i))~=0
        next_i=1;
    end
    dorun = 1;
else 
    disp('not printing melt');
    dorun = 0;    
end

if dorun
    depth = data.Approx_Depth(i);
    cfg.bedSearchRange = [depth-0.1 depth+0.1];
    if i==45
        cfg.bedSearchRange = [453 454.6];
    end

    cfg.maxDepthConfig = depth-30;
    cfg.maxRange = depth+30;
    f = string(meta.File_Name(i-next_i));
    g = string(meta.File_Name(i));
    data_out = fmcw_melt2(f,g);
end

    data_out = fmcw_melt2(f,g);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
data.Approx_Depth(2) = 711;
data.Approx_Depth(11) = 717.5;
data.Approx_Depth(14) = 614.6;
data.Approx_Depth(23) = 568.6;
data.Approx_Depth(41) = 433;
data.Approx_Depth(15) = 699.223;

%individual detailed search cross channel epoch 2.
i=15;

cfg.doPlotMelt = 1;
disp(i)
if string(meta.Site_ID(i-1))==string(meta.Site_ID(i))
    if abs(meta.Epoch(i-1)-meta.Epoch(i))~=0
        next_i=1;
    end
    dorun = 1;
else 
    disp('not printing melt');
    dorun = 0;    
end

if dorun
    depth = data.Approx_Depth(i);
    cfg.bedSearchRange = [depth-0.1 depth+0.1];
    if i==45
        cfg.bedSearchRange = [453 454.6];
    end

    cfg.maxDepthConfig = depth-30;
    cfg.maxRange = depth+30;
    f = string(meta.File_Name(i-next_i));
    g = string(meta.File_Name(i));
    data_out = fmcw_melt2(f,g);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%










data.meltRate(i) = data_out.meltRate;
data.meltRateError(i) = data_out.meltRateError;
data.bed(i) = data_out.bed.range;

%print(ya.meltRate)
disp('printed melt')

if 1
    disp('yes')
end



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
