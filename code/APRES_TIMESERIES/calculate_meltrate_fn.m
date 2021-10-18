function calculate_meltrate_fn(N)

    path = strcat('/Users/home/whitefar/DATA/APRES_TIMESERIES/kis2_apres_timeseries_meta_',int2str(N),'days.csv');
    data = readtable(path);

    global cfg
    cfg.polyorder=1;
    cfg = fmcw_process_config_vsr;

    cfg.doPlotMelt = 0;
    cfg.doPlotAll = 0;
    cfg.doSaveOutput = 0;
    cfg.verbose = 0;

    %bulk calculations
    for i = 1:height(data)

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

    path_out = strcat('/Users/home/whitefar/DATA/APRES_TIMESERIES/apres_timeseries_output_',int2str(N),'days.csv');
    writetable(data,path_out)
end
