%(* runstring ='matlab '+'-nodisplay '+'-nosplash '+"-r "+  '\"run(\'' +'/Users/home/whitefar/MODELLING/MITgcm/code/make_phi0surf.m'+ %"\');exit;\""
%        subprocess.call(runstring,  shell=True)
%        print('phi0surf written to '+self.model_input_folder+'phi0surf') *)


addpath(genpath('/Users/home/whitefar/DATA/code/APRES_TIMESERIES/fmcw_210118/'))
addpath(genpath('/Users/home/whitefar/DATA/code/APRES_TIMESERIES/'))
addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))
addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8622021/KIS Science/ApRES data files/KIS2'))

for i = 9:-1:1
    calculate_meltrate_fn(i)
    
end
disp('all done')

calculate_meltrate_fn(29)

for i = 1:31
    disp(i)
    
end