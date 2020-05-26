addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))

test = '2019-12-07_225600.dat'
test = '2019-12-07_231044.dat'
test = '2019-12-07_232530.dat'
test = '2019-12-07_233943.dat'
test = '2019-12-07_235231.dat'
test = '2019-12-08_000409.dat'
test = '2019-12-08_001936.dat'
test = '2019-12-08_003236.dat'
test = '2019-12-08_013943.dat'
test = '2019-12-08_015742.dat'
test = '2019-12-08_021102.dat'
test ='2019-12-08_022603.dat'
test ='2019-12-08_024034.dat'
test ='2019-12-08_025312.dat'
test ='2019-12-08_030609.dat'
test ='2019-12-22_031028.dat'
test ='2019-12-22_032236.dat'
test ='2019-12-22_033556.dat'
test ='2019-12-22_034921.dat'
test ='2019-12-22_040202.dat'
test ='2019-12-22_041156.dat'
test ='2019-12-22_042312.dat'
test ='2019-12-22_043429.dat'
test ='2019-12-22_044431.dat'
test ='2019-12-22_045518.dat'

tt =['/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/',test]


fmcw_plot(tt,'maxrange',800)

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
