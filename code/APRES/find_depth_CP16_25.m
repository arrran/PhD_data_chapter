addpath(genpath('/Volumes/arc_04/FIELD_DATA/K8621920/APRES/'))

cp16 = '2019-12-22_031028.dat'
cp17 = '2019-12-22_032236.dat'
cp18 = '2019-12-22_033556.dat'
cp19 = '2019-12-22_034921.dat'
cp20 = '2019-12-22_040202.dat'
cp21 = '2019-12-22_041156.dat'
cp22 = '2019-12-22_042312.dat'
cp23 = '2019-12-22_043429.dat'
cp24 = '2019-12-22_044431.dat'
cp25 = '2019-12-22_045518.dat'

f = strcat("/Volumes/arc_04/FIELD_DATA/K8621920/APRES/backup_apres_31dec/Survey/",cp25)

fmcw_plot(f,'maxrange',600)


