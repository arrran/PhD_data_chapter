# PreProcessing Radar

Two python commands will preprocess the radar data -  cutting and timesyncing
1. exec(open('/Users/home/whitefar/DATA/code/RADAR_PROCESSING/preprocess_KIS1_lines.py').read())
2. exec(open('/Users/home/whitefar/DATA/code/RADAR_PROCESSING/preprocess_KIS2_lines.py').read())

## path must include '/home/whitefar/DATA/code' for load_ppp.py

## This will:

- Load process_radar. This has code to load and preprocess raw pixie files, with radarsurvey and radarline object classes.

- Preprocess each survey (on-off of the pixie) line, clipping, timesyncing etc. This is done in two files; preprocess_KIS1_lines.py and preprocess_KIS2_lines.py


GNSS
1. use code/septemtrio_to_rinex on raw files.
2. Upload to https://webapp.geod.nrcan.gc.ca/geod/tools-outils/ppp.php to get processed precise point p...?
3. process_radar uses load_ppp.py to load the correct day for syncing location with radar traces. 