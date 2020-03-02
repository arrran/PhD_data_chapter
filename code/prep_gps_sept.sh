#!/bin/sh
# 
# Preprocessing SEPTENTRIO GPS data 
# Huw Horgan 
if [ $# -ne 3 ];then 
echo "# USAGE: prep_gps_auto.sh Site AntH expt"
echo "# RUN FROM expt/rinex/ directory"
echo "#"
echo "#	Inputs	Site		gps site name (e.g. arc01)"
echo "#		AntH		Antennae Height (m)"
echo "#         expt            experimets (e.g. ta04)"
exit 0
fi

site=$1
anth=$2
expt=$3

site=`echo $site | tr '[:lower:]'  '[:upper:]'`
echo "Processing $site" 

dir=`pwd`
cdir=`basename $dir`
echo "Current directory $cdir"
if [ $cdir = "rinex" ]; then
    int=30 # GPS sampling interval
    echo "Resampling to $int seconds"
    echo "Inserting antennae height $anth into obs header" 

    list=tmp.lis
    ls ./POLAR* > $list
    initposfile=./arc8_rxposition.txt
    n=1
    max=`wc -l $list | awk '{print $1}'`

    while (test $n -le $max)
        do
            infile=`awk '{if(NR==n) print($0)}' n=$n $list`
            echo $infile
	        if [ -f $infile ]; then             # check file exists
        		outdir=./
		        fname=`basename $infile .16_` # Year SPECIFIC
        		nname=arc8
                NNAME=`echo $nname | awk '{print toupper($0)}'`		        
                newname=`echo $nname | awk '{print tolower($0)}'`
        		yyyy=2016
		        yy=`echo $yyyy | cut -c3-4`
                jd=`echo $infile | cut -c21-23`
                mm=`doy $yy $jd | grep Date | awk '{print(substr($0,11,2))}'`
                dd=`doy $yy $jd | grep Date | awk '{print(substr($0,14,2))}'`
        		yyo=`echo $yy``echo o`
		        yyn=`echo $yy``echo n`
        		jd0=`echo $jd``echo 0`
		        gweek=`doy $yyyy $mm $dd | grep week | awk '{print($3)}'`
                echo $yyyy $mm $dd $yy $jd $yyo $yyn $gweek
                echo "------------------------------------"
# NOTE. PolaRx does not need runpkr conversion.
# BUT DOES NEED INTIAL POSITIONS CALCULATED
# ../../../TASMAN_KINEMATIC/bin/initial_pos.sh
#initposfile=../initpos_summary.txt
#if [ -f $initposfile ]; then
#    ix=`grep "${nname}\ ${yyyy}\ ${mm}\ ${dd}" $initposfile | awk '{print($5)}'` 
#    iy=`grep "${nname}\ ${yyyy}\ ${mm}\ ${dd}" $initposfile | awk '{print($6)}'`
#    iz=`grep "${nname}\ ${yyyy}\ ${mm}\ ${dd}" $initposfile | awk '{print($7)}'`
#    echo $ix $iy $iz
#else 
#    echo "NEED TO POPULATE INTIAL POSITION FIRST"
#    echo "RUN ~/bin/gps_initial_pos.sh"
#    exit 0
#fi	
	# Create obs file names 
        		obsfile=${outdir}$newname$jd0.$yyo
		        navfile=${outdir}$newname$jd0.$yyn
                qcfile=${outdir}${fname}.qc
		        echo $obsfile
		        echo $navfile
		        echo $qcfile
                echo "------------------------------------"
		        echo "converting to rinex"
                # THIS IS MESSY. FIRST RUN NEXT LINE TO CREATE NAV FILES. 
	            # teqc -sep sbf -R -S -E -C -J -I -week $gweek -O.int $int -O.dec 30s -O.at SEPPOLANT_X_MF -O.pe $anth 0 0 +nav $navfile $infile 
                # NOW CALCULATE INTIAL POSITION
                # RUN ~/bin/gps_initial_pos.sh
                ix=`grep "${NNAME}\ ${yyyy}\ ${mm}\ ${dd}" $initposfile | awk '{print($5)}'` 
                iy=`grep "${NNAME}\ ${yyyy}\ ${mm}\ ${dd}" $initposfile | awk '{print($6)}'`
                iz=`grep "${NNAME}\ ${yyyy}\ ${mm}\ ${dd}" $initposfile | awk '{print($7)}'`
                echo $ix $iy $iz
                teqc -sep sbf -E -S -C -J -I -week $gweek -O.int $int -O.dec 30s -O.at SEPPOLANT_X_MF \
                -O.px $ix $iy $iz \
                -O.rn 4501365 \
                -O.rt "SEPT POLARX5" \
                -O.mn "ARC8" \
                -O.rv 1.0 \
                -O.an SEPPOLANT_X_MF \
                -O.pe $anth 0 0 $infile > $obsfile
			    teqc +qc -report $obsfile > $qcfile
	fi	
        n=`expr $n + 1`
    done
else
	echo "Run from rinex directory only. Exiting"
	exit 0
fi
 


