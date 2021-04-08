#!/bin/sh
cd /mnt/usersVol2/sulc/MM_refitting/for_petr/CUDADNA/SIM_normalDNA/BULGE_1/RUN_0  
echo =========================================================   
echo PBS job: submitted  date = Wed Feb 20 18:52:11 GMT 2013      
date_start=`date +%s`
echo 1 cpus on 1 nodes \( SMP=1 \)        
#echo Executable file: /mnt/usersVol2/sulc/MM_refitting/for_petr/CUDADNA/SIM_normalDNA/BULGE_1/RUN_0/../../../Release/oxDNA                             
#echo $PBS_NODEFILE
#cat $PBS_NODEFILE
#ls -l $PBS_NODEFILE
echo 1 hosts used:                                 
echo -------------                    
cat $PBS_NODEFILE | cut -f 1 -d \. | sort  | fmt -w 30
echo =========================================================   
#cat $PBS_NODEFILE
echo Job output begins                                           
echo -----------------                                           
echo   
/mnt/usersVol2/sulc/MM_refitting/for_petr/CUDADNA/SIM_normalDNA/BULGE_1/RUN_0/../../../Release/oxDNA input
echo   
echo ---------------                                           
echo Job output ends                                           
date_end=`date +%s`
seconds=$((date_end-date_start))
minutes=$((seconds/60))
seconds=$((seconds-60*minutes))
hours=$((minutes/60))
minutes=$((minutes-60*hours))
echo =========================================================   
echo PBS job: finished   date = `date`   
echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
echo =========================================================


#echo Clearing up any stray processes...
#for oi in `cat $PBS_NODEFILE`
#do
#    echo $oi
#    rsh $oi /usr/local/bin/clearup
#done
#echo Done

# clean up machine file   
#rm -f $HOME/.gmpi/conf.$PBS_JOBID
