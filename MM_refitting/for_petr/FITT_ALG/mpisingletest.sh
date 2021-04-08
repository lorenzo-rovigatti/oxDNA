#!/bin/sh
cd /usersVol2/sulc/MM_refitting/for_petr/FITT_ALG  
echo =========================================================   
echo PBS job: submitted  date = Sun Jul  6 22:56:27 BST 2014      
date_start=`date +%s`
echo 16 cpus on 16 nodes \( SMP=1 \)        
#echo Executable file: /usersVol2/sulc/MM_refitting/for_petr/FITT_ALG/SRC_TWO/mpisingletest                             
#echo $PBS_NODEFILE
#cat $PBS_NODEFILE
#ls -l $PBS_NODEFILE
echo 16 hosts used:                                 
echo -------------                    
cat $PBS_NODEFILE | cut -f 1 -d \. | sort  | fmt -w 30
echo =========================================================   
cat $PBS_NODEFILE >>/usersVol2/jop/jobnodes/$PBS_JOBID
echo Job output begins                                           
echo -----------------                                           
echo   
/usr/local/shared/mpich3hydraintel/bin/mpiexec -bootstrap=ssh -f $PBS_NODEFILE -n 16 /usersVol2/sulc/MM_refitting/for_petr/FITT_ALG/SRC_TWO/mpisingletest mpiMALY5gTEST new_rescaled.txt H 1 1 maly5g
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
