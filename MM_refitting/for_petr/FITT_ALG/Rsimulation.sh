#!/bin/sh
cd /usersVol2/sulc/MM_refitting/for_petr/FITT_ALG  
echo =========================================================   
echo PBS job: submitted  date = Tue Mar 12 18:41:45 GMT 2013      
date_start=`date +%s`
echo 20 cpus on 20 nodes \( SMP=1 \)        
#echo Executable file: /usersVol2/sulc/MM_refitting/for_petr/FITT_ALG/SRC/Rsimulation                             
#echo $PBS_NODEFILE
#cat $PBS_NODEFILE
#ls -l $PBS_NODEFILE
echo 20 hosts used:                                 
echo -------------                    
cat $PBS_NODEFILE | cut -f 1 -d \. | sort  | fmt -w 30
echo =========================================================   
#cat $PBS_NODEFILE
echo Job output begins                                           
echo -----------------                                           
echo   
/usr/local/shared/mpich2hydra/bin/mpiexec -bootstrap=ssh -f $PBS_NODEFILE -n 20 /usersVol2/sulc/MM_refitting/for_petr/FITT_ALG/SRC/Rsimulation sekvence_fit SEKVENCE/old_param.last HST 1 1000000 SEKVENCE/seq_6_all_HTg0.txt.CORR ../data/6.FIT SEKVENCE/seq_8_2000_HT.txt.CORR ../data/8.FIT SEKVENCE/seq_10_2000_HT.txt.CORR ../data/10.FIT SEKVENCE/seq_12_2000_HT.txt.CORR ../data/12.FIT SEKVENCE/seq_18_HT_2000.txt.CORR ../data/18.FIT
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
