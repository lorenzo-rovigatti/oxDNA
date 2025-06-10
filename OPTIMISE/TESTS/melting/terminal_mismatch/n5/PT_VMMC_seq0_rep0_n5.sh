#!/bin/bash -l
#SBATCH -export=ALL

#SBATCH --partition=standard
#
#Specify project account
#SBATCH --account=henrich-oxdna3
#

#SBATCH --nodes=1
#SBATCH --ntasks=8     #max 40 per node
#SBATCH --cpus-per-task=1

#SBATCH --time 24:00:00 # 1 day, more than enough

#SBATCH --job-name=PTVMMC_s0r0_n5

path=/users/yqb22156/Sample_Melting/n5
path_oxDNA=/users/yqb22156/oxdna

l=0
m=0

path_rep=${path}/Seq${l}/Rep${m}
cd ${path_rep}
mpirun -n 8 ${path_oxDNA}/build_mpi/bin/oxDNA_mpi input_PT_VMMC > out
