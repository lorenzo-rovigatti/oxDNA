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

#SBATCH --time 48:00:00 # 2 days, more than enough

#SBATCH --job-name=PTVMMC_s0r0_n15

path=/users/yqb22156/Sample_Melting/n15
path_oxDNA=/users/yqb22156/oxdna

l=0
m=0

path_rep=${path}/Seq${l}/Rep${m}
cd ${path_rep}
mpirun -n 8 ${path_oxDNA}/build_mpi/bin/oxDNA_mpi input_PT_VMMC > out
