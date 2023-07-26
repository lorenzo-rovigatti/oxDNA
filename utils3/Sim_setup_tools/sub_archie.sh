#!/bin/bash -l
#SBATCH -export=ALL

#SBATCH --partition=standard
#
#Specify project account
#SBATCH --account=henrich-oxdna3
#

#SBATCH --nodes=1
#SBATCH --ntasks=4	#max 40 per node
#SBATCH --cpus-per-task=1

#SBATCH --time 48:00:00 # 12 hours

#SBATCH --job-name=jname

#SBATCH --mail-type=END
#SBATCH --mail-user=andrea.bonato@strath.ac.uk
#

module load gcc python/2.7



