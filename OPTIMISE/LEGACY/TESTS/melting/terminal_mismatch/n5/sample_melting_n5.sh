#!/bin/bash -l
#SBATCH -export=ALL

#SBATCH --partition=standard
#
#Specify project account
#SBATCH --account=henrich-oxdna3
#

#SBATCH --nodes=1
#SBATCH --ntasks=10     #max 40 per node
#SBATCH --cpus-per-task=1

#SBATCH --time 120:00:00 # 5 days

#SBATCH --job-name=melting_samp_n5

#SBATCH --mail-type=END
#SBATCH --mail-user=andrea.bonato@strath.ac.uk


path=/users/yqb22156/Sample_Melting/n5
oxDNA_path=/users/yqb22156/oxdna

nreps=1
box_size=20

seqs=("AACCT" "GAGGA" "GTGTC" "GTCGA" "GGTGG" "GACGG" "ACCCG" "CGCCA" "GCGCA" "GCGCG" )
OP_files=( "op_n5.txt" "op_n5.txt" "op_n5.txt" "op_n5.txt" "op_n5.txt" "op_n5.txt" "op_n5.txt" "op_n5.txt" "op_n5.txt" "op_n5.txt" )
W_files=( "wfile_n5.txt" "wfile_n5.txt" "wfile_n5.txt" "wfile_n5.txt" "wfile_n5.txt" "wfile_n5.txt" "wfile_n5.txt" "wfile_n5.txt" "wfile_n5.txt" "wfile_n5.txt"  )
Nseqs=${#seqs[@]}
nbps=( 5 5 5 5 5 5 5 5 5 5 )
mTs=( 9 12 15 19 22 25 29 32 35 42 )

for ((l=0; l < ${Nseqs}; l++)); do
	for ((m=0; m < ${nreps}; m++)); do

		rep_path=${path}/Seq${l}/Rep${m}

		mkdir -p ${rep_path}

		cp ${path}/input_VMMC ${rep_path}
		#cp ${path}/input_VMMC_relax ${rep_path}
		cp ${path}/${OP_files[$l]} ${rep_path}
		cp ${path}/${W_files[$l]} ${rep_path}
		cp ${path}/oxDNA_sequence_dependent_parameters_in.txt ${rep_path}
		cd ${rep_path}
		echo "double ${seqs[l]}" > gen.txt
		python3 ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt
		sed -i "s|seed = .*|seed = ${RANDOM}|g" input_VMMC
		sed -i "s|T = .*|T = ${mTs[$l]}C|g" input_VMMC
		sed -i "s|topology = .*|topology = ${rep_path}/generated.top|g" input_VMMC
		sed -i "s|conf_file = .*|conf_file = ${rep_path}/generated.dat|g" input_VMMC
		sed -i "s|trajectory_file = .*|trajectory_file = ${rep_path}/trajectory.dat|g" input_VMMC
		sed -i "s|name = .*|name = ${rep_path}/split_energy.dat|g" input_VMMC
		#sed -i "s|seed = .*|seed = ${RANDOM}|g" input_VMMC_relax
		#sed -i "s|T = .*|T = ${mTs[$l]}C|g" input_VMMC_relax
		extr_Ts="extrapolate_hist ="
		if [ $l -gt 0 ]; then
			for ((n=-5; n<=5; n++)); do
				ex_T=$(( $n*2 + ${mTs[$l]} ))
				if [ $n -lt 5 ]; then
					extr_Ts="${extr_Ts} ${ex_T}.0C,"
				else
					extr_Ts="${extr_Ts} ${ex_T}.0C"
				fi
			done
		else
			for ((n=-3; n<=8; n++)); do
				ex_T=$(( $n*2 + ${mTs[$l]} ))
				if [ $n -lt 5 ]; then
					extr_Ts="${extr_Ts} ${ex_T}.0C,"
				else
					extr_Ts="${extr_Ts} ${ex_T}.0C"
				fi
			done
		fi
		sed -i "s|extrapolate_hist = .*|${extr_Ts}|g" input_VMMC
		#${oxDNA_path}/build/bin/oxDNA input_VMMC_relax > out &

	done
done

wait

for ((l=0; l < ${Nseqs}; l++)); do
	for ((m=0; m < ${nreps}; m++)); do

		rep_path=${path}/Seq${l}/Rep${m}
		cd ${rep_path}
		${oxDNA_path}/build/bin/oxDNA input_VMMC > out &
	done
done

wait
