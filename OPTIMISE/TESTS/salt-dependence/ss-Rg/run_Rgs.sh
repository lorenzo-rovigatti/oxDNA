path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

Ns=(20 40 60 80 100)
Cs=(0.1 0.15 0.2 0.3 0.5 1.1)

Nlens=${#Ns[@]}
NCs=${#Cs[@]}
Ncores=12

counts=0

for((i = 0; i < ${NCs}; i++)); do
	for((j = 0; j < ${Nlens}; j++)); do
		echo "Running rep: "${counts}
		counts=$(($counts+1))
	
		sim_path=${path}/${i}/${j}
		mkdir -p ${sim_path}
		cp oxDNA3_sequence_dependent_parameters.txt ${sim_path}
		cp input_MD ${sim_path}
	
		cd ${sim_path}
	
		sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD
		sed -i "s|salt_concentration = .*|salt_concentration = ${Cs[$i]}|g" input_MD
	
		python3 ${oxdna_path}/OPTIMISE/Utils/gen_seqs_richness_ss.py ${Ns[${j}]} 1 1 1 1 #random sequence
	
		python3 ${oxdna_path}/oxdna/utils/generate-sa.py 50 seq.txt
		
		${oxdna_path}/build/bin/oxDNA input_MD > out &
	
		if [ $counts -eq ${Ncores} ]; then
			wait
		fi
	done
done

wait
