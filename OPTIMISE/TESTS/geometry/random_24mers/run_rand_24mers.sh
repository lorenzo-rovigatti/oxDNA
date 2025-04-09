#this script sets up and runs MDs of Nseqs random 24mers

path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

Nseqs=10
Ncores=10

counts=0
for((i = 0; i < ${Nseqs}; i++)); do

	echo "Running rep: "${counts}
	counts=$(($counts+1))
	
	sim_path=${path}/${i}
	mkdir -p ${sim_path}
	cp oxDNA3_sequence_dependent_parameters.txt ${sim_path}
	cp input_MD ${sim_path}
	
	cd ${sim_path}
	
	sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD

	python3 ${oxdna_path}/OPTIMISE/Utils/gen_seqs_richness.py 24 1 1 1 1 #random 24 mers
	
	python3 ${oxdna_path}/oxdna/utils/generate-sa.py 20 seq.txt
		
	${oxdna_path}/build/bin/oxDNA input_MD > out &
	
	if [ $counts -eq ${Ncores} ]; then
		wait
	fi
	
done
