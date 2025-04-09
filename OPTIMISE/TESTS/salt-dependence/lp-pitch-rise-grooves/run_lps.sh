path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

Cs=(0.1 0.15 0.2 0.3 0.5 1.1)

Ncores=6
Nseqs=6
counts=0

 #generate here so that the sequence is the same for all salt concentrations
python3 ${oxdna_path}/OPTIMISE/Utils/gen_seqs_richness.py 100 1 1 1 1 #random sequence

for((i = 0; i < ${Nseqs}; i++)); do

	echo "Running rep: "${counts}
	counts=$(($counts+1))
	
	sim_path=${path}/${i}
	mkdir -p ${sim_path}
	cp oxDNA3_sequence_dependent_parameters.txt ${sim_path}
	cp input_MD ${sim_path}
	cp seq.txt ${sim_path}
	
	cd ${sim_path}
	
	sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD
	
	sed -i "s|salt_concentration = .*|salt_concentration = ${Cs[$i]}|g" input_MD
	
	python3 ${oxdna_path}/oxdna/utils/generate-sa.py 50 seq.txt
		
	${oxdna_path}/build/bin/oxDNA input_MD > out &
	
	if [ $counts -eq ${Ncores} ]; then
		wait
	fi
	
done
