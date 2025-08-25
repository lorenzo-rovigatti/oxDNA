path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

poly_steps=('AA' 'AG' 'AC' 'AT' 'CG' 'CC')
richness=(10 30 50 70 90)

Ncores=11

counts=0
for((i = 0; i < 11; i++)); do

	echo "Running rep: "${counts}
	counts=$(($counts+1))
	
	sim_path=${path}/${i}
	mkdir -p ${sim_path}
	cp oxDNA3_sequence_dependent_parameters.txt ${sim_path}
	cp input_MD ${sim_path}
	
	cd ${sim_path}
	
	sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD
	
	if [ $i -lt 6 ]; then    #these are all poly repeats
	        sequence=""
		for((j = 0; j < 50; j++)); do
			sequence+=${poly_steps[$i]}
		done
		echo "double ${sequence}" > seq.txt
	
	else    #these are random sequences with various richness
	        entry=$(($i-6))
		python3 ${oxdna_path}/OPTIMISE/Utils/gen_seqs_richness.py 100 ${richness[$entry]} $((100-${richness[$entry]})) $((100-${richness[$entry]})) ${richness[$entry]}
	
	fi
	
	python3 ${oxdna_path}/oxdna/utils/generate-sa.py 50 seq.txt
		
	${oxdna_path}/build/bin/oxDNA input_MD > out &
	
	if [ $counts -eq ${Ncores} ]; then
		wait
	fi
	
done
