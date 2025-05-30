path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

poly_steps=('AA' 'AG' 'AC' 'AT' 'CG' 'CC' 'CT' 'GG' 'GT' 'TT')

Ncores=11
Nseqs=11

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
	
	if [ $i -lt 10 ]; then    #these are all poly repeats
	        sequence=""
		for((j = 0; j < 16; j++)); do
			sequence+=${poly_steps[$i]}
		done
		echo "${sequence}" > seq.txt #no double = single strand
	
	else    #these are random sequences with various richness
	        entry=$(($i-10))
		../GEN ${RANDOM}
		python3 ${oxdna_path}/oxdna/utils/generate-sa.py 20 seq.txt
	
	fi
	
	python3 ${oxdna_path}/oxdna/utils/generate-sa.py 20 seq.txt
		
	${oxdna_path}/build/bin/oxDNA input_MD > out &
	
	if [ $counts -eq ${Ncores} ]; then
		wait
	fi
	
done
