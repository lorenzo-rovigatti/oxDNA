path=$(pwd)

counts=0
for((i = 0; i < 5; i++)); do
	for((j = 0; j < 5; j++)); do
	if [ $i -eq 0 ] && [ $j -eq 0 ]; then
		continue
	fi
	
	echo ${counts}
	counts=$(($counts+1))
	
	sim_path=${path}/${i}${j}
	mkdir -p ${sim_path}
	cp oxDNA3_sequence_dependent_parameters.txt ${sim_path}
	cp input_MD ${sim_path}
	
	cd ${sim_path}
	
	sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD
	
	python3 ~/oxdna/OPTIMISE/Utils/gen_seqs_richness.py 100 ${i} ${j} ${j} ${i}
	python3 ~/oxdna/utils/generate-sa.py 50 seq.txt
	
	~/oxdna/build/bin/oxDNA input_MD > out &
	
	if [ $counts -eq 12 ]; then
		wait
	fi
	
	done
done
