main_path=$(pwd)
oxDNA_path=/home/users/yqb22156/oxDNA/oxdna/
box_size=50

Nreps=4
Nsteps=50
#other parameters

#copy optimisation code
cp ${main_path}/target_mean_and_covariance.txt ${main_path}/Step0/
cp ${main_path}/codenames.txt ${main_path}/Step0/
cp ${main_path}/oxdna_to_internal_wflip.py ${main_path}/Step0/
cp ${main_path}/continuity_constraints.py ${main_path}/Step0/
cp ${main_path}/oxDNA2_sequence_dependent_parameters.txt ${main_path}/Step0/
#sed -i "s|Nreps = *|Nreps = ${Nreps}|g" ${main_path}/main.py
#sed other relevant parameters
cp ${main_path}/main.py ${main_path}/Step0/

#setup and run repetitions for step 0
for ((j=0; j< ${Nreps};j++)); do
	mkdir -p ${main_path}/Step0/Rep${j}/	
	cp ${main_path}/Step0/gen.txt ${main_path}/Step0/Rep${j}/
	cp ${main_path}/Step0/input_MD ${main_path}/Step0/Rep${j}/
	cp ${main_path}/Step0/oxDNA2_sequence_dependent_parameters_in.txt ${main_path}/Step0/Rep${j}/
	cp ${main_path}/input1.an ${main_path}/Step0/Rep${j}/
	cp ${main_path}/input2.an ${main_path}/Step0/Rep${j}/
	
	cd ${main_path}/Step0/
	rep_path=$(pwd)/Rep${j}
	step_path=$(pwd)
	cd ${rep_path}
	sed -i "s|replace_me1|${rep_path}|g" input1.an	
	sed -i "s|replace_me1|${rep_path}|g" input2.an
	sed -i "s|replace_me2|${step_path}|g" input1.an	
	sed -i "s|replace_me2|${step_path}|g" input2.an
	python3 ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt
	sed -i "s|seed = *|seed = ${RANDOM}|g" input_MD
	${oxDNA_path}/build/bin/oxDNA input_MD > out &
done

wait

#optimise (5 steps)
cd ${main_path}/Step0
python3 main.py

mv oxDNA2_sequence_dependent_parameters_tmp.txt oxDNA2_sequence_dependent_parameters_fin.txt

for ((i=1; i < ${Nsteps}; i++)); do

	mkdir -p ${main_path}/Step${i}/
	
	#copy optimisation code
	cp ${main_path}/target_mean_and_covariance.txt ${main_path}/Step${i}/
	cp ${main_path}/codenames.txt ${main_path}/Step${i}/
	cp ${main_path}/oxdna_to_internal_wflip.py ${main_path}/Step${i}/
	cp ${main_path}/continuity_constraints.py ${main_path}/Step${i}/
	cp ${main_path}/oxDNA2_sequence_dependent_parameters.txt ${main_path}/Step${i}/
	#sed -i "s|Nreps = *|Nreps = ${Nreps}|g" ${main_path}/main.py
	#sed other relevant parameters
	cp ${main_path}/main.py ${main_path}/Step${i}/
	cp ${main_path}/Step$(($i-1))/oxDNA2_sequence_dependent_parameters_fin.txt ${main_path}/Step${i}/oxDNA2_sequence_dependent_parameters_in.txt
	
	for ((j=0; j< ${Nreps};j++)); do
		mkdir -p ${main_path}/Step${i}/Rep${j}/
		cp ${main_path}/Step0/gen.txt ${main_path}/Step${i}/Rep${j}/
		cp ${main_path}/Step0/input_MD ${main_path}/Step${i}/Rep${j}/
		cp ${main_path}/Step$(($i-1))/oxDNA2_sequence_dependent_parameters_fin.txt ${main_path}/Step${i}/Rep${j}/
		cp ${main_path}/input1.an ${main_path}/Step${i}/Rep${j}/
		cp ${main_path}/input2.an ${main_path}/Step${i}/Rep${j}/
		cd ${main_path}/Step${i}/Rep${j}/
		mv oxDNA2_sequence_dependent_parameters_fin.txt oxDNA2_sequence_dependent_parameters_in.txt
		
		cd ${main_path}/Step${i}/
		rep_path=$(pwd)/Rep${j}
		step_path=$(pwd)
		cd ${rep_path}
		sed -i "s|replace_me1|${rep_path}|g" input1.an	
		sed -i "s|replace_me1|${rep_path}|g" input2.an
		sed -i "s|replace_me2|${step_path}|g" input1.an	
		sed -i "s|replace_me2|${step_path}|g" input2.an
		python3 ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt
		sed -i "s|seed = *|seed = ${RANDOM}|g" input_MD
		${oxDNA_path}/build/bin/oxDNA input_MD > out &
	done
	
	wait
	
	cd ${main_path}/Step${i}
	python3 ${main_path}/main.py
	mv oxDNA2_sequence_dependent_parameters_tmp.txt oxDNA2_sequence_dependent_parameters_fin.txt
done
