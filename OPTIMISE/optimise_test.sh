main_path=$(pwd)
echo "Usage: "
echo "bash optimise.sh config_file"

if [ "$#" -ne 1 ]; then
	echo "Illegal number of parameters"
	exit 1
fi

config=$1 #configuration file

#parse input from conf file
oxDNA_path=$(awk '$1 == "OXDNA_PATH" {print $2}' $config)
if [[ "$oxDNA_path" == "" ]]; then
	echo "Ox DNA path not specified. Usage"
	echo "In config_file: "
	echo "OXDNA_PATH path"
	echo "Terminating."
	exit 1
fi
echo "oxDNA path:"
echo $oxDNA_path
opti_path=$(awk '$1 == "OPTI_PATH" {print $2}' $config)
if [[ "$opti_path" == "" ]]; then
	echo "Optimisation program path not specified. Usage"
	echo "In config_file: "
	echo "OPTI_PATH path"
	echo "Terminating."
	exit 1
fi
echo "Optimisation program path:"
echo $opti_path
box_size=$(awk '$1 == "BOX_SIZE" {print $2}' $config)
if [[ "$box_size" == "" ]]; then
	echo "Box size not specified. Setting it to 50. Usage:"
	echo "In config_file: "
	echo "BOX_SIZE box_size"
	box_size=50
fi
echo "Sim box size:"
echo $box_size
Nreps=$(awk '$1 == "REPS" {print $2}' $config)
if [[ "$Nreps" == "" ]]; then
	echo "Number of replicas not specified. Setting it to 1. Usage:"
	echo "In config_file: "
	echo "REPS nreplica"
	Nreps=1
fi
echo "# of replicas:"
echo $Nreps
Nsteps=$(awk '$1 == "NSTEPS" {print $2}' $config)
if [[ "$Nsteps" == "" ]]; then
	echo "Number of optimisation steps. Setting it to 20. Usage:"
	echo "In config_file: "
	echo "NSTEPS nsteps"
	Nsteps=20
fi
echo "# of optimisation steps:"
echo $Nsteps
timesteps=$(awk '$1 == "TIMESTEPS" {print $2}' $config)
if [[ "$timesteps" == "" ]]; then
	echo "Number of timesteps not specified. Setting it to 1e7. Usage:"
	echo "In config_file: "
	echo "TIMESTEPS timesteps"
	timesteps=1e7
fi
echo "# of sim timesteps:"
echo $timesteps

#other parameters

#initialise first step

#generate gen.txt file for oxDNA
seq=$(awk '$1 == "SEQ" {print $2}' $config)
if [[ "$seq" == "" ]]; then
	echo "Sequence not specified. Usage:"
	echo "In config_file: "
	echo "SEQ ACGTCGA..."
	echo "Terminating."
	exit 1
fi
echo "Sequence:"
echo $seq

echo "double "$seq > ${main_path}/gen.txt

#create directory
mkdir ${main_path}/Step0
cp ${main_path}/oxDNA_sequence_dependent_parameters.txt ${main_path}/Step0/
cp ${main_path}/oxDNA_sequence_dependent_parameters_in.txt ${main_path}/Step0/

#setup and run repetitions for step 0
for ((j=0; j< ${Nreps};j++)); do
	mkdir -p ${main_path}/Step0/Rep${j}/	
	cp ${main_path}/gen.txt ${main_path}/Step0/Rep${j}/
	cp ${opti_path}/input_MD ${main_path}/Step0/Rep${j}/
	cp ${main_path}/Step0/oxDNA_sequence_dependent_parameters_in.txt ${main_path}/Step0/Rep${j}/
	cp ${opti_path}/input1.an ${main_path}/Step0/Rep${j}/
	cp ${opti_path}/input2.an ${main_path}/Step0/Rep${j}/
	
	cd ${main_path}/Step0/
	rep_path=$(pwd)/Rep${j}
	step_path=$(pwd)
	cd ${rep_path}
	sed -i "s|replace_me1|${rep_path}|g" input1.an	
	sed -i "s|replace_me1|${rep_path}|g" input2.an
	sed -i "s|replace_me2|${step_path}|g" input1.an	
	sed -i "s|replace_me2|${step_path}|g" input2.an
	python3 ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt
	sed -i "s|seed = 1|seed = ${RANDOM}|g" input_MD
	sed -i "s|steps = 1e7|steps = ${timesteps}|g" input_MD
	${oxDNA_path}/build/bin/oxDNA input_MD > out &
done

wait

#optimise
cd ${main_path}/Step0
python3 ${opti_path}/optimise.py ../$config


mv oxDNA_sequence_dependent_parameters_tmp.txt oxDNA_sequence_dependent_parameters_fin.txt

for ((i=1; i < ${Nsteps}; i++)); do

	mkdir -p ${main_path}/Step${i}/
	
	#copy optimisation code
	cp ${main_path}/oxDNA_sequence_dependent_parameters.txt ${main_path}/Step${i}/
	cp ${main_path}/Step$(($i-1))/oxDNA_sequence_dependent_parameters_fin.txt ${main_path}/Step${i}/oxDNA_sequence_dependent_parameters_in.txt
	
	for ((j=0; j< ${Nreps};j++)); do
		mkdir -p ${main_path}/Step${i}/Rep${j}/
		cp ${main_path}/gen.txt ${main_path}/Step${i}/Rep${j}/
		cp ${opti_path}/input_MD ${main_path}/Step${i}/Rep${j}/
		cp ${main_path}/Step$(($i-1))/oxDNA_sequence_dependent_parameters_fin.txt ${main_path}/Step${i}/Rep${j}/
		cp ${opti_path}/input1.an ${main_path}/Step${i}/Rep${j}/
		cp ${opti_path}/input2.an ${main_path}/Step${i}/Rep${j}/
		cd ${main_path}/Step${i}/Rep${j}/
		mv oxDNA_sequence_dependent_parameters_fin.txt oxDNA_sequence_dependent_parameters_in.txt
		
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
		sed -i "s|steps = 1e7|steps = ${timesteps}|g" input_MD
		${oxDNA_path}/build/bin/oxDNA input_MD > out &
	done
	
	wait
	
	cd ${main_path}/Step${i}
	python3 ${opti_path}/optimise.py ../$config
	mv oxDNA_sequence_dependent_parameters_tmp.txt oxDNA_sequence_dependent_parameters_fin.txt
done
