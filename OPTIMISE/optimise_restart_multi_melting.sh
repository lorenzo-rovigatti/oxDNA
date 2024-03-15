main_path=$(pwd)
echo "Usage: "
echo "bash optimise.sh config_file restart_step"
echo "Note: restart_step is the first step of the restart run, not the last step of the interrupted run."

if [ "$#" -ne 2 ]; then
	echo "Illegal number of parameters"
	exit 1
fi

config=$1 #configuration file

rstep=$2
echo "Restarting from step: ${rstep}"

###########################
#INPUT FROM FILE
###########################

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

echo "Single strand concentration = $concentration mM"

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

#read seq

echo "Usage: "
echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"

st=$(awk '$1 == "SEQ" {print $2}' $config)
if [[ "$st" == "" ]]; then
	echo "Sequence not specified. Usage:"
	echo "In config_file: "
	echo "SEQ ACGTCGA..."
	echo "Terminating."
	exit 1
fi

seq=()
for word in $st; do
	seq+=("$word")
done
#seq=$(echo $st | tr " " "\n")

Nseq=${#seq[@]}

echo "Sequences:"

for (( k = 0; k < ${Nseq}; k++ )); do
	echo ${seq[$k]}
done
#read box_size

st=$(awk '$1 == "SEQ" {print $3}' $config)
if [[ "$st" == "" ]]; then
	echo "Box_size? Usage:"
	echo "In config_file: "
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
	exit 1
fi

box_size=()
for word in $st; do
	box_size+=("$word")
done

if [[ ${#box_size[@]} -eq ${Nseq} ]]; then 
	echo "Number check ok."
else
	echo "number of box sizes is not equal to number of sequences."
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
exit 1
fi
#seq=$(echo $st | tr " " "\n")

echo "Box sizes:"

for (( k = 0; k < ${Nseq}; k++ )); do
	echo ${box_size[$k]}
done

#Salt concentrations

st=$(awk '$1 == "SEQ" {print $4}' $config)
if [[ "$st" == "" ]]; then
	echo "salt? Usage:"
	echo "In config_file: "
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
	exit 1
fi

salt=()
for word in $st; do
	salt+=("$word")
done

if [[ ${#salt[@]} -eq ${Nseq} ]]; then 
	echo "Number check ok."
else
	echo "number of salt concentrations is not equal to number of sequences."
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
exit 1
fi


#seq=$(echo $st | tr " " "\n")

echo "Salt concentrations:"

for (( k = 0; k < ${Nseq}; k++ )); do
	echo ${salt[$k]}
done

#Temperatures

st=$(awk '$1 == "SEQ" {print $5}' $config)
if [[ "$st" == "" ]]; then
	echo "Temperature? Usage:"
	echo "In config_file: "
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
	exit 1
fi

Ts=()
for word in $st; do
	Ts+=("$word")
done
#seq=$(echo $st | tr " " "\n")

echo "Temperatures:"

for (( k = 0; k < ${Nseq}; k++ )); do
	echo ${Ts[$k]}
done

#order parameter files

st=$(awk '$1 == "SEQ" {print $8}' $config)
if [[ "$st" == "" ]]; then
	echo "Order parameter file? Usage:"
	echo "In config_file: "
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
	exit 1
fi

OP_files=()
for word in $st; do
	OP_files+=("$word")
done
#seq=$(echo $st | tr " " "\n")

if [[ ${#Ts[@]} -eq ${Nseq} ]]; then 
	echo "Number check ok."
else
	echo "number of OP files is not equal to number of sequences."
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
exit 1
fi

echo "OP files:"

for (( k = 0; k < ${Nseq}; k++ )); do
	echo ${OP_files[$k]}
done


# Weights files

st=$(awk '$1 == "SEQ" {print $9}' $config)
if [[ "$st" == "" ]]; then
	echo "Weights? Usage:"
	echo "In config_file: "
	echo "SEQ seq box_size salt_conc simT delta_extrapolated_T op_file weight_file"
	echo "Terminating."
	exit 1
fi

W_files=()
for word in $st; do
	W_files+=("$word")
done
#seq=$(echo $st | tr " " "\n")

if [[ ${#W_files[@]} -eq ${Nseq} ]]; then 
	echo "Number check ok."
else
	echo "number of weight files is not equal to number of sequences."
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
exit 1
fi

echo "W_files:"

for (( k = 0; k < ${Nseq}; k++ )); do
	echo ${W_files[$k]}
done


###########################
#OPTIMISATION
###########################


echo "Running ${Nsteps}-${rstep} optimisation steps."

#optimise-restart
cd ${main_path}/Step$(($rstep-1))


file_tmp=oxDNA_sequence_dependent_parameters_tmp.txt
file_fin=oxDNA_sequence_dependent_parameters_fin.txt

if [ ! -f "$file_fin" ]; then
	if [ ! -f "$file_tmp" ]; then
		mv oxDNA_sequence_dependent_parameters_tmp.txt oxDNA_sequence_dependent_parameters_fin.txt
	else
		echo "could not find $file_tmp nor $file_fin in ${main_path}/Step$(($rstep-1))"
		echo "cannot restart"
		exit 1
	fi
fi

for ((i=${rstep}; i < ${Nsteps}; i++)); do

	Tfile=${main_path}/Step$(($i-1))/final_rew_mTs.txt

	st=$(awk '{print $1}' $Tfile)
	if [[ "$st" == "" ]]; then
		echo "Rew Temperatures previous step missing. Optimisation not working? "
		echo "Terminating."
		exit 1
	fi

	rew_Ts=()
	for word in $st; do
		rew_Ts+=("$word")
	done
	
	echo "Rew temperatures:"

	for (( k = 0; k < ${Nseq}; k++ )); do
		echo ${rew_Ts[$k]}
	done

	mkdir -p ${main_path}/Step${i}/
	
	#copy optimisation code
	cp ${main_path}/oxDNA_sequence_dependent_parameters.txt ${main_path}/Step${i}/
	cp ${main_path}/Step$(($i-1))/oxDNA_sequence_dependent_parameters_fin.txt ${main_path}/Step${i}/oxDNA_sequence_dependent_parameters_in.txt
	for ((l=0; l < ${Nseq}; l++)); do
		for ((j=0; j< ${Nreps};j++)); do
			mkdir -p ${main_path}/Step${i}/Seq${l}/Rep${j}/
			cp ${main_path}/Step$(($i-1))/Seq${l}/Rep${j}/input_VMMC ${main_path}/Step${i}/Seq${l}/Rep${j}/
			cp ${main_path}/Step$(($i-1))/oxDNA_sequence_dependent_parameters_fin.txt ${main_path}/Step${i}/Seq${l}/Rep${j}/
			cp ${main_path}/Step$(($i-1))/Seq${l}/Rep${j}/${OP_files[$l]} ${main_path}/Step${i}/Seq${l}/Rep${j}/
			cp ${main_path}/Step$(($i-1))/Seq${l}/Rep${j}/${W_files[$l]} ${main_path}/Step${i}/Seq${l}/Rep${j}/
			cp ${opti_path}/input1_melting.an ${main_path}/Step${i}/Seq${l}/Rep${j}/
			cp ${opti_path}/input2_melting.an ${main_path}/Step${i}/Seq${l}/Rep${j}/
			cd ${main_path}/Step${i}/Seq${l}/Rep${j}/
			mv oxDNA_sequence_dependent_parameters_fin.txt oxDNA_sequence_dependent_parameters_in.txt
			
			cd ${main_path}/Step${i}/
			rep_path=$(pwd)/Seq${l}/Rep${j}
			step_path=$(pwd)
			cd ${rep_path}
			
			echo "double "${seq[$l]} > gen.txt
			
			sed -i "s|replace_me1|${rep_path}|g" input1_melting.an	
			sed -i "s|replace_me1|${rep_path}|g" input2_melting.an
			sed -i "s|replace_me2|${step_path}|g" input1_melting.an	
			sed -i "s|replace_me2|${step_path}|g" input2_melting.an
			sed -i "s|op.txt|${OP_files[$l]}|g" input1_melting.an
			sed -i "s|op.txt|${OP_files[$l]}|g" input2_melting.an
			python3 ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt
			sed -i "s|seed = *|seed = ${RANDOM}|g" input_VMMC
			sed -i "s|T = *|T = ${rew_Ts[$l]}C|g" input_VMMC
			sed -i "s|T = *|T = ${rew_Ts[$l]}C|g" input1_melting.an
			${oxDNA_path}/build/bin/oxDNA input_VMMC > out &
		done
	done
	
	wait
	
	cd ${main_path}/Step${i}
	python3 ${opti_path}/optimise_melting.py ../$config > OutOpti
	mv oxDNA_sequence_dependent_parameters_tmp.txt oxDNA_sequence_dependent_parameters_fin.txt
done
