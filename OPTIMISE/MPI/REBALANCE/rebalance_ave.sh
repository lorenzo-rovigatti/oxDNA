main_path=$(pwd)
echo "Usage: "
echo "bash optimise.sh config_file"

if [ "$#" -ne 1 ]; then
	echo "Illegal number of parameters"
	exit 1
fi

config=$1 #configuration file

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
Nbp=()
for word in $st; do
	seq+=("$word")
	Nbp+=(${#word})
done
#seq=$(echo $st | tr " " "\n")

Nseq=${#seq[@]}

echo "Sequences:"

for (( k = 0; k < ${Nseq}; k++ )); do
	echo ${seq[$k]}
	echo "Nbp: " ${Nbp[$k]}
done

#total number of cpus is number_of_sequences*number_of_reps
Nproc=$((${Nreps}*${Nseq}))

echo "Number of CPUs used: ${Nproc}"


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
#REBALANCE
###########################

#setup files

CRST=( 76 78.375 80.75 83.125 85.5 87.875 90.25 92.625 95 )

Ncrst=${#CRST[@]}

for ((l=0; l < ${Nseq}; l++)); do
	for ((j=0; j< ${Nreps};j++)); do

		cp ${opti_path}/MPI/input1_melting.an ${main_path}/Nbp${Nbp[l]}/Rep${j}/
		cp ${opti_path}/MPI/input2_melting.an ${main_path}/Nbp${Nbp[l]}/Rep${j}/
		
		
		rep_path=${main_path}/Nbp${Nbp[l]}/Rep${j}/
		step_path=${main_path}
		
		cd ${main_path}/Nbp${Nbp[l]}/Rep${j}/

		sed -i "s|replace_me1|${rep_path}|g" input1_melting.an	
		sed -i "s|replace_me1|${rep_path}|g" input2_melting.an
		sed -i "s|replace_me2|${rep_path}|g" input1_melting.an	
		sed -i "s|replace_me2|${step_path}|g" input2_melting.an
		sed -i "s|op.txt|${OP_files[$l]}|g" input1_melting.an
		sed -i "s|op.txt|${OP_files[$l]}|g" input2_melting.an
		sed -i "s|T = .*|T = ${Ts[$l]}C|g" input1_melting.an
		sed -i "s|T = .*|T = ${Ts[$l]}C|g" input2_melting.an
		sed -i "s|salt_concentration = .*|salt_concentration = ${salt[$l]}|g" input1_melting.an
		sed -i "s|salt_concentration = .*|salt_concentration = ${salt[$l]}|g" input2_melting.an
	done
done

cd ${main_path}/


for ((k=0; k < ${Ncrst}; k++)); do

	if [ $k -ge 0 ]; then
		sed -i "s| = ${CRST[$(($k-1))]}| = ${CRST[$k]}|g" oxDNA_sequence_dependent_parameters.txt
	fi

	mpirun -np ${Nproc} python3 ${opti_path}/MPI/REBALANCE/rebalance_strengths_melting_mpi.py $config 0 > OutOpti

	mv oxDNA_sequence_dependent_parameters_in.txt oxDNA_sequence_dependent_parameters_in_copy.txt
	mv oxDNA_sequence_dependent_parameters.txt oxDNA_sequence_dependent_parameters_copy.txt
	mv oxDNA_sequence_dependent_parameters_tmp.txt oxDNA_sequence_dependent_parameters_in.txt
	head -34 oxDNA_sequence_dependent_parameters_copy.txt >> oxDNA_sequence_dependent_parameters.txt
	tail -n 16 oxDNA_sequence_dependent_parameters_in.txt >> oxDNA_sequence_dependent_parameters.txt

	
	mv parameters_v_iter.txt PARAMS/parameters_v_iter${k}_stck.txt
	mv final_rew_mTs.txt final_PARAMS/rew_mTs${k}_stck.txt
	mv OutOpti PARAMS/OutOpti${k}_stck
	
	mpirun -np ${Nproc} python3 ${opti_path}/MPI/REBALANCE/rebalance_strengths_melting_mpi.py opti_hydr.txt 0 > OutOpti
	
	mv oxDNA_sequence_dependent_parameters_tmp.txt PARAMS/oxDNA_sequence_dependent_parameters_fin${k}.txt
	
	mv oxDNA_sequence_dependent_parameters_copy.txt oxDNA_sequence_dependent_parameters.txt	
	mv oxDNA_sequence_dependent_parameters_in_copy.txt oxDNA_sequence_dependent_parameters_in.txt
	
	mv parameters_v_iter.txt PARAMS/parameters_v_iter${k}.txt
	mv final_rew_mTs.txt PARAMS/final_rew_mTs${k}.txt
	mv OutOpti PARAMS/OutOpti${k}
done
