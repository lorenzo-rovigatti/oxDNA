main_path=$(pwd)
path_geo=${main_path}/GEO_SAMPLING
path_melting=${main_path}/MELTING_STEP1/Step0
oxdna_path=/home/users/yqb22156/oxDNA/oxdna

echo "Usage: "
echo "bash optimise.sh config_file_mT Nreps_geo Nsteps"
echo "This is the simple sequence version of the process"
echo "Each of the geo folders (AA AT AC ...) must contain a proper config file named opti_geometry.txt"
echo "The melting config file must be inside the path_melting"

if [ "$#" -ne 3 ]; then
	echo "Illegal number of parameters"
	exit 1
fi

config_melting=$1 #configuration file melting
Nreps_geo=$2	#number of reps for geoetry
Nsteps=$3 #number of iterations


echo "oxDNA path:"
echo $oxDNA_path
opti_path=$(awk '$1 == "OPTI_PATH" {print $2}' $path_melting/$config_melting)
if [[ "$opti_path" == "" ]]; then
	echo "Optimisation program path not specified. Usage"
	echo "In config_file: "
	echo "OPTI_PATH path"
	echo "Terminating."
	exit 1
fi

echo "Optimisation program path:"
echo $opti_path

#parameters geometry (simple sequences version)

Seq_geo=("AA" "AC" "AG" "AT" "CC" "CG")

Nseq_geo=${#Seq_geo[@]}

seq_geo_len=30

Nproc_geo=$((${Nreps_geo}*${Nseq_geo}))

echo "Number of CPUs used (geometry): ${Nproc_geo}"


#read melting parameters (setup analysis inputs)

Nreps_melting=$(awk '$1 == "REPS" {print $2}' $path_melting/$config_melting)
if [[ "$Nreps" == "" ]]; then
	echo "Number of replicas melting not specified. Setting it to 1. Usage:"
	echo "In config_file: "
	echo "REPS nreplica"
	Nreps_melting=1
fi


echo "Usage: "
echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"

st=$(awk '$1 == "SEQ" {print $2}' $path_melting/$config_melting)
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

Nseq_melting=${#seq[@]}

echo "Sequences:"

for (( k = 0; k < ${Nseq_melting}; k++ )); do
	echo ${seq[$k]}
done

#total number of cpus is number_of_sequences*number_of_reps
Nproc_melting=$((${Nreps_melting}*${Nseq_melting}))

echo "Number of CPUs used (melting): ${Nproc_melting}"

#Salt concentrations

st=$(awk '$1 == "SEQ" {print $4}' $path_melting/$config_melting)
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

if [[ ${#salt[@]} -eq ${Nseq_melting} ]]; then 
	echo "Number check ok."
else
	echo "number of salt concentrations is not equal to number of sequences."
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
exit 1
fi


#seq=$(echo $st | tr " " "\n")

echo "Salt concentrations:"

for (( k = 0; k < ${Nseq_melting}; k++ )); do
	echo ${salt[$k]}
done

#Temperatures

st=$(awk '$1 == "SEQ" {print $5}' $path_melting/$config_melting)
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

for (( k = 0; k < ${Nseq_melting}; k++ )); do
	echo ${Ts[$k]}
done

#order parameter files

st=$(awk '$1 == "SEQ" {print $8}' $path_melting/$config_melting)
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

if [[ ${#Ts[@]} -eq ${Nseq_melting} ]]; then 
	echo "Number check ok."
else
	echo "number of OP files is not equal to number of sequences."
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
exit 1
fi

echo "OP files:"

for (( k = 0; k < ${Nseq_melting}; k++ )); do
	echo ${OP_files[$k]}
done


# Weights files

st=$(awk '$1 == "SEQ" {print $9}' $path_melting/$config_melting)
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

if [[ ${#W_files[@]} -eq ${Nseq_melting} ]]; then 
	echo "Number check ok."
else
	echo "number of weight files is not equal to number of sequences."
	echo "SEQ seq box_size salt_conc simT Number_extrapolated_Ts delta_extrapolated_T op_file weight_file"
	echo "Terminating."
exit 1
fi

echo "W_files:"

for (( k = 0; k < ${Nseq_melting}; k++ )); do
	echo ${W_files[$k]}
done



##################################
# SETTING UP #####################
##################################

#melting

for ((l=0; l < ${Nseq_melting}; l++)); do
	for ((j=0; j < ${Nreps_melting};j++)); do


		rep_path=${path_melting}/Seq${l}/Rep${j}

		cp ${opti_path}/MPI/input1_melting.an ${rep_path}
		cp ${opti_path}/MPI/input2_melting.an ${rep_path}
		
		cd ${rep_path}

		sed -i "s|replace_me1|${rep_path}|g" input1_melting.an	
		sed -i "s|replace_me1|${rep_path}|g" input2_melting.an
		sed -i "s|replace_me2|${rep_path}|g" input1_melting.an	
		sed -i "s|replace_me2|${path_melting}|g" input2_melting.an
		sed -i "s|op.txt|${OP_files[$l]}|g" input1_melting.an
		sed -i "s|op.txt|${OP_files[$l]}|g" input2_melting.an
		sed -i "s|T = .*|T = ${Ts[$l]}C|g" input1_melting.an
		sed -i "s|T = .*|T = ${Ts[$l]}C|g" input2_melting.an
		sed -i "s|salt_concentration = .*|salt_concentration = ${salt[$l]}|g" input1_melting.an
		sed -i "s|salt_concentration = .*|salt_concentration = ${salt[$l]}|g" input2_melting.an
	done
done

#geometry

for ((l=0; l < ${Nseq_geo}; l++)); do
	for ((j=0; j < ${Nreps_geo};j++)); do

		seq_path=${path_geo}/${Seq_geo[$l]}
		rep_path=${path_geo}/${Seq_geo[$l]}/Seq0/Rep${j}

		cp ${opti_path}/MPI/input1.an ${rep_path}
		cp ${opti_path}/MPI/input2.an ${rep_path}
		
		cd ${rep_path}

		sed -i "s|replace_me1|${rep_path}|g" input1.an	
		sed -i "s|replace_me1|${rep_path}|g" input2.an
		sed -i "s|replace_me2|${rep_path}|g" input1.an	
		sed -i "s|replace_me2|${seq_path}|g" input2.an
		python3 ${opti_path}/MPI/GenOP.py ${seq_geo_len}	#generate op file (to check if sampled confs are melted)
	done
done

#####################################
# RUNNING CYCLE #####################
#####################################



mkdir -p ${main_path}/OPTI_FILES/PARAMS
mkdir -p ${main_path}/OPTI_FILES/OPTI_OUTPUTS


##############################
# STEP 0        ##############
##############################

i=0

#We start from the melting
cd ${path_melting}

#cp oxDNA_sequence_dependent_parameters_in.txt ${main_path}/OPTI_FILES/PARAMS/oxDNA_sequence_dependent_parameters_in.txt

#mpirun -np ${Nproc_melting} python3 ${opti_path}/MPI/optimise_melting_mpi.py $path_melting/$config_melting 0 > OutOpti

#wait

#cp oxDNA_sequence_dependent_parameters_tmp.txt ${main_path}/OPTI_FILES/PARAMS/oxDNA_sequence_dependent_parameters_fin_melting${i}.txt
#mv OutOpti ${main_path}/OPTI_FILES/OPTI_OUTPUTS/OutOpti_melting${i}
#mv parameters_v_iter.txt ${main_path}/OPTI_FILES/OPTI_OUTPUTS/parameters_v_iter_melting${i}.txt

#geometry
for ((l=0; l < ${Nseq_geo}; l++)); do
	
	#preparing par files
	cp ${path_melting}/oxDNA_sequence_dependent_parameters_tmp.txt ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_in.txt
	
	#fxed parameters
	
	cd ${path_geo}/${Seq_geo[$l]}/
	python3 ${opti_path}/Utils/prepare_fixed_pars_file_from_config.py ${path_geo}/${Seq_geo[$l]}/opti_geometry.txt oxDNA_sequence_dependent_parameters_in.txt
	
	for ((m=0; m < ${Nreps_geo}; m++)); do
		cp ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_in.txt ${path_geo}/${Seq_geo[$l]}/Seq0/Rep${m}/oxDNA_sequence_dependent_parameters_in.txt
		cd ${path_geo}/${Seq_geo[$l]}/Seq0/Rep${m}/
		${oxdna_path}/build/bin/oxDNA input_MD > out &
	done

done

wait


for ((l=0; l < ${Nseq_geo}; l++)); do
	cd ${path_geo}/${Seq_geo[$l]}/
	mpirun -np ${Nreps_geo} python3 ${opti_path}/MPI/optimise_geometry_mpi.py opti_geometry.txt > OutOpti &
done

wait

mkdir -p ${main_path}/OPTI_FILES/PLOTS/

for ((l=0; l < ${Nseq_geo}; l++)); do
	cd ${path_geo}/${Seq_geo[$l]}/
	cp Seq0_Cov_diag_propeller.pdf ${main_path}/OPTI_FILES/PLOTS/Cov_diag_propeller_${Seq_geo[$l]}_$i.pdf
	cp Seq0_Cov_diag_rise.pdf ${main_path}/OPTI_FILES/PLOTS/Cov_diag_rise_${Seq_geo[$l]}_$i.pdf
	cp Seq0_Cov_diag_roll.pdf ${main_path}/OPTI_FILES/PLOTS/Cov_diag_roll_${Seq_geo[$l]}_$i.pdf
	cp Seq0_Cov_diag_twist.pdf ${main_path}/OPTI_FILES/PLOTS/Cov_diag_twist_${Seq_geo[$l]}_$i.pdf
	cp Seq0_propeller.pdf ${main_path}/OPTI_FILES/PLOTS/propeller_${Seq_geo[$l]}_$i.pdf
	cp Seq0_rise.pdf ${main_path}/OPTI_FILES/PLOTS/rise_${Seq_geo[$l]}_$i.pdf
	cp Seq0_roll.pdf ${main_path}/OPTI_FILES/PLOTS/roll_${Seq_geo[$l]}_$i.pdf
	cp Seq0_twist.pdf ${main_path}/OPTI_FILES/PLOTS/twist_${Seq_geo[$l]}_$i.pdf
done

#assemble par file for melting

for ((l=0; l < ${Nseq_geo}; l++)); do

	cd ${path_geo}/${Seq_geo[$l]}/

	if [ $l -eq 0 ]; then
		#head -n 1 ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_tmp.txt > ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt #stck_eps
		cat ${main_path}/Fixed_SD_pars.txt > ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt
		python3 ${opti_path}/Utils/get_opti_pars_from_config.py ${path_melting}/$config_melting ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_tmp.txt
		cat ${path_geo}/${Seq_geo[$l]}/optimised_parameters.txt >> ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt
	fi
	
	python3 ${opti_path}/Utils/get_opti_pars_from_config.py ${path_geo}/${Seq_geo[$l]}/opti_geometry.txt ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_tmp.txt
	cat ${path_geo}/${Seq_geo[$l]}/optimised_parameters.txt >> ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt
	
	mv ${path_geo}/${Seq_geo[$l]}/OutOpti ${main_path}/OPTI_FILES/OPTI_OUTPUTS/OutOpti_geo${i}_${Seq_geo[$l]}
	#mv ${path_geo}/${Seq_geo[$l]}/parameters_v_iter.txt ${main_path}/OPTI_FILES/OPTI_OUTPUTS/parameters_v_iter_geo${i}_${Seq_geo[$l]}.txt
done

cp ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt ${main_path}/OPTI_FILES/PARAMS/oxDNA_sequence_dependent_parameters_fin_geo${i}.txt

#cp ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt ${path_melting}/oxDNA_sequence_dependent_parameters_in.txt
#python3 ${opti_path}/Utils/prepare_fixed_pars_file_from_config.py $config_melting oxDNA_sequence_dependent_parameters_in.txt


##############################
# CYCLE         ##############
##############################


for ((i=1; i<${Nsteps}; i++)); do

	cd ${path_melting}

	#melting
	cp ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt ${path_melting}/oxDNA_sequence_dependent_parameters_in.txt
	python3 ${opti_path}/Utils/prepare_fixed_pars_file_from_config.py $path_melting/$config_melting oxDNA_sequence_dependent_parameters_in.txt

	mpirun -np ${Nproc_melting} python3 ${opti_path}/MPI/optimise_melting_mpi.py $path_melting/$config_melting 0 > OutOpti

	wait

	cp oxDNA_sequence_dependent_parameters_tmp.txt ${main_path}/OPTI_FILES/PARAMS/oxDNA_sequence_dependent_parameters_fin_melting${i}.txt
	mv OutOpti ${main_path}/OPTI_FILES/OPTI_OUTPUTS/OutOpti_melting${i}
	mv parameters_v_iter.txt ${main_path}/OPTI_FILES/OPTI_OUTPUTS/parameters_v_iter_melting${i}.txt

	#geometry
	for ((l=0; l < ${Nseq_geo}; l++)); do
		
		#preparing par files
		cp ${path_melting}/oxDNA_sequence_dependent_parameters_tmp.txt ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_in.txt
		
		#fxed parameters
		
		cd ${path_geo}/${Seq_geo[$l]}/
		python3 ${opti_path}/Utils/prepare_fixed_pars_file_from_config.py ${path_geo}/${Seq_geo[$l]}/opti_geometry.txt oxDNA_sequence_dependent_parameters_in.txt
		
		for ((m=0; m < ${Nreps_geo}; m++)); do
			cp ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_in.txt ${path_geo}/${Seq_geo[$l]}/Seq0/Rep${m}/oxDNA_sequence_dependent_parameters_in.txt
			cd ${path_geo}/${Seq_geo[$l]}/Seq0/Rep${m}/
			${oxdna_path}/build/bin/oxDNA input_MD > out &
		done
	done
	
	wait	
		
		
	for ((l=0; l < ${Nseq_geo}; l++)); do
		cd ${path_geo}/${Seq_geo[$l]}/
		mpirun -np ${Nreps_geo} python3 ${opti_path}/MPI/optimise_geometry_mpi.py opti_geometry.txt > OutOpti &
	done
	
	for ((l=0; l < ${Nseq_geo}; l++)); do
		cd ${path_geo}/${Seq_geo[$l]}/
		cp Seq0_Cov_diag_propeller.pdf ${main_path}/OPTI_FILES/PLOTS/Cov_diag_propeller_${Seq_geo[$l]}_$i.pdf
		cp Seq0_Cov_diag_rise.pdf ${main_path}/OPTI_FILES/PLOTS/Cov_diag_rise_${Seq_geo[$l]}_$i.pdf
		cp Seq0_Cov_diag_roll.pdf ${main_path}/OPTI_FILES/PLOTS/Cov_diag_roll_${Seq_geo[$l]}_$i.pdf
		cp Seq0_Cov_diag_twist.pdf ${main_path}/OPTI_FILES/PLOTS/Cov_diag_twist_${Seq_geo[$l]}_$i.pdf
		cp Seq0_propeller.pdf ${main_path}/OPTI_FILES/PLOTS/propeller_${Seq_geo[$l]}_$i.pdf
		cp Seq0_rise.pdf ${main_path}/OPTI_FILES/PLOTS/rise_${Seq_geo[$l]}_$i.pdf
		cp Seq0_roll.pdf ${main_path}/OPTI_FILES/PLOTS/roll_${Seq_geo[$l]}_$i.pdf
		cp Seq0_twist.pdf ${main_path}/OPTI_FILES/PLOTS/twist_${Seq_geo[$l]}_$i.pdf
	done
	
	
	

	#assemble par file for melting

	for ((l=0; l < ${Nseq_geo}; l++)); do
	
		cd ${path_geo}/${Seq_geo[$l]}/

		if [ $l -eq 0 ]; then
			#head -n 1 ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_tmp.txt > ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt #stck_eps
			cat ${main_path}/Fixed_SD_pars.txt > ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt
			python3 ${opti_path}/Utils/get_opti_pars_from_config.py ${path_melting}/$config_melting ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_tmp.txt
			cat ${path_geo}/${Seq_geo[$l]}/optimised_parameters.txt >> ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt
		fi
		
		python3 ${opti_path}/Utils/get_opti_pars_from_config.py ${path_geo}/${Seq_geo[$l]}/opti_geometry.txt ${path_geo}/${Seq_geo[$l]}/oxDNA_sequence_dependent_parameters_tmp.txt
		cat ${path_geo}/${Seq_geo[$l]}/optimised_parameters.txt >> ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt
		
		mv ${path_geo}/${Seq_geo[$l]}/OutOpti ${main_path}/OPTI_FILES/OPTI_OUTPUTS/OutOpti_geo${i}_${Seq_geo[$l]}
		#mv ${path_geo}/${Seq_geo[$l]}/parameters_v_iter.txt ${main_path}/OPTI_FILES/OPTI_OUTPUTS/parameters_v_iter_geo${i}_${Seq_geo[$l]}.txt
	done

	cp ${path_geo}/oxDNA_sequence_dependent_parameters_fin.txt ${main_path}/OPTI_FILES/PARAMS/oxDNA_sequence_dependent_parameters_fin_geo${i}.txt

done

echo "ALL DONE"
