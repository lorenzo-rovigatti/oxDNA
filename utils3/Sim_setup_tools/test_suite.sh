#!/bin/bash

#module load anaconda/python-3.8.3/2020.07  #for archie!

###########################################################
#PARAMETERS - SET ACCORDINGLY!


template_path="/home/andrea/Desktop/oxDNA3_explo_analysis/template" #template folder path
sim_path="/home/andrea/Desktop/oxDNA3_explo_analysis/prova" #where the simulation folder is created

Nseqs=5 #number of simulated sequences
Nreps=1 #number of repetitions for each sequence
NmainPars=2 #number of main parameters to change (not accounting for continuity)
Npars=2 #total number of parameters to change (accounting for continuity)
NparVals=4 #number of sets of values

declare -a input_script=("input_MD_w_stacking_prop" "input_VMMC" "input_VMMC" "input_VMMC" "input_VMMC") #see template/oxdna_input_scripts. seed updated and set to random
#Each sequence has its personalised input script (e.g. VMMC or MD)

declare -a umbrella=(0 1 1 1 1) #1 if umbrella is active. If umbrella==1, copy order parameter and weights and update input accordingly
declare -a op_file=("" "op_Nbp17.txt" "op_Nbp17_plus_mismatched.txt" "op_Nbp17_plus_dangling.txt" "op_hairpin20_6.txt")
declare -a w_file=("" "w_Nbp17.txt" "w_Nbp17.txt" "w_Nbp17.txt" "w_hairpin20_6.txt")

declare -a input_conf_method=("32bp_all_steps_twice" "16bp_all_steps" "16bp_plus_mismatched" "16bp_plus_dangling" "hairpin")  #Each sequence is generated in the specified way. Array of length Nseqs. See template/in_confs/${input_conf_method[$i]}
declare -a NgenPar=(1 1 1 1 1) #number of parameters of the specified method
declare -A gen_parNames #matrix with the parameters of the input conf generation. See template/in_confs/${input_conf_method[$i]}/gen.sh. Note that template_path is always updated
gen_parNames[0,0]="BOX_SIZE"
gen_parNames[1,0]="BOX_SIZE"
gen_parNames[2,0]="BOX_SIZE"
gen_parNames[3,0]="BOX_SIZE"
gen_parNames[4,0]="TYPE"
declare -A gen_parVals #values of the parameters
gen_parVals[0,0]="50."
gen_parVals[1,0]="50."
gen_parVals[2,0]="50."
gen_parVals[3,0]="50."
gen_parVals[4,0]="pur_pyr"


#strating oxdna version
oxdna_path="/home/andrea/Desktop/oxDNA3_explo_analysis/CROSS_STK/oxdna_asymm_crstk"
oxdna_version="oxdna_asymm_crstk"
model_default="model_ox2.h" #base model to change

#folders names, parameters values
main_folder_name="CRST_K_33__CRST_K_55/TEST_SUITE"

declare -a parNames=("CRST_K_33" "CRST_K_55") #parameters to change
#the first NmainPars are the main parameters, the others the continuity parameters. XXXX NAMES MUST MATCH oxdna/src/model.h!!!
declare -A parVals #matrix with the values -> first row = values for par1, second row, value for par2, ...
parVals[0,0]="0.f"; parVals[0,1]="47.5f"; parVals[0,2]="95.f"; parVals[0,3]="0.f"
parVals[1,0]="0.f"; parVals[1,1]="47.5f"; parVals[1,2]="0.f"; parVals[1,3]="95.f"
#parVals[2,0]="1.4f"; parVals[2,1]="1.8f"

#echo "${parVals[0,1]}"



#############################################################
#SET UP SIMULATION

mkdir -p ${sim_path}/${main_folder_name}

#setup archie script (part1)
cp ${template_path}/sub_archie.sh ${sim_path}/${main_folder_name}
echo "simpath=${sim_path}/${main_folder_name}" >> ${sim_path}/${main_folder_name}/sub_archie.sh
echo "" >> ${sim_path}/${main_folder_name}/sub_archie.sh
ntsk=$((${NparVals}*${Nseqs}*${Nreps}))
sed -i "s/ntasks=4/ntasks=${ntsk}/g" ${sim_path}/${main_folder_name}/sub_archie.sh
jname="ab_"${main_folder_name}
sed -i "s/jname/${jname}/g" ${sim_path}/${main_folder_name}/sub_archie.sh



for (( i=0; i<${NparVals}; i++ )); do

	#create folder
	echo "creating directory i=${i}"
	folder_name=${parNames[0]}"_"${parVals[0,$i]}
	#echo ${folder_name}
	for (( j=1; j<${NmainPars}; j++ )); do
		echo
		folder_name=${folder_name}"__"${parNames[$j]}"_"${parVals[$j,$i]}
		echo ${folder_name}
	done
	mkdir -p ${sim_path}/${main_folder_name}/${folder_name}
	
	echo "copying oxdna i=${i}\n"
	#copy oxdna (each different set of parameters requires its own copy of oxdna)
	rm -r ${sim_path}/${main_folder_name}/${folder_name}/${oxdna_version}
	cp -r ${oxdna_path} ${sim_path}/${main_folder_name}/${folder_name}
	
	#modify model
	echo "updating model i=${i}"
	cp ${template_path}/models/${oxdna_version}/${model_default} ${sim_path}/${main_folder_name}/${folder_name}
	for (( j=0; j<${Npars}; j++ )); do
		sed -i "s/.*${parNames[$j]}.*/#define ${parNames[$j]} ${parVals[$j,$i]}/" ${sim_path}/${main_folder_name}/${folder_name}/${model_default}
	done
	
	#update current model in template/models (it is needed to generate the initial configuration)
	cp ${sim_path}/${main_folder_name}/${folder_name}/${model_default} ${template_path}/models/${oxdna_version}/model.h
	#save copy in models/template for future use
	echo "saving model i=${i} in template/models"
	model_name="model_"${parNames[0]}"_"${parVals[0,$i]}
	for (( j=1; j<${NmainPars}; j++ ));  do
		model_name=${model_name}"__"${parNames[$j]}"_"${parVals[$j,$i]}
	done
	model_name=${model_name}".h"
	
	#save current parameters in text file
	for (( j=0; j<${Npars}; j++ )); do
		line=${parNames[$j]}" "${parVals[$j,$i]}
		echo $line >> ${sim_path}/${main_folder_name}/${folder_name}/parameters.txt
	done
	
	cp ${template_path}/models/${oxdna_version}/model.h ${template_path}/models/${oxdna_version}/${model_name}
	
	#update model.h in oxdna directory
	mv ${sim_path}/${main_folder_name}/${folder_name}/${model_default} ${sim_path}/${main_folder_name}/${folder_name}/${oxdna_version}/src/model.h
	
	#build oxdna
	echo "compiling oxdna i=${i}"
	cd ${sim_path}/${main_folder_name}/${folder_name}/${oxdna_version}/
	mkdir build
	cd build/
	cmake ..
	make -j8
	
	#initial configurations and input scripts
	echo "creating input_confgurations and oxdna input scripts oxdna i=${i}"

	for (( k=0; k<${Nseqs}; k++ ));  do
		for (( l=0; l<${Nreps}; l++ )); do
			current_rep_path=${sim_path}"/"${main_folder_name}"/"${folder_name}"/Seq"${k}"/Rep"${l}
			mkdir -p ${current_rep_path} #make rep directory

			cd ${current_rep_path}			
			
			#generate initial configurations
			
			cp ${template_path}/in_confs/${input_conf_method[$k]}/gen.sh .
			sed -i "s|template_path=.*|template_path=${template_path}|g" gen.sh
			
			for(( z=0; z<${NgenPar[$k]}; z++ ));  do
				sed -i "s|${gen_parNames[$k,$z]}=.*|${gen_parNames[$k,$z]}=${gen_parVals[$k,$z]}|g" gen.sh
			done
			bash gen.sh
			
			#copy input script (and change seed)
			cp ${template_path}/oxdna_input_scripts/${input_script[$k]} .
			sed -i "s/seed = .*/seed = ${RANDOM}/g" ./${input_script[$k]}
			sed -i "s|seq_dep_file = .*|seq_dep_file = ${template_path}/${oxdna_version}/oxDNA2_sequence_dependent_parameters.txt|g" ./${input_script[$k]}
			#if umbrella, copy order p and weights files and update input with correct names
			if [ ${umbrella[$k]} -eq 1 ]; then
				cp ${template_path}/oxdna_input_scripts/order_parameters/${op_file[$k]} .
				sed -i "s/op_file = .*/op_file = ${op_file[$k]}/g" ./${input_script[$k]}
				cp ${template_path}/oxdna_input_scripts/order_parameters/${w_file[$k]} .
				sed -i "s/weights_file = .*/weights_file = ${w_file[$k]}/g" ./${input_script[$k]}
			fi
			
			#append to archie script
			echo "cd \${simpath}/${folder_name}/Seq${k}/Rep${l} ; \${simpath}/${folder_name}/${oxdna_version}/build/bin/oxDNA ${input_script[$k]} &" >> ${sim_path}/${main_folder_name}/sub_archie.sh
			
			
		done
	done
	
	echo "finalising archie script i=${i}"
	
done

echo "" >> ${sim_path}/${main_folder_name}/sub_archie.sh
echo 'wait' >> ${sim_path}/${main_folder_name}/sub_archie.sh
