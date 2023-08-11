#!/bin/bash

module load anaconda/python-3.8.3/2020.07  #for archie/python3

###########################################################
#PARAMETERS - SET ACCORDINGLY!


template_path="/home/andrea/Desktop/oxDNA3_explo_analysis/template" #template folder path
sim_path="/home/andrea/Desktop/oxDNA3_explo_analysis/prova" #where the simulation folder is created
sim_type="MD" #MD or VMMC

script_gnuplot=""
compare=1 #1 to generate gnuplot scripts comparing interbnal coordinates for different sets of parameters

Nseqs=1 #number of simulated sequences
Nreps=1 #number of repetitions for each sequence
NmainPars=2 #number of main parameters (not accounting for continuity)
NparVals=1 #number of values for each parameter

#folders names, parameters values
main_folder_name="CRST_K__THETA4_T0"


declare -a parNames=("CRST_K" "CRST_THETA4_T0") #parameters to change
#the first NmainPars are the main parameters, the others the continuity parameters. XXXX NAMES MUST MATCH oxdna/src/model.h!!!
declare -A parVals #matrix with the values -> first row = values for par1, second row, value for par2, ...
parVals[0,0]="50.5f"; parVals[0,1]="51.5f"
parVals[1,0]="0.3f"; parVals[1,1]="0.6f"


declare -a colors=("black" "blue" "grey" "red") #at least as many as NparVals

#echo "${parVals[0,1]}"

for (( i=0; i<${NparVals}; i++ )); do

	#work out directory name
	folder_name=${parNames[0]}"_"${parVals[0,$i]}
	for (( j=1; j<${NmainPars}; j++ )); do
		folder_name=${folder_name}"__"${parNames[$j]}"_"${parVals[$j,$i]}
		#echo ${folder_name}
	done
	echo "entering directory ${sim_path}/${main_folder_name}/${folder_name}"
	
	data_path=${sim_path}/${main_folder_name}/${folder_name}	
	
	for (( k=0; k<${Nseqs}; k++ ));  do
	
		script_gnuplot_compare=${sim_path}/${main_folder_name}/Analysis/Internal_coord/"gnu_Seq"${k}
		output_appendix_compare="_Seq"${k}
		
		for (( j=0; j<${NmainPars}; j++ )); do 
			script_gnuplot_compare=${script_gnuplot_compare}"__"${parNames[$j]}
			output_appendix_compare=${output_appendix_compare}"__"${parNames[$j]}
		done
		script_gnuplot_compare=${script_gnuplot_compare}".pl"
		output_appendix_compare=${output_appendix_compare}".pdf"
		
		for (( l=0; l<${Nreps}; l++ )); do
			current_rep_path=${data_path}/Seq${k}/Rep${l}_${sim_type}
			
			cd ${current_rep_path}			
			
			#compute internal coordinates
			python3 ${template_path}/mapping/oxdna_to_internal_wflip.py trajectory.dat generated.top
			
		done
		
		mkdir -p ${sim_path}/${main_folder_name}/Analysis/Internal_coord/${folder_name}/Seq${k}
		cd ${sim_path}/${main_folder_name}/Analysis/Internal_coord/
		
		#average internal coordinates over repetitions for each sequence
		
		ave_seq_name="av_int_coord_Seq"${k}
		gnu_script_file=${sim_path}/${main_folder_name}/Analysis/Internal_coord/${folder_name}/Seq${k}/"gnu_Seq"${k}
		output_appendix_name="_Seq"${k}
		label=""
		
		
		for (( j=0; j<${NmainPars}; j++ )); do
			ave_seq_name=${ave_seq_name}"__"${parNames[$j]}"_"${parVals[$j,$i]}
			gnu_script_file=${gnu_script_file}"__"${parNames[$j]}"_"${parVals[$j,$i]}
			output_appendix_name=${output_appendix_name}"__"${parNames[$j]}"_"${parVals[$j,$i]}
			label=${label}"__"${parNames[$j]}"_"${parVals[$j,$i]}
			
		done
		ave_seq_name=${ave_seq_name}".txt"		
		ave_seq_file=${data_path}/Seq${k}/${ave_seq_name}
		gnu_script_file=${gnu_script_file}".pl"
		output_appendix_name=${output_appendix_name}".pdf"
		
		echo ${ave_seq_file}
		
		awk -v repn=$i 'BEGIN{
			tot_files = 0
	            }
		    {
			if(FNR == 1) {
				tot_files++;
				if(NR == FNR) {
					print $0 
				}
			}
			if(FNR != 1) {
				if(NR==FNR) bases[FNR-1]=$1
				for(i=2;i<=13;i++) a[FNR-1,i-1]+=$i
			}
		     }
		     END{
		     	for(i=1;i<FNR-1;i++) {
		     		printf "%s", bases[i]
		     		for(j=1;j<=12;j++) printf " %f", a[i,j]/tot_files
		     		printf "\n"
		     	}
		     	printf "%s", bases[FNR-1]
	     		for(j=1;j<=6;j++) printf " %f", a[i,j]/tot_files
	     		printf "\n"
	     		
	     		gnuticks="set xtics ( "
	     		for(i=1;i<FNR-1;i++) gnuticks=gnuticks "\"" bases[i] "\" " i-1 ", "
	     		gnuticks=gnuticks "\"" bases[FNR-1] "\" " FNR-2 " )"
	     		gnuxrange="set xrange [-1:"FNR-1"]"	
	     		print gnuticks > "temp.txt"
	     		print gnuxrange >> "temp.txt"
	     		
	     		if(repn == 0) {
	     			print gnuticks > "temp_compare.txt"
	     			print gnuxrange >> "temp_compare.txt"     		
		     	}
		     }' ${data_path}/Seq${k}/*/av_int_coord_mapv4.txt > ${ave_seq_file}
		     
		#cat temp_compare.txt
		#produce gnuplot script for each sequence
		  
		cp ${template_path}/gnuplot_scripts/plot_seq.pl ${gnu_script_file}
		
		echo -e "input_file=\"${ave_seq_name}\"" >> temp.txt
		echo -e "output_file=\"${output_appendix_name}\"\n" >> temp.txt
		     
		awk '{print $0}' ${gnu_script_file} >> temp.txt 
		    
		mv temp.txt ${gnu_script_file}
		
		file_final=${sim_path}/${main_folder_name}/Analysis/Internal_coord/${folder_name}/Seq${k}/${ave_seq_name}	
		cp ${ave_seq_file} ${file_final}
		
		if [ $i -eq 0 ]; then			
		     
		        echo -e "output_file=\"${output_appendix_compare}\"" >> temp_compare.txt
		     	echo -e "Nfiles=${NparVals}" >> temp_compare.txt
		     	echo -e "array files[${NparVals}]" >> temp_compare.txt
		     	echo -e "array titles[${NparVals}]" >> temp_compare.txt
		     	echo -e "array colors[${NparVals}]\n" >> temp_compare.txt			
			
			#echo -e "array files[${NparVals}]\n$(cat ${gnu_script_file_compare})" > ${gnu_script_file_compare}
		fi		
		
		echo -e "files[$(($i+1))]=\"${file_final}\"" >> temp_compare.txt
		echo -e "titles[$(($i+1))]=\"${label}\"" >> temp_compare.txt
		echo -e "colors[$(($i+1))]=\"${colors[$i]}\"\n" >> temp_compare.txt
		
		if [ $i -eq $((${NparVals}-1)) ]; then
			
			cp ${template_path}/gnuplot_scripts/plot_compare.pl ${script_gnuplot_compare}
			awk '{print $0}' ${script_gnuplot_compare} >> temp_compare.txt		    
			echo ${script_gnuplot_compare}
			mv temp_compare.txt ${script_gnuplot_compare}
		fi		  
		  
		 
	done
	
done

#produce comparative gnuplot script
