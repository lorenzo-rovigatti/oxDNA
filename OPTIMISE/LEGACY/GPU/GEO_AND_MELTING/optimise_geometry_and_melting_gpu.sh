#when running this optimisation, there must be two directories inside the main folder:
#Melting and Geometry.
#Melting must contain the sampled trajectiories inside the 3 folders n5, n8 and n15, in the usual format (i.e. in separated folders /Seq${$
#In Melting there should be also the initial oxDNA_sequence_dependent_parameters_in.txt
#both Melting and Geometry must contain their optim.txt configuration files, which are the arguments of this script


main_path=$(pwd)
echo "Usage: "
echo "bash optimise.sh config_file conf_file_melting"

if [ "$#" -ne 2 ]; then
	echo "Illegal number of parameters"
	exit 1
fi

config=$1 #configuration file geometry
config_melt=$2 #configuration file melting

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
opti_path_geo=$(awk '$1 == "OPTI_PATH_GEO" {print $2}' $config)
if [[ "$opti_path_geo" == "" ]]; then
	echo "Geometry optimisation program path not specified. Usage"
	echo "In config_file: "
	echo "OPTI_PATH_GEO path"
	echo "Terminating."
	exit 1
fi

opti_path_melt=$(awk '$1 == "OPTI_PATH_MELT" {print $2}' $config)
if [[ "$opti_path_melt" == "" ]]; then
    echo "Melting optimisation program path not specified. Usage"
    echo "In config_file: "
    echo "OPTI_PATH_MELT path"
    echo "Terminating."
    exit 1
fi

echo "Optimisation program paths:"
echo $opti_path_geo
echo $opti_path_melt

box_size=$(awk '$1 == "BOX_SIZE" {print $2}' $config)
if [[ "$box_size" == "" ]]; then
	echo "Box size not specified. Setting it to 20. Usage:"
	echo "In config_file: "
	echo "BOX_SIZE box_size"
	box_size=20
fi

echo "Sim box size:"
echo $box_size

concentration=$(( 2687/$box_size/$box_size/$box_size ))

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
Nsteps_geo=$(awk '$1 == "NSTEPS_GEO" {print $2}' $config)
if [[ "$Nsteps_geo" == "" ]]; then
    echo "Number of optimisation steps. Setting it to 1. Usage:"
    echo "In config_file: "
    echo "NSTEPS_GEO nsteps"
    Nsteps_geo=1
fi
echo "# of optimisation steps:"
echo $Nsteps_geo

timesteps=$(awk '$1 == "TIMESTEPS" {print $2}' $config)
if [[ "$timesteps" == "" ]]; then
	echo "Number of timesteps not specified. Setting it to 1e7. Usage:"
	echo "In config_file: "
	echo "TIMESTEPS timesteps"
	timesteps=1e7
fi
echo "# of sim timesteps:"
echo $timesteps

temperature=$(awk '$1 == "TEMPERATURE" {print $2}' $config)
if [[ "$temperature" == "" ]]; then
	echo "teperature not specified. Setting it to 300K. Usage:"
	echo "In config_file: "
	echo "TIMESTEPS timesteps"
	temperature=300K
fi
echo "# of sim timesteps:"
echo $timesteps

#other parameters

#generate gen.txt file for oxDNA

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

echo "Sequences:"
echo ${seq[0]}
echo ${seq[1]}

Nseq=${#seq[@]}


#total number of cpus is number_of_sequences*number_of_reps
Nproc=$((${Nreps}*${Nseq}))

echo "Number of CPUs used: ${Nproc}"

#initialise first step

#create directories for storing opti files

pars_path=${main_path}/OPTI_PARS/
logs_path=${main_path}/OPTI_LOGS/
figs_path=${main_path}/OPTI_FIGS/

mkdir ${pars_path}
mkdir ${logs_path}
mkdir ${figs_path}

geo_path=${main_path}/Geometry
melt_path=${main_path}/Melting

#run first melting optimisation and store parsed trajectories/enrgies

#store tensors (trajectories and initial energies) to pytorch file; avoid overhead and need for 3 parameters file (1 extra for en0)
sed -i "s|print_energy_to_file .*|print_energy_to_file True|g" ${main_path}/${config_melt}
sed -i "s|print_coords_to_file .*|print_coords_to_file True|g" ${config_melt}
sed -i "s|read_energy_from_file .*|read_energy_from_file False|g" ${config_melt}
sed -i "s|read_coords_from_file .*|read_coords_from_file False|g" ${config_melt}

#sed -i "s|print_energy_to_file .*|print_energy_to_file False|g" ${config_melt}
#sed -i "s|print_coords_to_file .*|print_coords_to_file False|g" ${config_melt}
#sed -i "s|read_energy_from_file .*|read_energy_from_file True|g" ${config_melt}
#sed -i "s|read_coords_from_file .*|read_coords_from_file True|g" ${config_melt}

cd ${melt_path}

python ${opti_path_melt}/optimise_melting_gpu.py ${main_path}/$config_melt > OutOpti.log

cp oxDNA_sequence_dependent_parameters_in.txt ${pars_path}/oxDNA_sequence_dependent_parameters_in0.txt
cp oxDNA_sequence_dependent_parameters_fin.txt ${pars_path}/oxDNA_sequence_dependent_parameters_M0.txt
cp OutOpti.log ${logs_path}/OutOpti_M0.log

cd ..

#switch to read mode in the melting for all other steps.
sed -i "s|print_energy_to_file .*|print_energy_to_file False|g" ${config_melt}
sed -i "s|print_coords_to_file .*|print_coords_to_file False|g" ${config_melt}
sed -i "s|read_energy_from_file .*|read_energy_from_file True|g" ${config_melt}
sed -i "s|read_coords_from_file .*|read_coords_from_file True|g" ${config_melt}

#switch to geometry and run cycle

for((s=0; s <${Nsteps}; s++)); do

    #Geometry step 0 of step s

    stepid=$(($s+1))

    cp ${melt_path}/oxDNA_sequence_dependent_parameters_fin.txt ${geo_path}/oxDNA_sequence_dependent_parameters_in.txt
    #maybe separating step0 and the other is not necessary. Don't have time to adapt things
    for ((l=0; l < ${Nseq}; l++)); do
    	for ((j=0; j < ${Nreps};j++)); do
    		mkdir -p ${geo_path}/Step0/Seq${l}/Rep${j}/
    		cp ${opti_path_geo}/input_MD ${geo_path}/Step0/Seq${l}/Rep${j}/
            cp ${geo_path}/oxDNA_sequence_dependent_parameters_in.txt ${geo_path}/Step0/
    		cp ${geo_path}/Step0/oxDNA_sequence_dependent_parameters_in.txt ${geo_path}/Step0/Seq${l}/Rep${j}/

    		rep_path=${geo_path}/Step0/Seq${l}/Rep${j}
    		step_path=${geo_path}/Step0/
    		cd ${rep_path}

    		echo "double "${seq[$l]} > gen.txt	#used for generating in conf and topology

    		python ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt	#gen topology and in conf
    		sed -i "s|seed = 1|seed = ${RANDOM}|g" input_MD		#randomise seed
    		sed -i "s|steps = 1e7|steps = ${timesteps}|g" input_MD
    		sed -i "s|T = 300K|T = ${temperature}|g" input_MD
    		sed -i "s|conf_file = .*|conf_file = generated.dat|g" input_MD
    		${oxDNA_path}/build/bin/oxDNA input_MD > out_main &
    	done
    done

    wait


    #optimise
    cd ${geo_path}/Step0
    python ${opti_path_geo}/optimise_geometry_gpu.py ${main_path}/$config > OutOpti.log

    cp OutOpti.log ${logs_path}/OutOpti_G${stepid}_S0.log
    cp oxDNA_sequence_dependent_parameters_fin.txt ${pars_path}/oxDNA_sequence_dependent_parameters_G${stepid}_S0.txt
    mkdir -p ${figs_path}/G${stepid}/Step0/
    cp *.pdf ${figs_path}/G${stepid}/Step0/


    #same as above, but for all the other geometry steps
    #copy all necessary files form previous step and
    #turn output (param file) of optimise into input param file for next step

    for ((i=1; i < ${Nsteps_geo}; i++)); do

    	mkdir -p ${geo_path}/Step${i}/

    	#copy optimisation code
    	cp ${geo_path}/oxDNA_sequence_dependent_parameters.txt ${geo_path}/Step${i}/
    	cp ${geo_path}/Step$(($i-1))/oxDNA_sequence_dependent_parameters_fin.txt ${geo_path}/Step${i}/oxDNA_sequence_dependent_parameters_in.txt	#copy output of optimse at previous step and make it the starting parameter files of this step
    	for ((l=0; l < ${Nseq}; l++)); do
    		for ((j=0; j< ${Nreps};j++)); do
    			mkdir -p ${geo_path}/Step${i}/Seq${l}/Rep${j}/
    			cp ${opti_path_geo}/input_MD ${main_path}/Step${i}/Seq${l}/Rep${j}/
    			cp ${geo_path}/Step$(($i-1))/oxDNA_sequence_dependent_parameters_fin.txt ${main_path}/Step${i}/Seq${l}/Rep${j}/
    			#cp ${main_path}/Step$(($i-1))/Seq${l}/Rep${j}/op.txt ${main_path}/Step${i}/Seq${l}/Rep${j}/
    			cd ${geo_path}/Step${i}/Seq${l}/Rep${j}/
    			mv oxDNA_sequence_dependent_parameters_fin.txt oxDNA_sequence_dependent_parameters_in.txt

    			cd ${geo_path}/Step${i}/
    			rep_path=$(pwd)/Seq${l}/Rep${j}
    			step_path=$(pwd)
    			cd ${rep_path}

    			echo "double "${seq[$l]} > gen.txt

    			python ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt
    			sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD
    			sed -i "s|steps = 1e7|steps = ${timesteps}|g" input_MD
    			sed -i "s|T = 300K|T = ${temperature}|g" input_MD
    			sed -i "s|conf_file = .*|conf_file = generated.dat|g" input_MD
    			${oxDNA_path}/build/bin/oxDNA input_MD > out_main &
    		done
    	done

        wait


    	cd ${geo_path}/Step${i}
    	python ${opti_path_geo}/optimise_geometry_gpu.py ${main_path}/$config > OutOpti.log

        cp OutOpti.log ${logs_path}/OutOpti_G${stepid}_S${i}.log
        cp oxDNA_sequence_dependent_parameters_fin.txt ${pars_path}/oxDNA_sequence_dependent_parameters_G${stepid}_S${i}.txt
        mkdir -p ${figs_path}/G${stepid}/Step${i}/
        cp *.pdf ${figs_path}/G${stepid}/Step${i}/

    done

    #switch to melting.
    #trajectories and initial energies are read from pytorch tensor files this time.

    cp ${geo_path}/Step$((${Nsteps_geo}-1))/oxDNA_sequence_dependent_parameters_fin.txt ${melt_path}/oxDNA_sequence_dependent_parameters_in.txt

    cd ${melt_path}

    python ${opti_path_melt}/optimise_melting_gpu.py ${main_path}/$config_melt > OutOpti.log

    cp oxDNA_sequence_dependent_parameters_fin.txt ${pars_path}/oxDNA_sequence_dependent_parameters_M${stepid}.txt
    cp OutOpti.log ${logs_path}/OutOpti_M${stepid}.log

done
