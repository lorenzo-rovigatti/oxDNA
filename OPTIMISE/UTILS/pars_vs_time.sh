#Read and save how parameters change during an optimisation cycle. Postprocessing.
#Only sampled steps are considered (for chenges in the parameters while reweighting, see output of optimise.py)
#Usage:
#in sim folder (where the Steps forders are located)
# bash pars_vs_time.sh input steps delta
#steps = how many sampled steps
#delta = save only every selta steps

echo "Usage: "
echo "bash pars_vs_time.sh input steps delta"

if [ "$#" -ne 3 ]; then
	echo "Illegal number of parameters"
	exit 1
fi

name_file=$1 #configuration file
steps=$2
delta=$3

path=$(pwd)

out=${path}/pars_vs_time.txt
in0=oxDNA_sequence_dependent_parameters_in.txt
in=oxDNA_sequence_dependent_parameters_fin.txt

names=($(ps | awk '$3=="OPTIMISE" { print $1 }' ${path}/${names_file}))

for i in ${names[@]}; do

echo $i

done

for ((i=0; i < ${steps}; i+=${delta})); do
	declare -a values=()
	
	if [ $i -eq 0 ]; then
		for j in ${names[@]}; do
			values+=($(ps | awk -v na="${j}" '$1==na { print $3 }' ${path}/Step${i}/${in0}))
	done
	
	outstr=""
	for j in ${values[@]}; do
		outstr+=" "$j
	done
	echo ${outstr} >> ${out}
	#echo "\n" >> ${out}
	
	fi
	
	declare -a values=()
	
	for j in ${names[@]}; do
		values+=($(ps | awk -v na="${j}" '$1==na { print $3 }' ${path}/Step${i}/${in}))
	done
	outstr=""
	for j in ${values[@]}; do
		outstr+=" "$j
	done
	echo ${outstr} >> ${out}
	#echo "\n" >> ${out}
	
done
