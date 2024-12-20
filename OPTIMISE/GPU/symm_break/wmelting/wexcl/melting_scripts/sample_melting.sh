path=$(pwd)

if [ "$#" -ne 2 ]; then
	echo "Illegal number of parameters in melting bash script"
	exit 1
fi

nbp=$1 #configuration file
oxDNA_path=$2

box_size=20
Ncores=15
#nbp=15
#T = 5C for 5 and 8bps, 15C for 15bps
T=278.15
seq_file=gen_seqs_n${nbp}.txt

st=$(awk '{print $1}' ${seq_file})
if [[ "$st" == "" ]]; then
	echo "Sequence not specified. Usage:"
	echo "In config_file: "
	echo "SEQ ACGTCGA..."
	echo "Terminating."
	exit 1
fi

seqs=()
for word in $st; do
	seqs+=("$word")
done
#seq=$(echo $st | tr " " "\n")

Nseqs=${#seqs[@]}

echo "Sequences:"

for (( k = 0; k < ${Nseqs}; k++ )); do
	echo ${seqs[$k]}
done

Nseq_geo=${#Seq_geo[@]}


Nbatches=$((${Nseqs}/${Ncores}))
echo "batches: ${Nbatches}"

mkdir -p Melting_data/nbp${nbp}

echo "Running"

for ((m=0; m < ${Ncores}; m++)); do
	mkdir -p ${path}/core${m}
	cp ${path}/input_MD_melting ${path}/core${m}
	sed -i "s|T = .*|T = ${T}K|g" ${path}/core${m}/input_MD_melting
	sed -i "s|seed = *|seed = ${RANDOM}|g" ${path}/core${m}/input_MD_melting
	cp ${path}/oxDNA_sequence_dependent_parameters_in.txt ${path}/core${m}
done

for ((l=0; l < ${Nbatches}; l++)); do
	for ((m=0; m < ${Ncores}; m++)); do

		rep_path=${path}/core${m}
		sid=$((${l}*${Ncores}+${m}))
		cd ${rep_path}
		echo "double ${seqs[sid]}" > gen.txt
		python3 ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt
		${oxDNA_path}/build/bin/oxDNA input_MD_melting > out &
	done

	wait

	for ((m=0; m < ${Ncores}; m++)); do
		rep_path=${path}/core${m}
		sid=$((${l}*${Ncores}+${m}))
		cd ${rep_path}
		mv split_energy.dat ../Melting_data/nbp${nbp}/Energy_seq${sid}.txt
		mv trajectory.dat ../Melting_data/nbp${nbp}/Melting_traj_seq${sid}.txt
	done

done
