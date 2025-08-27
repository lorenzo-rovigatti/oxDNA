path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

#read sequences
Ncores=10

q0=0.815 #inital charge
strength=0.0817494 #dh prefactor/q0**2

dq=0.02

forward_steps=20
backward_steps=30

declare -a seqs Cs Ts0
i=0

declare -a Ts

while read -r a b c; do
    seqs[i]=$a
    Cs[i]=$b
    Ts0[i]=$c
    ((i++))
done < <(tail -n +2 seqs.dat)  # Skip header with tail -n +2

Nseqs=${#seqs[@]} #number of sequences

Nbs=() #number of bases
for str in "${seqs[@]}"; do
    Nbs+=( $((${#str}*2)) ) #2 strands in both bn and un -> *2
done

Ts=${Ts0[@]}
q=$q0

#forward integration

for((j=0; j<${forward_steps};j++)); do
	counts=0
	for((i=0; i<${Nseqs};i++)); do

		counts=$(($counts+2))

		sim_path=${path}/$i

                if [ $j -eq 0 ]; then
                	mkdir -p $sim_path
			cp input_MD_bn $sim_path
			cp input_MD_un $sim_path
			
			cd $sim_path
			echo "double "${seqs[$i]} > gen.txt
			python ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt	#gen topology and in conf for bound state
			
			cp oxDNA3_sequence_dependent_parameters_bn.txt ${sim_path}
			cp oxDNA3_sequence_dependent_parameters_un.txt ${sim_path} #here we switch off hydrogen bonding = unbond state
			
			sed -i "s|salt_concentration = .*|salt_concentration = ${Cs[$i]}|g" input_MD_bn	
			sed -i "s|salt_concentration = .*|salt_concentration = ${Cs[$i]}|g" input_MD_un
			
			echo "${q} ${Ts[$i]}" >> integrated_mT_vs_q.dat
			
		fi
		
		cd $sim_path
		
		#set temperature and charge
		
		# Read last line of integrated_mT_vs_q.dat (this stores previus step final q and Tm)
		last_line=($(tail -n 1 integrated_mT_vs_q.dat})

		q=${last_line[0]}
		Ts[$i]=${last_line[1]}
		
		sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD_bn
		sed -i "s|T = .*|T = ${Ts[$i]}|g" input_MD_bn
		sed -i "s|dh_strength = .*|dh_strength = $(($q*$q*$strength))|g" input_MD_bn
		
		sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD_un
		sed -i "s|T = .*|T = ${Ts[$i]}|g" input_MD_un
		sed -i "s|dh_strength = .*|dh_strength = $(($q*$q*$strength))|g" input_MD_un

                ${oxdna_path}/build/bin/oxDNA input_MD_bn > out &
		${oxdna_path}/build/bin/oxDNA input_MD_un > out &

		if [ $counts -eq ${Ncores} ]; then
			counts=0
			wait
		fi

	done
	
	#integrate step
	counts=0
	for((i=0; i<${Nseqs};i++)); do

		counts=$(($counts+1))

		sim_path=${path}/$(($counts-1))
		
		cd $sim_path
		
		if [ $counts -eq ${Ncores} ]; then
			counts=0
			wait
		fi
		
		python3 ../integrate_step.py ${Ts[$i]} $q $dq ${Nbs[$i]} &

	done

done

wait

#backwards integration

for((j=0; j<${backwards_steps};j++)); do
	counts=0
	for((i=0; i<${Nseqs};i++)); do

		counts=$(($counts+1))

		sim_path=${path}/$i
		
		cd $sim_path
		
		#set temperature and charge
		
		# Read first line of integrated_mT_vs_q.dat (this stores previus step final q and Tm, when going backwards)
		first_line=($(head -n 1 integrated_mT_vs_q.dat))

		q=${first_line[0]}
		Ts[$i]=${first_line[1]}
		
		sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD_bn
		sed -i "s|T = .*|T = ${Ts[$i]}|g" input_MD_bn
		sed -i "s|dh_strength = .*|dh_strength = $(($q*$q*$strength))|g" input_MD_bn
		
		sed -i "s|seed = .*|seed = ${RANDOM}|g" input_MD_un
		sed -i "s|T = .*|T = ${Ts[$i]}|g" input_MD_un
		sed -i "s|dh_strength = .*|dh_strength = $(($q*$q*$strength))|g" input_MD_un

                ${oxdna_path}/build/bin/oxDNA input_MD_bn > out &
		${oxdna_path}/build/bin/oxDNA input_MD_un > out &

		if [ $counts -eq ${Ncores} ]; then
			counts=0
			wait
		fi

	done
	
	#integrate step
	counts=0
	for((i=0; i<${Nseqs};i++)); do

		counts=$(($counts+1))

		sim_path=${path}/$(($counts-1))
		
		cd $sim_path
		
		if [ $counts -eq ${Ncores} ]; then
			counts=0
			wait
		fi
		
		python3 ../integrate_step.py ${Ts[$i]} $q -$dq ${Nbs[$i]} &

	done

done

wait

mkdir -p ${path}/Results/${Cs[$i]}/Averages

cd ${path}/Results

for((i=0; i<${Nseqs};i++)); do
	sim_path=${path}/$i
	cp ${sim_path}/integrated_mT_vs_q.dat ${path}/Results/${Cs[$i]}/integrated_mT_vs_q_${i}.dat
	cp ${sim_path}/integrated_averages.dat ${path}/Results/${Cs[$i]}/Averages/integrated_averages_${i}.dat
	
	cd ${path}/Results/${Cs[$i]}
	
	echo "# q_eff mT; seq=${seqs[$i]} cs=${Cs[$i]} T0=${Ts[$i]}" | cat - integrated_mT_vs_q_${i}.dat > temp.txt && mv temp.txt integrated_mT_vs_q_${i}.dat[1][2] #prepend
	cd ${path}/Results/${Cs[$i]}/Averages/
	echo "# q_eff mT <Vdh_bn> <V_bn> <V_dh_un> <V_un>; seq=${seqs[$i]} cs=${Cs[$i]} T0=${Ts[$i]}" | cat - integrated_averages_${i}.dat > temp.txt && mv temp.txt integrated_averages_${i}.dat[1][2]
done

