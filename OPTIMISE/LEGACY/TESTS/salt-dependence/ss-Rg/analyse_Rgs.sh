path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

Ns=(20 40 60 80 100)
Cs=(0.1 0.15 0.2 0.3 0.5 1.1)

Nlens=${#Ns[@]}
NCs=${#Cs[@]}
Ncores=12

counts=0
for((i = 0; i < ${NCs}; i++)); do
	for((j = 0; j < ${Nlens}; j++)); do
		echo "Analysing rep: "${counts}
		counts=$(($counts+1))
		
		sim_path=${path}/${i}
		
		cd ${sim_path}
			
		python3 ../compute_Rg_ss.py > out_analysis &
		
		if [ $counts -eq ${Ncores} ]; then
			wait
		fi
	done
done

wait
cd ${path}
python3 plot_Rg_ete.py
