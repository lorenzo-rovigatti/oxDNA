path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

#Ns = [20, 40, 60, 80, 100]

Nseqs=5
Ncores=5

counts=0
for((i = 0; i < ${Nseqs}; i++)); do
	echo "Analysing rep: "${counts}
	counts=$(($counts+1))
	
	sim_path=${path}/${i}
	
	cd ${sim_path}
		
	python3 ../compute_Rg_ss.py > out_analysis &
	
	if [ $counts -eq ${Ncores} ]; then
		wait
	fi
	
done

cd ${path}
python3 plot_Rg_ete.py
