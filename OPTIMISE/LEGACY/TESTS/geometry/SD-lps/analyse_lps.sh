path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

poly_steps=('AA' 'AG' 'AC' 'AT' 'CG' 'CC')
richness=(10 30 50 70 90)

Ncores=11

counts=0
for((i = 0; i < 11; i++)); do
	echo "Analysing rep: "${counts}
	counts=$(($counts+1))
	
	sim_path=${path}/${i}
	
	cd ${sim_path}
		
	python3 ../compute_elastic_moduli.py > out_analysis &
	
	if [ $counts -eq ${Ncores} ]; then
		wait
	fi
	
done

cd ${path}
python3 plot_lps.py
