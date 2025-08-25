path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

Nseqs=10
Ncores=10

counts=0
for((i = 0; i < ${Nseqs}; i++)); do
	echo "Analysing rep: "${counts}
	counts=$(($counts+1))
	
	sim_path=${path}/${i}
	
	cd ${sim_path}
		
	python3 ../plot_internal_w_comparison.py trajectory.dat generated.top ${i} > out_analysis &
	
	if [ $counts -eq ${Ncores} ]; then
		wait
	fi
	
done

mkdir -p ${path}/Plots

for((i = 0; i < ${Nseqs}; i++)); do
	sim_path=${path}/${i}
	cp ${sim_path}/*.pdf ${path}/Plots
	
done
