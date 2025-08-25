path=$(pwd)
oxdna_path=/home/users/yqb22156/oxdna

Nseqs=6
Ncores=12

counts=0
for((i = 0; i < ${Nseqs}; i++)); do
	echo "Analysing rep: "${counts}
	counts=$(($counts+2))
	
	sim_path=${path}/${i}
	
	cd ${sim_path}
		
	python3 ../compute_elastic_moduli.py > out_analysis &
	python3 ../compute_pitch.py generated.top trajectory.dat 20 > out_pitch &
	
	if [ $counts -eq ${Ncores} ]; then
		wait
	fi
	
done

echo "Analysis done"
echo "Remember to create a file with the limit values of the persistence lengths, before running plot_lps.py"
echo "These values can be found in the out_analysis files inside the sim directories"
echo "This file should be called limt_lps.txt"
echo "The structure should be: salt_c lb err lt/2 err"
