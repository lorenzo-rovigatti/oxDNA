path=$(pwd)

counts=0
for((i = 0; i < 5; i++)); do
	for((j = 0; j < 5; j++)); do
	if [ $i -eq 0 ] && [ $j -eq 0 ]; then
		continue
	fi
	
	echo ${counts}
	counts=$(($counts+1))
	sim_path=${path}/${i}${j}
	cd ${sim_path}

	python3 compute_pitch.py generated.top trajectory.dat 20 > out_pitch &
	
	if [ $counts -eq 12 ]; then
		wait
	fi
	
	done
done
