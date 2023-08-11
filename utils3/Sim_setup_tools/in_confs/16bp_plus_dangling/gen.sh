#!/bin/bash

#random 16-mer double helix with all bp steps, plus a dangling base at the 3' end

template_path=/home/andrea/Desktop/oxDNA3_explo_analysis/template
BOX_SIZE=20.
Nbp=18 #do not change
Bp="$((${Nbp}*2-1))"

${template_path}/in_confs/16bp_all_steps/GEN ${RANDOM}
rnd=$(( ${RANDOM} % 4 ))

if [[ ${rnd} -eq 0 ]]; then
	sed -i 's/$/A/' seq.txt
elif [[ ${rnd} -eq 1 ]]; then
	sed -i 's/$/C/' seq.txt
elif [[ ${rnd} -eq 2 ]]; then
	sed -i 's/$/G/' seq.txt
else
	sed -i 's/$/T/' seq.txt
fi

python3 ${template_path}/in_confs/generate-sa.py ${BOX_SIZE} seq.txt

sed -i '$ d' generated.top
sed -i '$ d' generated.dat

tac generated.top | awk '
		NR==1{
			$4="-1"
		}1' | tac >tmp.txt
mv tmp.txt generated.top

awk -v Bp=${Bp} '
	NR==1{
		$1=Bp
	}1' generated.top > tmp.txt
	
mv tmp.txt generated.top

