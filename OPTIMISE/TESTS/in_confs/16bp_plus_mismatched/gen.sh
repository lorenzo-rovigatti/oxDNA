#!/bin/bash

#random 18-mer double helix with all bp steps plus a mismatched basepair at the 3' end

template_path=/home/andrea/Desktop/oxDNA3_explo_analysis/template
BOX_SIZE=20.

${template_path}/in_confs/16bp_all_steps/GEN ${RANDOM}

#add extra bp
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

tac generated.top | awk '
		NR==1{
			if($2=="A") $2="G"
			else if($2=="G") $2="A"
			else if($2=="T") $2="C"
			else $2="T"
		}1' | tac >tmp.txt
mv tmp.txt generated.top

