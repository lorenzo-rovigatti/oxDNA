#!/bin/bash

#random N-mer double helix with a dangling base at the 3' end

template_path=/home/andrea/Desktop/oxDNA3_explo_analysis/template
BOX_SIZE=20.
Nbp=12 #nuber of bp of the double helix+1 (we remove a base from a double helix with Nbp base pairs)
Bp="$((${Nbp}*2-1))"

${template_path}/in_confs/random_mismatched_bp/GEN ${Nbp} ${RANDOM}
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

