#!/bin/bash

#random N-mer double helix with a mismatched basepair at the 3' end

template_path=here
BOX_SIZE=20.
Nbp=12

${template_path}/in_confs/random_mismatched_bp/GEN ${Nbp} ${RANDOM}
python3 ${template_path}/in_confs/generate-sa.py ${BOX_SIZE} seq.txt

tac generated.top | awk '
		NR==1{
			if($2=="A") $2="G"
			else if($2=="G") $2="A"
			else if($2=="T") $2="C"
			else $2="T"
		}1' | tac >tmp.txt
mv tmp.txt generated.top

