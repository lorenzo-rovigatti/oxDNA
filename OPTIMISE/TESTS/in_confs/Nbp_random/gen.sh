#!/bin/bash

template_path=here
BOX_SIZE=20.
Nbp=12

${template_path}/in_confs/Nbp_random/GEN ${Nbp} ${RANDOM}
python3 ${template_path}/in_confs/generate-sa.py ${BOX_SIZE} seq.txt
