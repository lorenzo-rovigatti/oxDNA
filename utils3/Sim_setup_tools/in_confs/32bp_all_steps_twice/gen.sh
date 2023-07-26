#!/bin/bash

template_path=/home/andrea/Desktop/oxDNA3_explo_analysis/template
BOX_SIZE=20.

${template_path}/in_confs/32bp_all_steps_twice/GEN ${RANDOM}
python3 ${template_path}/in_confs/generate-sa.py ${BOX_SIZE} seq.txt
