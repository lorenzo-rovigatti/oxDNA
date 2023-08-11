#!/bin/bash

template_path=here
BOX_SIZE=20.
FILE=file

cp ${template_path}/in_confs/examples/${FILE} .
python3 ${template_path}/in_confs/generate-sa.py ${BOX_SIZE} ${FILE}
