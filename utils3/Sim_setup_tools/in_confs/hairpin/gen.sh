#!/bin/bash

template_path=here
TYPE=pur_pur

cp ${template_path}/in_confs/hairpin/generated_hairpin.dat .
cp ${template_path}/in_confs/hairpin/generated_hairpin_${TYPE}.top .
mv generated_hairpin.dat generated.dat
mv generated_hairpin_${TYPE}.top generated.top
