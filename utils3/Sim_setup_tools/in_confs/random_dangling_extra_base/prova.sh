#!/bin/bash

Nbp=12
Bp="$((${Nbp}*2-1))"

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
