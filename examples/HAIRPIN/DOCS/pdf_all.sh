#!/bin/bash

for f in $(ls *agr)
do
	base=${f%.*}
	xmgrace -hardcopy -hdevice EPS $f -printfile $base.eps
	epstopdf $base.eps --outfile=$base.pdf
	rm $base.eps
done
