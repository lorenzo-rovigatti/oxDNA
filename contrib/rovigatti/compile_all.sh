#!/bin/bash

for f in $(ls -1 *cpp)
do
	base=${f%.*}
	echo "Compiling $base"
	bash compile.sh $base
done
