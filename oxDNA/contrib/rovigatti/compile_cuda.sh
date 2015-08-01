#!/bin/bash

if [ $# -lt 1 ]
then
	echo "Usage is $0 name [debug]"
	exit
fi

f=$1

compilation_flags="-arch=sm_20 -Xcompiler -fpic -c CUDA$f.cu -I../../src/CUDA/Interactions -I../../src/Interactions"
linking_flags="-shared -o CUDA$f.so CUDA$f.o $f.o -lgsl"
if [ $# -eq 1 ]
then
	nvcc --use_fast_math -O3 $compilation_flags
	g++ $linking_flags
else
	nvcc -O3 --use_fast_math -g $compilation_flags
	g++ -O3 -ffast-math -g2 $linking_flags
fi
