#!/bin/bash

if [ $# -lt 1 ]
then
	echo "Usage is $0 name [debug]"
	exit
fi

f=$1

compilation_flags="-Wall -fpic -c $f.cpp -I../../src/Observables -I../../src/Interactions"
linking_flags="-shared -o $f.so $f.o"
if [ $# -eq 1 ]
then
	g++ -O3 $compilation_flags
	g++ $linking_flags
else
	g++ -g2 $compilation_flags
	g++ -g2 $linking_flags
fi
