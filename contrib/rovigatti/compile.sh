#!/bin/bash

if [ $# -lt 1 ]
then
	echo "Usage is $0 name [debug]"
	exit
fi

f=$1


if [ $# -eq 1 ]
then
	g++ -O3 -fpic -c $f.cpp -I../../src/Observables -I../../src/Interactions
	g++ -shared -o $f.so $f.o
else
	g++ -g2 -fpic -c $f.cpp -I../../src/Observables -I../../src/Interactions
	g++ -g2 -shared -o $f.so $f.o
fi
