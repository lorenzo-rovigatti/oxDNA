#!/bin/bash

if [ $# -ne 1 ]
then
	echo "Usage is $0 obs_name"
	exit
fi

f=$1

g++ -fpic -c $f.cpp -I../../src/Observables
g++ -shared -o $f.so $f.o
