#!/bin/bash

if [ $# -gt 0 ]
then
    filepath="--exec_filepath=$1"
    if [ ! -e $1 ]
    then
	echo "Executable $1 not found. Try to first compile oxDNA and then run tests"
	exit
    fi
else
    filepath=""
fi

python RunTest.py $filepath RUNNER_LIST_BUILD
