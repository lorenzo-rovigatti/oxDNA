#!/bin/bash

CODEDIR=../..
OXDNA_BIN=$CODEDIR/bin/oxDNA

# generate the initial configuration files (initial.conf and initial.top)
echo GCGTTGCTTCTCCAACGC > sequence.sqs
python $CODEDIR/UTILS/generate-sa.py 50 sequence.sqs
mv generated.dat initial.conf
mv generated.top initial.top

if [ $# -gt 0 ]
then
    if [ $1 = "--generate-only" ]
    then
	echo "Configuration generated, exiting"
	exit
    fi
fi

echo "Simulating a 18-base hairpin (6-base stem + 6-base loop) at 334 K"

if [ -e $OXDNA_BIN ]
then
    # prepare the inputMD_seq_dep file
    sed "s;seq_dep_file = ../../sequence_dependent_parameters.txt;seq_dep_file = $CODEDIR\/sequence_dependent_parameters.txt;" inputMD_seq_dep > input.tmp
    mv -f input.tmp inputMD_seq_dep

    # run the examples
    echo "Starting MD simulation with the sequence-averaged version of the model"
    $OXDNA_BIN inputMD
    
    echo "Starting MD simulation with the sequence-dependent version of the model"
    $OXDNA_BIN inputMD_seq_dep

    echo "Starting MC simulation with the sequence-averaged version of the model and traps acting between nucleotides (see hairpin_forces.dat for details of the traps)"
    $OXDNA_BIN inputTRAP
else
    echo "Can't find $OXDNA_BIN, did you compile oxDNA?"
    exit
fi
