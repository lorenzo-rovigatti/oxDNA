#!/bin/bash

CODEDIR=../build
OXDNA_BIN=$CODEDIR/bin/oxDNA

if [ -e $OXDNA_BIN ]
then
    # prepare the inputMD_seq_dep file
    # sed "s;seq_dep_file = ../oxDNA1_sequence_dependent_parameters.txt;seq_dep_file = $CODEDIR\/oxDNA1_sequence_dependent_parameters.txt;" input_seq_dep > input.tmp
    # mv -f input.tmp input_seq_dep

    # run the samples
    echo "Starting VMMC simulation with the sequence-averaged version of the model"
    $OXDNA_BIN input
    
    echo "Starting VMMC simulation with the sequence-dependent version of the model"
    $OXDNA_BIN input_seq_dep

    echo "Starting VMMC simulation with the sequence-averaged version of the model and traps acting between nucleotides (see hairpin_forces.dat for details of the traps)"
    $OXDNA_BIN input_trap
else
    echo "Can't find $OXDNA_BIN, did you compile oxDNA?"
    exit
fi
