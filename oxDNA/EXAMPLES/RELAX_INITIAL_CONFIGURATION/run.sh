#!/bin/sh

# run the relaxation procedure

../../bin/oxDNA input_weak_constant
mv last_conf.dat relaxed_weak.conf
../../bin/oxDNA input_strong_constant
mv last_conf.dat relaxed.conf
