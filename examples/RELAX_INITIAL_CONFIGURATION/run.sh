#!/bin/sh

# run the relaxation procedure

../../build/bin/oxDNA input_weak_constant
mv last_conf.dat relaxed_weak.conf
../../build/bin/oxDNA input_strong_constant
mv last_conf.dat relaxed.conf
