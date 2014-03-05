#!/bin/sh

previous_grooving=$OXDNA_GROOVE
export OXDNA_GROOVE=1
../../UTILS/traj2chimera.py last_conf.dat generated.top
export OXDNA_GROOVE=$previous_grooving
