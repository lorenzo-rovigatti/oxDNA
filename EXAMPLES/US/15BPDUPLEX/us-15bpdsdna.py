#!/usr/bin/env python
#
# Example file for the set up of an umbrella sampling simulation of the free
# energy profile for bonding (number of hydrogen bonds >= 1) for a 15bp dsDNA at
# T=343K (Figure 6 in Ouldridge et. al., JCP, 2011). This setup does NOT include
# the first jump. It takes a couple of days on a single CPU to get reliable
# results, maybe a week to get publication quality results.
#

import sys

# if you don't have the base.py and readers.py in your PYTHONPATH,
# the following will work only if you move this file in the same
# directory as base.py and readers.py
import base, readers
import order_parameters as op
from umbrella_sampling import *

# we need an initial system to build the order parameters on
system = op.load_system ('input', 'inizio.dat', 'prova.top')

# op.complementary_bonds_between_strands in an helper function to build a
# dictionary
mydict = op.complementary_bonds_between_strands(system, 'COR_HB_0_1', 0, 1, start=0, include_all_parameters=False)

# we use the dictionary to build the order parameter myop
myop = op.OrderParameters (mydict)

# simple array of weights
myarr = [1. for i in xrange (0, 16)]

# helper functions to define weights: one from an array, one from an
# evaluateable string
mydict1 = op.define_weight ('W1', 'COR_HB_0_1', myarr)
mydict2 = op.define_weight ('W2', '$(COR_HB_0_1) >= 1', myop)

# we build the Weights objects to use in the umbrella sampling. myw is the
# weight that the simulation actually uses, myforb will be useful later for the
# adaptive weights. Please note that myw contains both dictionaries just built
# with the above helper functions
myw = op.Weights (mydict1, mydict2)
myforb = op.Weights (mydict2)

# Here we build the UmbrellaSampling object. It requires a directory, an
# executable, an imput file, etc.
sim = UmbrellaSampling (
    dir = 'SIM', # directory
    exe = './CUDADNA',
    input = 'input',
    ordparams = myop,
    weights = myw,
    debug = False,
    output_interval = 10
    )

# If something goes wrong, the first thing to check is if the simulation
# actually runs in the directory. run_test() is useful to run a quick test
# sim.run_test ()

# Here we run the simulation.
outer = 10000
inner = 1000
print "STARTING ..."
for i in xrange (outer):
    for dummy in xrange (inner):
        # run the given number of steps
        sim.run (20)

        # print information about the umbrella sampling object. Useful to check
        # for exchange rates, acceptance rate of trajectories, python overhead
        sim.print_info ()

        # we print the histogram both to a file and to the screen
        # columns are 
        # [value of order parameters 1, ..., N], biased prob, unbiased_prob
        # so you will have M + 2 columns, where M is the number of your order
        # parameters (1 in this case)
        sim.print_histo ("15bphisto.dat")
        sim.print_histo ()
        print

    # every time the inner cycle is finished, we automatically adjust the weights
    # given one histogram. op.adaptive_weights returns a Weights instance where
    # the weights are the inverse of the probabilities found in the histogram.
    # Note that:
    #    a) the weights are not changed if the state has not been visited
    #    b) at the maximum, the weight is changed by a factor 2.
    #    c) the last argument is a list of "carried over" weights, that we want
    #       to keep active. It can be empty. It is useful to keep, as done in
    #       this example, a weight that accounts only for forbidden states.
    print '--> changing weights...'
    myw = op.adaptive_weights (sim.histo, myw, [myforb])

    # you may or may not want to reset the histogram. Please note that data
    # gathered before will be accounted for unless you do a histo.deep_reset()
    # it is useful to reset the histogram to have a feeling of how many states
    # the simulation has visited between each weight update.
    sim.histo.reset ()

# speaks for itself :)
print "FINISHED"

# Flavio

