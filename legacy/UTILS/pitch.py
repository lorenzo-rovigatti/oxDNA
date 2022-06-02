#!/usr/bin/env python

# find pitch of a double strand
# WARNING: assumes that 
import readers
import base
import sys
try:
    import numpy as np
except:
    import mynumpy as np
import utils
import os
import pickle
import subprocess

# We ignore a quantity TRIM nucleotides at the beginning and end of the strand when doing this analysis
TRIM = 5
PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")
CONF_SKIP = 0
DEBUG = False
CROSS_RANGE = 200

def find_pitch_angle(system, nuc1i, nuc3i, interaction_list):
    nuc2i = interaction_list[nuc1i]
    nuc4i = interaction_list[nuc3i]
    nuc1 = system._nucleotides[nuc1i]
    nuc2 = system._nucleotides[nuc2i]
    nuc3 = system._nucleotides[nuc3i]
    nuc4 = system._nucleotides[nuc4i]
    r1 = nuc1.get_pos_base()
    r2 = nuc2.get_pos_base()
    r3 = nuc3.get_pos_base()
    r4 = nuc4.get_pos_base()
    n = (r1+r2)/2 - (r3+r4)/2
    n /= np.sqrt(np.dot(n,n))

    b1 = nuc1.get_pos_back()
    b2 = nuc2.get_pos_back()
    b3 = nuc3.get_pos_back()
    b4 = nuc4.get_pos_back()
    l1 = b1 - b2
    l2 = b3 - b4
    l1prime = l1 - np.dot(l1, n) * n
    l1prime /= np.sqrt(np.dot(l1prime,l1prime))
    l2prime = l2 - np.dot(l2, n) * n
    l2prime /= np.sqrt(np.dot(l2prime,l2prime))

    theta = np.arccos(np.dot(l1prime, l2prime))
    return theta

if len(sys.argv) < 5:
    base.Logger.log("Usage is %s configuration topology input_file outfile [-up_to quantity] [-repeat number]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

up_to_conf = False
repeat_num = False
if len(sys.argv) > 3:
    if "-up_to" in sys.argv:
        up_to_conf = int(sys.argv[sys.argv.index("-up_to") + 1])

    if "-repeat" in sys.argv:
        repeat_num = int(sys.argv[sys.argv.index("-repeat") + 1])

outfile = sys.argv[4]

l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()

if CONF_SKIP:
    base.Logger.log("skipping %d configurations" % CONF_SKIP, base.Logger.INFO)
while CONF_SKIP:
    s = l.get_system()
    CONF_SKIP -= 1
    
s.map_nucleotides_to_strands()
infile = sys.argv[3]
try:
    open(infile)
except:
    base.Logger.log("unable to find file %s, exit" % infile, base.Logger.CRITICAL)
    sys.exit()
launchargs = [PROCESSDIR + 'output_bonds',infile,sys.argv[1],str(CONF_SKIP)]
myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
s.read_H_bonds(myinput.stdout.readlines())
interaction_list = [-1 for i in range(s._N)]
for strand in s._strands:
    for nucleotide in strand._nucleotides:
        for i in nucleotide.interactions:
            interaction_list[nucleotide.index] = i

counter = 0
for i in interaction_list:
    if i == -1:
        counter += 1
if counter > 0:
    base.Logger.log("%d unhybridised nucleotides found, is this intended?" % counter, base.Logger.WARNING)

base.Logger.log("Using trim %s" % (TRIM), base.Logger.INFO)

w_len = s._N / 2 - 2 * TRIM - 1 # length of the part of the strand we will work with
angles = [0 for x in range(w_len)] # angle averaged over all configs for a particular pair of adjacent nucleotide pairs (npp)
conf_counter = CONF_SKIP
conf_av_angle = [] # angle averaged over all relevant pairs for a particular configuration
while s:
    conf_av_angle.append(0)
    base.Logger.log("Working on conf %s" % str(conf_counter + 1), base.Logger.INFO)

    for i in range(w_len):
        this_angle = find_pitch_angle(s, i + TRIM, i + TRIM + 1, interaction_list)
        angles[i] += this_angle
        conf_av_angle[-1] += this_angle

    s = l.get_system()
    conf_counter += 1
    if up_to_conf and conf_counter >= up_to_conf:
        break

# find statistics
confs_used = conf_counter - CONF_SKIP

angles = [x/confs_used for x in angles]
av_angle = sum(angles) / len(angles)

av_angle2 = 0
conf_av_angle = [x/w_len for x in conf_av_angle]
for angle in conf_av_angle:
    av_angle2 += angle*angle # sum of angle^2 where angle is for averaged over npp for a particular configuration
av_angle2 /= confs_used

cc_i = [0 for i in range(CROSS_RANGE)]
for i in range(CROSS_RANGE):
    for j in range(confs_used - CROSS_RANGE):
        cc_i[i] += conf_av_angle[j] * conf_av_angle[j + i] - av_angle*av_angle

standard_deviation = np.sqrt(av_angle2 - av_angle*av_angle)
cc_i = [x/((confs_used-CROSS_RANGE)*standard_deviation*standard_deviation) for x in cc_i]

pitch = 2 * np.pi / av_angle

# print out everything
r2d = 180/np.pi
print "pitch: %f bp/turn" % pitch
print av_angle * r2d, av_angle2 * r2d*r2d, standard_deviation, standard_deviation/np.sqrt(confs_used)

f_out = open("angle_series.dat", "w")
for i in range(len(conf_av_angle)):
    f_out.write("%f\n" % conf_av_angle[i])

f_out.close()
f_out = open(outfile, "w")
for i in range(len(angles)):
    f_out.write("%f\n" % angles[i])

f_out.close()
f_out = open("autocor.dat", "w")
for i in range(len(cc_i)):
    f_out.write("%f\n" % cc_i[i])
