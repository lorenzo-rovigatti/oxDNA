#!/usr/bin/env python

# detect weave pattern in 2D origami
# NB consider what happens when bp are not hybridised

import readers
import base
import sys
import numpy as np
import subprocess
import os
import origami_utils as oru

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")
CONF_SKIP = 0
IGNORE_UNBONDED = True # ignore base pair that aren't hydrogen bonded - if you have a nice starting configuration with all base pairs hybridised you can set to false to speed up the program. Otherwise you can set to true and the program refreshes the interaction list for every configuration.

def get_bb_midpoint(system, strand, n_index):
    # get midpoint vector between 2 hybridised bases
    r1 = strand._nucleotides[n_index].get_pos_base()
    r2 = system._nucleotides[interaction_list[strand._nucleotides[n_index].index]].get_pos_base()
    vec = (r1+r2)/2
    return vec

if (len(sys.argv) < 5):
  print 'Usage %s trajectory topology input output' % sys.argv[0]
  sys.exit()

conffile = sys.argv[1]
topologyfile = sys.argv[2]
infile = sys.argv[3]
outfile = sys.argv[4]
myreader = readers.LorenzoReader(conffile,topologyfile)
s = myreader.get_system()
scaf_index = oru.get_scaffold_index(s)

if not os.path.isfile(PROCESSDIR + "output_bonds"):
	print "Cannot execute output_bonds program. Please make sure to go to process_data/ directory and type make"
	sys.exit(1)

if not os.path.isfile(infile):
    base.Logger.log("unable to open file %s" % infile, base.Logger.CRITICAL)
    sys.exit()
    
# skip configurations
for counter in range(CONF_SKIP):
    s = myreader.get_system()

# get vhelices data from file - format is either <auto,[number of vhelices]> or, if a region of an origami is to be analysed, <[region width],[list of starting nucleotide index for the region for each vhelix]>
vhelix_def_file = open("data.vhd", "r")
data = [x for x in vhelix_def_file.readline().replace("\n","").replace(" ","").split(",")]
if data[0] == "auto":
    num_vh = int(data[1])
    origami_width = s._strands[scaf_index].get_length() / num_vh
    vhelix_indices = []
    start_nuc_ind = -1
    for i in range(num_vh):
        if i % 2 == 0:
            start_nuc_ind += 1
        else:
            start_nuc_ind += origami_width*2 - 1
        vhelix_indices.append(start_nuc_ind)
else:
    origami_width = int(data[0])
    vhelix_indices = [int(x) for x in data[1:]]
    num_vh = len(vhelix_indices)
base.Logger.log("using file data.vhd, %d virtual helices found" % num_vh, base.Logger.INFO)

dist = [[0 for y in range(origami_width)] for x in range(num_vh-1)]
data_counter = [[0 for y in range(origami_width)] for x in range(num_vh-1)]
conf_counter = 0
while s != False:
    if IGNORE_UNBONDED or conf_counter == 0:
        launchargs = [PROCESSDIR + 'output_bonds',infile,conffile,str(CONF_SKIP+conf_counter)]
        s.map_nucleotides_to_strands()
        myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        s.read_H_bonds(myinput.stdout.readlines())
        interaction_list = [-1 for i in range(s._N)]
        for strand in s._strands:
            for nucleotide in strand._nucleotides:
                for i in nucleotide.interactions:
                    interaction_list[nucleotide.index] = i
        if conf_counter == 0:
            counter = 0
            for i in interaction_list:
                if i == -1:
                    counter += 1
            if counter > 0:
                base.Logger.log("%d unhybridised nucleotides found, is this intended?" % counter, base.Logger.WARNING)
    for i in range(num_vh-1):
        for nucleotide in range(origami_width):
            if i % 2 == 0:
                dir = 1
            else:
                dir = -1
            n1 = vhelix_indices[i] + nucleotide * dir
            n2 = vhelix_indices[i+1] + nucleotide * dir * -1
            if (interaction_list[n1] != -1 and interaction_list[n2] != -1) or not IGNORE_UNBONDED: # can choose to use only h-bonded sites
                try:
                    vecdist = get_bb_midpoint(s, s._strands[scaf_index], n1) - get_bb_midpoint(s, s._strands[scaf_index], n2)
                except IndexError:
                    base.Logger.log("IndexError while getting distance between bps - debug output:", base.Logger.CRITICAL)
                    print i, nucleotide, dir, vhelix_indices[i], vhelix_indices[i+1]
                    print n1, n2
                    sys.exit()
                dist[i][nucleotide] += np.sqrt(np.dot(vecdist,vecdist))
                data_counter[i][nucleotide] += 1

    conf_counter += 1
    s = myreader.get_system()

fout = open(outfile, "w")
for j in range(len(dist[0])):
    for i in range(len(dist)):
        if data_counter[i][j] != conf_counter:
            print data_counter[i][j], i, j
        if data_counter[i][j] != 0:
            out_ij = dist[i][j]/data_counter[i][j]
        else:
            out_ij = dist[i][j]
        fout.write(str(out_ij) + " ")
    fout.write("\n")
fout.close()
