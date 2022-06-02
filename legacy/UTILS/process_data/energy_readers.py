#!/usr/bin/env python27


import base_mod as base
import numpy as np
import os.path
import sys
#import readers 
class LorenzoReader:
    def __init__(self, configuration, topology, check_overlap=False):
        self._conf = False
        
        if not os.path.isfile(configuration):
            base.Logger.die("Configuration file '%s' is not readable" % configuration)

        if not os.path.isfile(topology):
            base.Logger.die("Topology file '%s' is not readable" % topology)
        
        self._check_overlap = check_overlap
        self._conf = open(configuration, "r")

        f = open(topology, "r") 
        f.readline()
        self._top_lines = f.readlines()

    def __del__(self):
        if self._conf: self._conf.close()
        
    def _read(self, only_strand_ends=False, skip=False):
        timeline = self._conf.readline()  
        time = 0.
        if  len(timeline) == 0:
            return False
        else:
            time = float(timeline.split()[2])
		
        box = np.array([float(x) for x in self._conf.readline().split()[2:]])
        [E_tot, E_pot, E_kin] = [float(x) for x in self._conf.readline().split()[2:5]]

        if skip:
            for tl in self._top_lines:
                self._conf.readline()

            return False

        system = base.System(box, time=time, E_pot=E_pot, E_kin=E_kin)
        base.Nucleotide.index = 0
        base.Strand.index = 0 
        s = base.Strand()

        for tl in self._top_lines:
            tls = tl.split()
            n3 = int(tls[2])
            n5 = int(tls[3])
            b = base.base_to_number[tls[1]]
            ls = self._conf.readline().split()
            cm = [float(x) for x in ls[0:3]]
            a1 = [float(x) for x in ls[3:6]]
            a3 = [float(x) for x in ls[6:9]]
            v = [float(x) for x in ls[9:12]]
            L = [float(x) for x in ls[12:15]]
            if not only_strand_ends or n3 == -1 or n5 == -1:
                s.add_nucleotide(base.Nucleotide(cm , a1, a3, b, v, L))

            if n5 == -1:
                system.add_strand(s, self._check_overlap)
                s = base.Strand()

        return system
        

    # if only_strand_ends == True then a strand will contain only the first and the last nucleotide
    # useful for some analysis like csd for coaxial interactions
    def get_system(self, only_strand_ends=False, N_skip=0):
        for i in range(N_skip):
            self._read(skip=True)
 
        return self._read(only_strand_ends=only_strand_ends, skip=False)




print "Starting" 

if (len(sys.argv) < 3):
  print 'Usage %s input_file configuration_file' % sys.argv[0]
  sys.exit()
  
#now get topology file name:
inputfile = sys.argv[1]
conffile = sys.argv[2]
topologyfile = ""
fin = open(inputfile)
for line in fin:
  if "topology" in line:
    topologyfile = line.split('=')[1].replace(' ','').replace('\n','')

myreader = LorenzoReader(conffile,topologyfile)
mysystem = myreader.get_system()
mysystem.map_nucleotides_to_strands()

if not os.path.isfile("output_bonds"):
	print "Cannot execute output_bonds program. Please make sure to go to process_data/ directory and type make"
	sys.exit(1)

myinput = os.popen('./output_bonds ' + str(inputfile ) + " " + ' 2>/dev/null' )
mysystem.read_H_bonds(myinput.readlines())
myinput.close()
mysystem.show_H_interactions()
