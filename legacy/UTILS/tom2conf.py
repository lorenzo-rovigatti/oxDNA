#!/usr/bin/env python

import sys
try:
    import numpy as np
except:
    import mynumpy as np
import base
import generators as gen

def read_tom_strand_from_string(strandline,box):
	dvals = strandline.split('(')
	if(len(dvals) != 2):
		raise IOError
	vals = dvals[1].replace('\"',' ').split(',')
	length = int(vals[0])
	s = base.Strand()
	valindex = 1
	for i in range(length):
		#valindex += i*9+1
		posback = []
		posbackbase = []
		posalign = []
		
		nuctype = int(vals[valindex])
		valindex += 1
		for j in range(3):
			posback.append(float(vals[valindex]))
			valindex  += 1
		for j in range(3):
			posbackbase.append(float(vals[valindex]))
			valindex += 1
		for j in range(3):
			posalign.append(float(vals[valindex]))
			valindex += 1
		poscm = np.array(posback) - np.array(posbackbase) *  base.POS_BACK
                #poscm = poscm % box
		n = base.Nucleotide(poscm,posbackbase,posalign,nuctype)
		s.add_nucleotide(n)
	return s	

def read_tom_strands_from_file(filename,box = 20):
	infile = open(filename,'r')
	mysys = base.System([box,box,box])
	for line in infile.readlines():
		if len(line) > 2:
			strand = read_tom_strand_from_string(line,box)
			mysys.add_strand(strand, check_overlap=False)
	
	return mysys
	

if (len(sys.argv) < 2):
	print "Usage: %s <file_with_strands_in_tom_format> [boxsize]" % (sys.argv[0])
	print "File with strands has to have one strand per line in .cpp format (DNA strand(....);) as produced by the Tom's code. Boxsize has to be set to the cell size, by default it is 20"
        sys.exit(1)
    
fname = sys.argv[1]    

box = 20
if len(sys.argv) >= 3:
    box = float(sys.argv[2])
    print >> sys.stderr, "## using box size =", box

s = read_tom_strands_from_file(fname,box)
s.print_lorenzo_output(fname+".conf", fname+".top")

print "Loaded ", s.get_N_strands(), " strands"
print "Output written to " , fname + ".conf", " ", fname + ".top"

'''else:
    if sys.argv[1] == "ds": double_strands()
    elif sys.argv[1] == "t":  tetramers()
    else:
        base.Logger.die("Unrecognized option '%s'" % sys.argv[1])'''

