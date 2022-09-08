#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import fileinput
import getopt

input_options = ['seq_base_colour','colour_by_domain=']

longArgs=input_options
shortArgs='sc:'

try:
    args, files = getopt.getopt(sys.argv[1:], shortArgs, longArgs)
except:
    base.Logger.log( "Wrong usage. Aborting",base.Logger.CRITICAL)
    sys.exit (-2)


if len(files) < 2:
    base.Logger.log("Usage is %s [-s] [-c domain_file] configuration topology [output]" % sys.argv[0] +"\n-s colours bases by their type\n-c allows colouring of backbones by domain\ndomain_file should assign domains to all bases, in the format:\n#a,#b,#c,#d,\n#e,#f,#g,#h,\netc\nwhere strand #a, bases #b->#c have a domain value of #d;\nstrand #e, bases #f->#g have a domain value of #h etc.\nAll bases of a certain domain will be the same colour, more than\none region of DNA can be tagged as a certain domain (eg., #d=#h is allowed).\nTo label two domains #d and #h as complementary, set #h=-#d.\nSimilar colours will be used for two complementary domains.\nThe system can handle 8 distinctly coloured domains and 8 complementary domains,\nusing any more will lead to repeated colours.\nAny base not assigned to a domain will be coloured sandy brown.", base.Logger.CRITICAL)
    sys.exit()

if len(files) >2:
    output = files[2]
else: output = files[0] + ".pdb"

l = readers.LorenzoReader(files[0], files[1])
s = l.get_system()
cdm20 = True
append = False
domain=[]
domains=False
domain_file=''
colour_by_seq=False
for opt,arg in args:
    if opt in('-c','colour_by_domain'):
        domains=True
        domain_file=arg
    elif opt in('-s','seq_base_colour'):
        colour_by_seq=True
    else:
        base.Logger.log("Option '%s' not recognised" %opt, base.Logger.INFO)    

if domains == True:
    for i in range(len(s._strands)):
        new=[]
        for j in range(len(s._strands[i]._nucleotides)):
            new.append(0)
        domain.append(new)
        
    
    fi = fileinput.FileInput([domain_file])
    try:
        for line in fi:
            d=line.split(",")
            #print d # WDEBUG
            strandid = int(d[0])
            firstbid =int(d[1])
            finalbid=int(d[2])
            domid=int(d[3])
            #print strandid, firstbid, finalbid, domid # WDEBUG    - Ben and Will fixed an error 29/05/12 - incorrect Logger syntax
            #print len(s._strands) # WDEBUG 
            if strandid >= len(s._strands):
                base.Logger.log ("strand id %d is too large: domain ignored."%strandid, base.Logger.INFO)
            elif (firstbid <0 or firstbid>= len(s._strands[strandid]._nucleotides)):
                base.Logger.log ("base id %d is not allowed for strand %d: domain ignored."%(firstbid, strandid), base.Logger.INFO)
            elif (finalbid <0 or finalbid>= len(s._strands[strandid]._nucleotides)):
                base.Logger.log ("base id %d is not allowed for strand %d: domain ignored."%(finalbid, strandid), base.Logger.INFO)
            else:
                for i in range(finalbid-firstbid+1):
                    domain[strandid][firstbid+i] = domid
    except:
        base.Logger.log( "Could not parse domains file %s, aborting" %domain_file,base.Logger.CRITICAL)
        sys.exit(-3)
    #print str(domain)
while s:
    #s.bring_in_box_nucleotides()
    # azzeriamo il centro di massa
    if cdm20:
        base.Logger.log ("tracking particle 0", base.Logger.INFO)
        cdm = (s._strands[0]._nucleotides[0].cm_pos).copy()
        for strand in s._strands:
            strand.translate (-cdm)
        base.Logger.log ("setting cdm to 0", base.Logger.INFO)
        cdm = np.array ([0.,0.,0.])
        for strand in s._strands:
            for n in strand._nucleotides:
                cdm += n.cm_pos
        cdm = np.array ([0.,0.,0.])
        cdm = cdm / float (s.get_N_Nucleotides())
        for strand in s._strands:
            strand.translate (-cdm)

    s.print_pdb_output_chimera(output, append=append, domain=domain,colour_by_seq=colour_by_seq)
    s = l.get_system()
    append = True

base.Logger.log("Output printed on '%s'" % output, base.Logger.INFO)
base.Logger.log("Chimera instructions printed to ./chimera.com '%s'" % output, base.Logger.INFO)
