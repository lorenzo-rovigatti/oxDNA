#!/usr/bin/env python
'''
Prints a topology file of a double strand in which each base is a dummy that can hydrogen bind only the
base opposite to it.
'''
import sys
usage="USAGE:"

def generate_dsDNA_nomismatch_top( file_name, nbp ):
  '''Generate a topology file for a dsDNA in which there can be no mismatch'''
  if nbp < 1:
    print 'ERROR: in function generate_dsDNA_nomismatch_top, nbp =',nbp,'but it should be nbp > 0'
    sys.exit(-1)
  try:
    f = open(file_name,'w')
  except:
    print 'ERROR: In function generate_dsDNA_nomismatch_top, could not open file',file_name,'in write mode.'
    sys.exit(-1)
  print >> f ,2*nbp,2
  for i in range(0,nbp-1):
    print >> f, 1, 13 + i, i - 1, i + 1 
  print >> f, 1, 13 + nbp - 1, nbp - 2, -1
  print >> f, 2, -10 - nbp + 1, -1, nbp + 1 
  for i in range(1,nbp-1):
    print >> f, 2, -10 -nbp + 1 + i, nbp + i - 1, nbp + i + 1 
  print >> f, 2, -10, nbp*2-2, -1
  f.close()

out_topology_filepath = 'prova.top'
if len(sys.argv) < 2:
	print __doc__ 
	print "USAGE:",sys.argv[0]," <N_bp> [output_file (defaults to", out_topology_filepath,")]"
N = int(sys.argv[1])
if len(sys.argv) >2:
	out_topology_filepath = sys.argv[2]
generate_dsDNA_nomismatch_top(out_topology_filepath,N)
