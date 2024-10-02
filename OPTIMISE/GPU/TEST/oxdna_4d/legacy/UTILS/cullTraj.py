#! /usr/bin/env python

"""
This program reads an oxDNA trajectory, and, if the file NO_DISCARD
is absent in the directory, only keeps one trajectory every N. This is
used to save space.
"""
import sys
import os.path

NO_DISCARD_filename = "NO_DISCARD"
# handle arguments and print USAGE if the syntax is wrong.
if len(sys.argv) < 3:
	print "USAGE:",sys.argv[0]," <traj_file> <N>"
	sys.exit(0)
if len(sys.argv) >1:
	traj_filename = sys.argv[1]
	N = int(sys.argv[2])
	if N < 2:
		print "(cullTraj.py) ERROR: N was",N,"but this should be 2 or greater."
		sys.exit(-1)

	
# TODO: check that no NO_DISCARD file is present and that a new one can be written.
if os.path.isfile(NO_DISCARD_filename):
	print "WARNING: The",NO_DISCARD_filename,"file is already in the directory, meaning that the trajectory"
	print "         has already been culled. If you want to cull it further just"
	print "         remove it and launch the script again."
	print "         No action has been performed."
	sys.exit(-1)
try:
	ff = open(NO_DISCARD_filename,'w')
except: 
	print "ERROR: I can't write the file",NO_DISCARD_filename,"."
	print "       I won't do anything unless I can write the file."
	sys.exit(-1)
	

 # read the file line by line
try:
  f = open(traj_filename,'r')
except:
  print "(cullTraj.py) ERROR: trajectory file",traj_filename," could not be opened."
  sys.exit(-1)
traj_by_lines = f.readlines()
f.close()
# count the number of configurations by counting the lines starting in t
N_confs = 0 
for i,l in enumerate(traj_by_lines):
	if l[0] == 't':
		N_confs += 1
print "I've got",N_confs,"configurations."

# keep track of the lines that start the configurations
t_lines = [0]*(N_confs)
N_confs=0
for i,l in enumerate(traj_by_lines):
	if l[0] == 't':
		t_lines[N_confs] = i 
		N_confs += 1
# we add one more line so that we know how to handle the last configuration.
t_lines += [len(traj_by_lines)]
# print them on a file.
culled_traj = []
new_N_confs = 0
for i in range(0,len(t_lines)-1):
	if i % N == 0:
		culled_traj += traj_by_lines[t_lines[i]:t_lines[i+1]]
		new_N_confs += 1
		
#print culled_traj

try:
	f = open(traj_filename,'w')
except:
	print "(cullTraj.py) ERROR: couldn't overwrite trajectory file",traj_filename,"."
	sys.exit(-1)

f.writelines(culled_traj)
ff.write("The presence of this file in this directory likely means that the script\n"
+"cullTraj.py has been launched in this directory. If you want to launch it again\n"
+"you'll have to remove this file.")

print "After the cull the file",traj_filename,"has",new_N_confs,"trajectories."
