#!/usr/bin/env python
import base
try:
    import numpy as np
except:
    import mynumpy as np
import math
import sys
import utils

def my_norm(vec):
	return math.sqrt(np.dot(vec,vec))

# need to add this to avoid stupid warning about invalid division at line dir /= np.sqrt(np.dot(dir,dir))
try:
    np.seterr(invalid='ignore')
except:
    pass

BP = "bp"
DEGREES = "degrees"
DEG = "degrees"
RAD = "radiants"
RADIANTS = "radiants"
BASE_PAIRS = "bp"


base.CM_CENTER_DS -= 0.2

class StrandGenerator (object):
    def generate(self, bp, sequence=None, start_pos=np.array([0, 0, 0]), dir=np.array([0, 0, 1]), perp=False, double=True, rot=None,base_base_distance = None):
	base.CM_CENTER_DS -= 0.
	if(base_base_distance == None):
		base_base_distance = base.BASE_BASE
	if(rot == None):
		rot = 35.9
        # we need numpy array for these
        start_pos = np.array(start_pos, dtype=float)
        dir = np.array(dir, dtype=float)
        if sequence == None:
            sequence = np.random.randint(0, 4, bp)
        elif len(sequence) != bp:
            n = bp - len(sequence)
            sequence += np.random.randint(0, 4, n)
            base.Logger.log("sequence is too short, adding %d random bases" % n, Logger.WARNING)
            
        # create the sequence of the second strand as made of complementary bases
        sequence2 = [3-s for s in sequence]
        sequence2.reverse()
        
        # we need to find a vector orthogonal to dir
        dir_norm = np.sqrt(np.dot(dir,dir))
        if dir_norm < 1e-10:
            base.Logger.log("direction must be a valid vector, defaulting to (0, 0, 1)", base.Logger.WARNING)
            dir = np.array([0, 0, 1])
        else: dir /= dir_norm

        if perp is None or perp is False:
            v1 = np.random.random_sample(3)
            v1 -= dir * (np.dot(dir, v1))
            v1 /= np.sqrt(sum(v1*v1))
        else:
            v1 = perp;
        
        # and we need to generate a rotational matrix
        R = utils.get_rotation_matrix(dir, np.deg2rad(rot))
        #R = get_rotation_matrix(dir, np.deg2rad(35.9))
        #R = utils.get_rotation_matrix(dir, [1, BP])
        
        ns1 = base.Strand()

        a1 = v1
        #a1 = np.dot (R0, a1)
        rb = np.array(start_pos)
        a3 = dir
        for i in range(bp):
            ns1.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence[i]))
            if i != bp-1:
                a1 = np.dot(R, a1)
                rb += a3 * base_base_distance
            
        if double == True:
            a1 = -a1
            a3 = -dir
            R = R.transpose()
            ns2 = base.Strand()
            for i in range(bp):
                ns2.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence2[i]))
                a1 = np.dot(R, a1)
                rb += a3 * base_base_distance 
                
            return ns1, ns2
        else: return ns1



    def a_generate(self, bp, sequence=None, start_pos=np.array([0, 0, 0]), dir=np.array([0, 0, 1]), perp=False, double=True, rot=None,base_base_distance = None,inclination=None,diameter=2.35):
	if inclination == None:
		inclination = 15.5

	bp_backback_distance = 2.0
	cord = math.cos(np.deg2rad(inclination)) * bp_backback_distance
	center_to_cord = math.sqrt((diameter/2.)**2 - (cord/2.)**2)
	
	if(base_base_distance == None):
		base_base_distance = 0.3287  #a-rna, source: neidle book
	if(rot == None):
		rot = 32.7 #a-dna, source: wiki

	# we need numpy array for these
	start_pos = np.array(start_pos, dtype=float)
	dir = np.array(dir, dtype=float)
	if sequence == None:
	    sequence = np.random.randint(0, 4, bp)
	elif len(sequence) != bp:
	    n = bp - len(sequence)
	    sequence += np.random.randint(0, 4, n)
	    base.Logger.log("sequence is too short, adding %d random bases" % n, Logger.WARNING)
	    
	# create the sequence of the second strand as made of complementary bases
	sequence2 = [3-s for s in sequence]
	sequence2.reverse()

	# we need to find a vector orthogonal to dir
	dir_norm = np.sqrt(np.dot(dir,dir))
	if dir_norm < 1e-10:
	    base.Logger.log("direction must be a valid vector, defaulting to (0, 0, 1)", base.Logger.WARNING)
	    dir = np.array([0, 0, 1])
	else: dir /= dir_norm

	#x1, y1, z1 = center_to_cord,  cord / 2., -(bp_backback_distance / 2.) * np.sin (np.deg2rad(inclination))
	#x2, y2, z2 = center_to_cord, -cord / 2., +(bp_backback_distance / 2.) * np.sin (np.deg2rad(inclination))
	
	x2, y2, z2 = center_to_cord,  -cord / 2., +(bp_backback_distance / 2.) * np.sin (np.deg2rad(inclination))
	x1, y1, z1 = center_to_cord, +cord / 2., -(bp_backback_distance / 2.) * np.sin (np.deg2rad(inclination))

	r1 = np.array ([x1, y1, z1])
	r2 = np.array ([x2, y2, z2])
	r1_to_r2 = r2 - r1
	
	ZAXIS = np.array ([0., 0., 1.])

	RISE = base_base_distance

	R = utils.get_rotation_matrix (ZAXIS, np.deg2rad(rot))
	
	ns1 = base.Strand()
	
	for i in xrange (len(sequence)):
		r1 = np.dot (R, r1) + RISE * ZAXIS
		r2 = np.dot (R, r2) + RISE * ZAXIS
		r1_to_r2 = r2 - r1	
		a1 = r1_to_r2 / math.sqrt(np.dot(r1_to_r2,r1_to_r2))
		a1proj = np.array([a1[0],a1[1],0.])
		a1projnorm = math.sqrt(np.dot(a1proj, a1proj))
		a3 = -(-math.cos(np.deg2rad(inclination)) * ZAXIS + math.sin(np.deg2rad(inclination))* a1proj / a1projnorm)
		#a3 = math.cos(np.deg2rad(inclination)) * ZAXIS - math.sin(np.deg2rad(inclination))* a1proj / a1projnorm
		ns1.add_nucleotide(base.Nucleotide(r1 + base.CM_CENTER_DS * a1, a1, a3, sequence[i]))
        
	if double == True:
		a1 = -a1
		a3 = -dir
		R = R.transpose()
		ns2 = base.Strand()
		for i in xrange (len(sequence)):
			r1_to_r2 = r2 - r1	
			a1 = -r1_to_r2 / math.sqrt(np.dot(r1_to_r2,r1_to_r2))
			a1proj = np.array([a1[0],a1[1],0.])
			a1projnorm = math.sqrt(np.dot(a1proj, a1proj))
			a3 = -(math.cos(np.deg2rad(inclination)) * ZAXIS + math.sin(np.deg2rad(inclination))* a1proj / a1projnorm)
			#a3 = -math.cos(np.deg2rad(inclination)) * ZAXIS - math.sin(np.deg2rad(inclination))* a1proj / a1projnorm
			ns2.add_nucleotide(base.Nucleotide(r2 + base.CM_CENTER_DS * a1, a1, a3, sequence2[i]))
			r1 = np.dot (R, r1) - RISE * ZAXIS
			r2 = np.dot (R, r2) - RISE * ZAXIS

		#display_info(ns1,ns2,np.deg2rad(rot),RISE,diameter/2.)
		return ns1, ns2
		
	else:
		return ns1

   
import random
import numpy as np
 
def do_strands(s,sequence,double_type=True,rotation=None,rise_per_bp=None):
    g = StrandGenerator()
    added = False
    #s = base.System([box, box, box])
    if double_type == False:
	oN = s._N
	#iprint oN
	strand =  g.a_generate(len(sequence), sequence,double=double_type,rot=rotation,base_base_distance=rise_per_bp)

	while not added:
        	pos = 10. * np.array([random.random(),random.random(),random.random()])
		#print pos
		strand.translate(np.array([0.2*s._box[0]*random.random(),0.2*s._box[1]*random.random(), 0.2*s._box[2]*random.random()]))
		#strand.translate(pos)
        	added = s.add_strands(strand,check_overlap=True)	
		if not added:
			print >> sys.stderr, 'Overlap detected, tryinhg to place the strand elsewhere'
	if oN == s._N:
		print 'Adding strand ', sequence, 'failed', s._N
    else:
        s1,s2 = g.a_generate(len(sequence), sequence,double=double_type,rot=rotation,base_base_distance=rise_per_bp)
	#pos =  np.array(s1._nucleotides[0].cm_pos) - np.array(s2._nucleotides[1].cm_pos)
	#print 'Strands for duplex generated, distance between COM is ', np.sqrt(np.dot(pos,pos))
	oN = s._N
	while not added:
		pos = np.array([0.2*s._box[0]*random.random(),0.2*s._box[1]*random.random(), 0.2*s._box[2]*random.random()])
        	#pos = 10. * np.array([random.random(),random.random(),random.random()])
		#print pos
		s1.translate(pos)
		s2.translate(pos)
		#strand.translate(pos)
        	added = s.add_strands([s1,s2],check_overlap=True)	
		if not added:
			print >> sys.stderr, 'Overlap detected, tryinhg to place the strands elsewhere'
	if oN == s._N:
		print 'Adding strand ', sequence, 'failed', s._N
    #s.print_lorenzo_output (prefix+".conf", prefix+".top")
    #s.print_pdb_output_chimera(prefix+".conf.pdb")
   
		
def process_sequence_file(fname,box):
	s = base.System([box,box,box])
	fin = open(fname,'r')
	for line in fin.readlines():
		vals = line.split()
		if(vals[0] == 'DOUBLE' ):
			seq = vals[1].strip()	
			track = [base.base_to_number[x] for x in seq]
			print >> sys.stderr, 'Adding duplex 3\'-', seq,'-5\''
			do_strands(s,track,True)	
		else:
			seq = vals[0].strip()
			track = [base.base_to_number[x] for x in seq]
			do_strands(s,track,False)
      			print >> sys.stderr, 'Adding sequence 3\'-', seq,'-5\''
	s.print_lorenzo_output (prefix+".conf", prefix+".top")
	return s

	
if (len(sys.argv) < 3):
	print "Usage: %s <files_with_sequence> prefix [box = 20]" % (sys.argv[0])
	sys.exit(1)

if len(sys.argv) >= 4:
	box = float(sys.argv[3])
else: 
	box = 20
file_seq = sys.argv[1]
prefix = sys.argv[2]


random.seed(4)
process_sequence_file(file_seq,box)

