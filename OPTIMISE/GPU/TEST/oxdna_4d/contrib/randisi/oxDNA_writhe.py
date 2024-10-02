#!/usr/bin/env python
'''
Compute the twist and writhe for each configuration of an oxDNA trajectory.

Assumes the system only contains the duplex.
'''
import numpy as np
import base, readers
import matplotlib.pyplot as plt
import sys 


def get_linking_number_from_system(system):
	'''
	Given a system (as returned from a LorenzoReader) containing two strands in a duplex,
	return its linking number. If one of the strands is circular, the duplex will be treated
	as if it were circular.
	'''
	s = system._strands
	circular = s[0]._circular or s[1]._circular
	
	# extract the positions
	positions = [np.array([nuc.cm_pos for nuc in strand._nucleotides]) for strand in s] 
	# from these, extract the un-normalised tangent vectors
	tangents = []
	if circular:
		for strand in positions:
			ahead = np.vstack([strand[1:],strand[0]])
			tan = ahead - strand
			tangents += [tan]
	# now compute the Linking number
	# The folowing is a bit obfuscated, but at least 2 order of magnitudes faster than doing it otherwise
	# First you compute cross product matrix, where cross[i][j] is the cross product of the i-th tangent vector of
  # the first strand with the j-th tangent vector of the second
	cross_m = np.cross( tangents[0][:,None,:], tangents[1]) # incredibly faster than doing it in the for loop
	# then compute the separation matrix between every particle pair of the two strands
	sep_m = positions[0][:, None, :] - positions[1]
	# divide it by the cubic norm
	dist = np.linalg.norm(sep_m,axis = 2)
	dist [ dist == 0] = np.inf
	sep_m /= dist[:,:,None]**3
	# finally, return the sum of all these terms
	return (sep_m * cross_m).sum()/(4.*np.pi)
	# The code above has been proven to be equivalent to this
#	for i,x, tx in zip(range(len(positions[0])),positions[0], tangents[0]):
#		for j,y, ty in zip(range(len(positions[1])),positions[1], tangents[1]):
#			cross = np.dot(tx,ty)#surprisingly much faster than np.cross(tx,ty)
#			sep = x - y
#			sep /= np.linalg.norm(sep)**3
#			lk += np.dot( cross, sep )
#	lk /= (4 * np.pi)
#
#	
#	return lk
def get_twist_from_system(system, first_bp = 0, last_bp = -1):
	'''
	Given a duplex configuration stored as a system returned by a LorenzoReader, return the twist.
	
	The twist is computed as described on pag.26 of Christian Matek's thesis. For each base-pair
	step in the system:
	1) We compute the vectors parallel to each base-pair in the base-pair step.
	2) The vectors are projected on a plane perpendicular to the line connecting the base-pair centers.
	3) The local twist angle tau_ij = arcos(r_i * r_i+1) is computed
	Then the local twist angles are summed together.
	'''
	pass

def get_bp_pos_tan_vec_from_system(s, first = 0, last_i = []):
	'''
	Given a duplex configuration stored as a system returned by a LorenzoReader, 
	return an array containing the base-pair position, tangent, and vector.

	'''
	circular = s._strands[0]._circular or s._strands[1]._circular 
	
	# get number of base-pairs and arguments
	n_bp_on_strand_A = len(s._strands[0]._nucleotides)
	n_bp_on_strand_B = len(s._strands[1]._nucleotides)
	if last_i == []: last = len(s._strands[0]._nucleotides) - 1
	else: 
		circular = False
		last = int(last_i)

	# perform some basic input sanity check
	if n_bp_on_strand_A != n_bp_on_strand_B:
		print 'ERROR (',__name__,'): function can only be used if the system is a duplex without sticky ends.'
		sys.exit()
	n_bp = n_bp_on_strand_A

	if first < 0 or last > n_bp - 1 or first > last:
		print 'ERROR (',__name__,'): it must be 0 <= first <= last <= n_bp - 1. That was not the case.'
		sys.exit()

	# compute the position of the base-pairs as the average position of the bases and return it
	pos = np.zeros([last - first + 1, 3])
	tan = np.zeros([last - first + 1, 3])
	vec = np.zeros([last - first + 1, 3])
	for array_line,i in enumerate(range(first, last+1)):
		i_c = n_bp - 1 - i
		base_1_pos = s._strands[0]._nucleotides[i].get_pos_base()
		base_2_pos = s._strands[1]._nucleotides[i_c].get_pos_base()
		pos[array_line] = (base_1_pos + base_2_pos) / 2.0
		vec[array_line] = (base_1_pos - base_2_pos) / 2.0
		vec[array_line] /= np.linalg.norm(vec[array_line])
		if i != first:
			tangent = pos[array_line] - pos[array_line - 1]
			tangent -= s._box * np.rint(tangent / s._box)
			tangent /= np.linalg.norm(tangent)
			tan[array_line - 1] = tangent
	if circular:
		tangent = pos[0] - pos[-1]
		tangent -= s._box * np.rint(tangent / s._box)
		tangent /= np.linalg.norm(tangent)
		tan[ - 1] = tangent
		return pos, tan, vec
	else:
		return pos, tan[:-1], vec

def get_bp_pos(s, first = 0, last_i = []):
	'''
	Given a duplex as a system from a LorenzoReader, extract the position of the base-pairs in some interval.
	
	If no interval is given, the base-pairs from first to last will be extracted.
	'''
	# get number of base-pairs and arguments
	n_bp_on_strand_A = len(s._strands[0]._nucleotides)
	n_bp_on_strand_B = len(s._strands[1]._nucleotides)
	if last_i == []: last = len(s._strands[0]._nucleotides) - 1
	else: last = int(last_i)

	# perform some basic input sanity check
	if n_bp_on_strand_A != n_bp_on_strand_B:
		print 'ERROR (',__name__,'): function can only be used if the system is a duplex without sticky ends.'
		sys.exit()
	n_bp = n_bp_on_strand_A

	if first < 0 or last > n_bp - 1 or first > last:
		print 'ERROR (',__name__,'): it must be 0 <= first <= last <= n_bp - 1. That was not the case.'
		print '                      first =',first,', last =',last,', n_bp =',n_bp
		sys.exit()

	# compute the position of the base-pairs as the average position of the bases and return it
	bp_array = np.zeros([last - first + 1,3])
	for array_line,i in enumerate(range(first, last+1)):
		i_c = n_bp - 1 - i
		base_1 = s._strands[0]._nucleotides[i]
		base_2 = s._strands[1]._nucleotides[i_c]
		bp_array[array_line] = (base_1.get_pos_base() + base_2.get_pos_base()) / 2.0
	return bp_array
	
def get_tangent_vectors_from_bp_pos(bp_pos, box, norm = False):
	'''
	Given an array containing the position of some base-pairs, return an array containing the tangent vectors.
	
	If norm equals True, the vectors returned will already be normalised. Defaults to False
	'''
	n_bp = bp_pos.shape[0]
	tangents_array = np.zeros([n_bp-1, 3])
	for i in xrange(n_bp-1):
		tangents_array[i] = (bp_pos[i+1] - bp_pos[i]) 
		tangents_array[i] -= box * np.rint(tangents_array[i] / box)
		if norm: tangents_array[i] /= np.linalg.norm(tangents_array[i])
	return tangents_array
	
def get_twist_writhe_from_pos_tan_vec(pos,tan,vec, is_circular = False):
	'''
	Given the base-pair position vector, the base-pair orientation vectors, and the tangent vectors of a duplex,
	compute the twist and writhe.
	
	This assumes a duplex, but you can treat a circular configuration assuming that the first and last nucleotides are linked (i.e., defining tan[-1] and extending the for loop in the twist).
	'''
	# compute twist
	twist = 0
	for i in range(len(tan)):
		i_p = i + 1
		if i_p == len(vec): i_p = 0
		# project the base-pair vectors on the plane perpendicular to their tangent vector
		if abs(np.linalg.norm(tan[i]) -1) > 1e-4:
			print 'error tan[',i,'] = ',tan[i]
		if abs(np.linalg.norm(vec[i]) -1) > 1e-4:
			print 'error vec[',i,'] = ',vec[i]
		if abs(np.linalg.norm(vec[i_p]) -1) > 1e-4:
			print 'error vec[',i_p,'] = ',vec[i_p]
		rp = vec[i]    - vec[i]   * np.dot(vec[i],   tan[i])
		rp /= np.linalg.norm(rp)
		rpp = vec[i_p] - vec[i_p] * np.dot(vec[i_p], tan[i])
		rpp /= np.linalg.norm(rpp)
		# add the number of turns separating them.
		t = np.arccos(np.dot(rp, rpp)) / (2 * np.pi)
		if np.isnan(t):
			print rp, rpp, np.dot(rp,rpp)
		twist += t

	# fast writhe computation
	# get the un-normalised tangents
	tan_unnorm = []
	if is_circular:
		ahead = np.vstack([pos[1:],pos[0]])
		tan_unnorm = ahead - pos 
		# I tried using this one too, but it doesn't help much
		#pos_avg = (pos + ahead)*0.5
	else:
		print 'expand the code so that it works on non-circular strands'
		print 'just make sure that the last tangent is undefined, plus make sure not to use the last position vector'
		sys.exit()

	cross_m = np.cross( tan_unnorm[:,None,:], tan_unnorm) # incredibly faster than doing it in the for loop
	# then compute the separation matrix between every particle pair of the two strands
	sep_m = pos[:, None, :] - pos
	# divide it by the cubic norm
	dist = np.linalg.norm(sep_m,axis = 2)
	dist [ dist == 0] = np.inf
	sep_m /= dist[:,:,None]**3
	# finally, return the sum of all these terms
	writhe = (sep_m * cross_m).sum()/(4.*np.pi)
	return twist, writhe
###### OLD version
def get_twist_slow_writhe_from_pos_tan_vec(pos,tan,vec, is_circular = False):
	'''
	Given the base-pair position vector, the base-pair orientation vectors, and the tangent vectors of a duplex,
	compute the twist and writhe.
	
	This assumes a duplex, but you can treat a circular configuration assuming that the first and last nucleotides are linked (i.e., defining tan[-1] and extending the for loop in the twist).
	'''
	# compute twist
	twist = 0
	for i in range(len(tan)):
		i_p = i + 1
		if i_p == len(vec): i_p = 0
		# project the base-pair vectors on the plane perpendicular to their tangent vector
		if abs(np.linalg.norm(tan[i]) -1) > 1e-4:
			print 'error tan[',i,'] = ',tan[i]
		if abs(np.linalg.norm(vec[i]) -1) > 1e-4:
			print 'error vec[',i,'] = ',vec[i]
		if abs(np.linalg.norm(vec[i_p]) -1) > 1e-4:
			print 'error vec[',i_p,'] = ',vec[i_p]
		rp = vec[i]    - vec[i]   * np.dot(vec[i],   tan[i])
		rp /= np.linalg.norm(rp)
		rpp = vec[i_p] - vec[i_p] * np.dot(vec[i_p], tan[i])
		rpp /= np.linalg.norm(rpp)
		# add the number of turns separating them.
		t = np.arccos(np.dot(rp, rpp)) / (2 * np.pi)
		if np.isnan(t):
			print rp, rpp, np.dot(rp,rpp)
		twist += t
	# compute writhe
	writhe = 0
	for i in range(len(tan)):
		for j in range(i):
			cross = np.cross( tan[j], tan[i])
			sep = pos[j] - pos[i]
			sep /= np.linalg.norm(sep)**3
			writhe += np.dot(cross,sep) / (2 * np.pi)
	return twist, writhe
######
		
	

def get_main_output(traj_file, out_filename, first_bp = 0, last_bp = -1):
	'''
	Take the trajectory file_name and the first and last base-pair to consider, then compute writhe and twist.
	'''
	# process input and do a basic sanity check
	l = readers.LorenzoReader(traj_file, 'prova.top')
	s = l.get_system()
	while s:
		base_pairs_pos = get_bp_pos(s, first_bp, last_bp)
		tangents = get_tangent_vectors_from_bp_pos(base_pairs_pos, s._box, True)
		#writhe = get_writhe_from_system(t_array)
		#twist = get_twist_from_system(t_array)
		circular = s._strands[0]._circular or s._strands[1]._circular
		pos, tan, vec = get_bp_pos_tan_vec_from_system(s, first_bp, last_bp)
		for b,p,t,tt in zip(base_pairs_pos[:-1], pos,tangents[:-1],tan):
			if np.linalg.norm(b-p) > 1e-3:
				print 'errorbp:',b,p
				sys.exit()
			if np.linalg.norm(t-tt) > 1e-3:
				print 'errorttt:',t,tt
				sys.exit()
		# compute twist and writhe
		tw, wr = get_twist_writhe_from_pos_tan_vec(pos,tan,vec,circular)
		print '%e\t%.4lf\t\t%.4lf\t\t%.4lf\t\t%.4lf'%(s._time, tw, wr, tw + wr, get_linking_number_from_system(s))
		s = l.get_system()

if __name__ == '__main__':
	first_bp = 0
	last_bp = -1
	out_filename = 'twist_writhe.dat'
	if len(sys.argv) < 2:
		print 'USAGE:',sys.argv[0],' <traj_file> [first_bp] [last_bp]'
		print '       reads a duplex trajectory, then writes a file',out_filename,
		print 'with twist and writhe for each configuration'
		print '       <traj_file> is the trajectory filename (mandatory)'
		print '       [first_bp] is the first base-pair to consider (defaults to',first_bp,')'
		print '       [last_bp] is the last base-pair to consider (defaults to',last_bp,')'
		sys.exit()
	traj_file = sys.argv[1]
	if len(sys.argv) > 2:
		first_bp = int(sys.argv[2])
	if len(sys.argv) > 3:
		last_bp = int(sys.argv[3])
	if last_bp == -1: last_bp = []
	print '#time		twist		writhe		twist+writhe		Lk'
	output = get_main_output(traj_file, out_filename, first_bp, last_bp)

