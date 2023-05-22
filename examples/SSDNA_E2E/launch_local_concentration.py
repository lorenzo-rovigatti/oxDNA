#! /usr/bin/env python
'''
Script that launches simulations of ssDNA and measures the e2e distance, as well as the distance between nucleotide 1 and nucleotide N-2, between nucleotide 2 and nucleotide N-3, etc.
'''
import sys, os, re
import math
import numpy as np
import shutil
import stat
import subprocess
import tempfile


def is_valid_DNA_sequence(sequence):
	'''
	Returns true if the string only contains A,T,C,G, and false otherwise.

	Uses a magic code snippet taken from here
	http://stackoverflow.com/questions/1323364/in-python-how-to-check-if-a-string-only-contains-certain-characters
	'''
	search = re.compile(r'[^ATCG.]').search
	return not bool(search(sequence))


def assert_is_dir(path):
	''' Assert that <path> is a valid directory. Notify the user if the assertion fails.'''
	if not os.path.isdir(path):
		if os.path.exists(path):
			print("ERROR:", path, "exists but is not a directory.")
			sys.exit(-1)
		else:
			print("ERROR:", path, "doesn't exist.")
			sys.exit(-1)


def replace(file_path, pattern, subst):
    # Create temp file
    fh, abs_path = tempfile.mkstemp()
    new_file = open(abs_path, 'w')
    old_file = open(file_path)
    for line in old_file:
        new_file.write(re.sub(pattern, subst, line))
    # close temp file
    new_file.close()
    os.close(fh)
    old_file.close()
    # Remove original file
    os.remove(file_path)
    # Move new file
    shutil.move(abs_path, file_path)


possible_seq_id_description = \
"---polyX_N, where X is any letter among 'ATGC',\n   and N is the number of nucleotides.\n   This signifies a strand with N identical nucleotides.\n   Example: polyA_33\n\n" + \
"---N, where N is the number of nucleotides.\n   This signifies a strand with N identical nucleotides\n" + "   of the average oxDNA model.\n   Example: 7\n\n" + \
"---direction-sequence, where direction is either '5p-3p'\n" + "   or 'oxdna', depending if the sequence is to be read \n" + "   from the 5\' end to the 3\' end (as is custumary in\n" + "   biology and nanotechnology) or from the 3\' end to the\n" + "   5\', as it's done in oxDNA.\n   Example: 5p3p-CACATA\n\n"


def process_seq_id(seq_id):
	'''
	Takes a sequences identifier and spits out a sequence and the value of the use_average_seq argument to be used.
	'''
	# if it's an integer, just run a single strand of some length
	try:
		N = int(seq_id)
		return N * 'A', True
	except:
		pass
	# handle arguments like polyA_30, polyG44, polyT_50, etc.
	if seq_id[:4] == 'poly':
		if seq_id[4] in 'ATGC':
			base = seq_id[4]
			start_n = 5
			if seq_id[5] in '-_':
				start_n += 1
			try:
				N = int(seq_id[start_n:])
				return N * base, False
			except Exception as e:
				print(e)
				print(seq_id[start_n:])
				print('ERROR (process_seq_id): sequence identifier starting with', seq_id[:5], 'must end with a number.')
				sys.exit(-1)
		
		else:
			print('ERROR (process_seq_id): base must be any of A,T,G,C, but was', seq_id[4])
			sys.exit(-1)
	# handle textual sequences with the 5p3p- or oxdna- prefix (e.g. ATATAGACGA)
	elif seq_id[:5] == '5p3p-' or seq_id[:6].lower() == 'oxdna-':
		if seq_id == '5p3p-' or seq_id == 'oxdna-':
			print('ERROR: sequence must not be empty.')
			sys.exit(-1)
		if seq_id[:5] == '5p3p-':
			sequence = seq_id[5:].upper()
			if not is_valid_DNA_sequence(sequence):
				print('ERROR: sequence', seq_id[5:], 'is not a valid sequence since it contains characters other than ATCGatcg.')
				sys.exit(-1)
		elif seq_id[:6].lower() == 'oxdna-':
			sequence = seq_id[6:].upper()
			if not is_valid_DNA_sequence(sequence):
				print('ERROR: sequence', seq_id[6:], 'is not a valid sequence since it contains characters other than ATCGatcg.')
				sys.exit(-1)
		else:
			print('CODING ERROR: this condition should never be reached!')
			print('              if you see this there\'s an error in the code, please contact the programmer.')
			sys.exit(-1)
		return sequence, False
	else:
		print('ERROR: (process_seq_id): couldn\'t process sequence identifier', seq_id)
		if is_valid_DNA_sequence(seq_id.upper()):
			print('WARNING: if the sequence identifier is a DNA sequence, it must start\nwith either 5p3p- or oxdna- depending on whether the sequence is given\nfrom the 5\' end to the 3\' end (which is the standard in biology and\nnanotechnology) or from the 3\' end to the 5\' end (which is oxDNA\'s convention).\nTherefore e.g. 5p3p-' + seq_id + ' and oxdna-' + seq_id + '\nare both valid sequence identifiers, though they mean different things.')
		sys.exit(-1)


def all_bonds_op_file(N, op_file_path):
	'''
	Writes an order parameter file where the order parameter is given by the number of all the hydrogen bonds in the configuration.
	
	N is the number of nucleotides, and op_file_path is the path of the order parameter.
	'''

	with open(op_file_path, 'w') as f:
		f.write('{\norder_parameter = bond\nname = every_single_bond\n')
		pair_index = 0
		for i in range(0, N):
			for j in range(i + 1, N):
				f.write('pair%d = %d, %d\n' % (pair_index, i, j))
				pair_index += 1
		

# Some scripting parameters
oxDNA_path = '/users/randisif/tep/CURRENToxdna-code/oxDNA/build/bin/oxDNA'
oxDNA_extra_arguments = ' '
src_input_dir = 'src_input'
input_file_name = 'input'
sequence_file_name = 'seq.txt'

n_distances_str = str(8)
# handle the arguments
if len(sys.argv) < 2:
	print('USAGE:', sys.argv[0], '<sequence_identifier> [n_distances_to_measure]')
	print('       where:')
	print('       <sequence_identifier> can take several values, as described below.')
	print('       [n_distances_to_measure] (optional) specifies how many nucleotides')
	print('       at the end of each chain have to be monitored (default is', n_distances_str, ').\n')
	print('OUTCOME: This program creates a directory with name <sequence_identifier>')
	print('         in which everything is ready to simulate a ssDNA with oxDNA and')
	print('         measure the distances between the nucleotides at the respective')
	print('         ends of the chain, e.g. betweer the first and the last nucleotides,')
	print('         between the second and the first-but-last, etc., for a total of')
	print('         [n_distances_to_measure] pairs.\n')
	print('POSSIBLE SEQUENCE IDENTIFIERS:\n\n', possible_seq_id_description)
	sys.exit(-1)
if len(sys.argv) > 2:
	n_distances_str = sys.argv[2]

seq_id = sys.argv[1]
n_distances = int(n_distances_str)

# figure out the sequence from the seq_id
sequence, use_average_seq = process_seq_id(seq_id)
N = len(sequence)

# make sure the source directory exists
assert_is_dir(src_input_dir)
# make sure the sequence id directory does NOT exist
if os.path.exists(seq_id):
	print('ERROR: directory', seq_id, 'already exists. Delete it if you want and then call me againg.')
	sys.exit(-1)
# create the simulation directory and copy the input file there
os.mkdir(seq_id)
shutil.copyfile(src_input_dir + '/' + input_file_name, seq_id + '/' + input_file_name)

os.chdir(seq_id)
# create the sequence file
with open("prova.top", "w") as f:
	print(len(sequence), 1, file=f)
	for i, c in enumerate(sequence):
		n3 = i - 1
		n5 = i + 1
		if n5 == len(sequence):
			n5 = -1
		print(1, c, n3, n5, file=f)
		
with open(sequence_file_name, 'w') as f:
	f.write(sequence)
# generate the configuration and topology file
# generate.read_strands(sequence_file_name,len(sequence)*0.3+10)
# edit the input file
# set the value of use_average_seq
replace(input_file_name, '^use_average_seq.*', 'use_average_seq =' + str(int(use_average_seq)))
# decide which distances to measure
columns = ''
for i in range(1, n_distances + 1):
	p_1 = i - 1
	p_2 = N - i
	if p_1 >= p_2: break
	columns += '	col_' + str(i) + ' = {\n' + \
						'		type = distance\n' + \
						'		particle_1 = ' + str(p_1) + '\n' + \
						'		particle_2 = ' + str(p_2) + '\n' + \
						'		PBC = false\n	}\n'
replace(input_file_name, 'DISTANCE_COLUMNS', columns)
# when the average model is being used, measure stacking as well.
if use_average_seq:
	with open(input_file_name, 'a') as f:
		unstacked_obs = 'data_output_2 = {\n' + \
										'	name = unstacked_profile.dat\n' + \
										'	print_every = 1e3\n' + \
										'	col_1 = {\n' + \
										'		type = unstacked_list\n' + \
										'	}\n+}\n'
		print(unstacked_obs, file=f)
# print the order parameter file		
all_bonds_op_file(N, 'op.txt')
# print the weight file
with open('wfile.txt', 'w') as f:
	f.write('0 1\n')

print("########## ALL DONE - EVERYTHING WENT OK")
print("To change directory, generate the initial configuration and run your simulation, run:")
print("cd", seq_id)
print("confGenerator input", len(sequence))
print("oxDNA input")
print("where oxDNA is the path to your oxDNA executable.")

