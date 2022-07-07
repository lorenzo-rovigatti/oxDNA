#!/usr/bin/env python

'''
Forward Flux sampling: shoot from one interface to the next

Flavio
'''

####################################
# edit here
####################################
# number of successes, can be changed via command line
desired_success_count = 100

# executable set-up
executable = '../oxDNA'
input = 'input'
logfilename = 'ffs.log'
starting_conf_pattern = '../I0I1/success*'
success_pattern = './success_'
keep_undetermined = True            # if True, the program saves undetermined confs
undetermined_pattern = './undefin_'

prune_prob = 0.75
#prune_prob = 0.0

# interfaces 
# interface lambda_{-1}
lambda_f_name = 'dist'
lambda_f_value = 4
lambda_f_compar = '>='

# interface lambda_{prune}
lambda_pf_name = 'dist'
lambda_pf_value = 1
lambda_pf_compar = '>='

# interface lambda_{n}
lambda_n_name = 'native'
lambda_n_value = 3

# interface lambda_{n+1}
lambda_m_name = 'native'
lambda_m_value = 8
lambda_m_compar = '>='
####################################

import os, sys, getopt
import subprocess as sp
import time, random as rnd, tempfile as tf
import shutil, glob
from multiprocessing import Process, Lock, JoinableQueue, Value, Array

def usage():
	print >> sys.stderr, 'usage: %s %s' % (sys.argv[0], '[-n <num_sucesses>] [-s <seed>] [-c <ncpus>] [-k <success_count>]') 

try:
	opts, files = getopt.gnu_getopt(sys.argv[1:], "n:s:c:k:v")
except getopt.GetoptError as err:
	usage()
	sys.exit(2)

initial_seed = time.time()
Verbose = False
ncpus = 1
initial_success_count = 0

try:
	for o in opts:
		k, v = o
		if k == '-n':
			desired_success_count = int(v)
		elif k == '-s':
			initial_seed = int(v)
		elif k == '-c':
			ncpus = int(v)
		elif k == '-v':
			Verbose = True
		elif k == '-k':
			initial_success_count = int(v)
		else:
			print >> sys.stderr, "Warning: option %s not recognized" % (k)
except:
	print >> sys.stderr, "Error parsing options"
	sys.exit(3)

log_lock = Lock()
log_file = open (logfilename, 'w', 1)
def log(text):
	log_lock.acquire()
	log_file.write(text + '\n')
	if Verbose:
		print >> sys.stdout, text
	log_lock.release()

# starting configurations
starting_confs = glob.glob(starting_conf_pattern)
log ("Main: Found %d configurations with pattern: %s" % (len(starting_confs), starting_conf_pattern))
if len(starting_confs) < 1:
	print >> sys.stderr, "0 starting configurations! aborting"
	sys.exit(2)

# check that we can write to the success pattern
try:
	checkfile = open (success_pattern + '0', 'w')
	checkfile.close()
	os.remove (success_pattern + '0')
except:
	print >> sys.stderr, "could not write to success_pattern", success_pattern
	sys.exit(3)
	
success_lock = Lock()
success_count = Value ('i', initial_success_count)
attempt_count = Value ('i', 0)
success_from = Array ('i', len(starting_confs)) # zeroed by default
success_pruning_count = Value ('i', 0)
attempt_from = Array ('i', len(starting_confs)) # zeroed by default
undetermined_lock = Lock()
undetermined_count = Value ('i', 0)

# timings
time_succeeding = Value ('f', 0)

# write the condition file
if os.path.exists('conditions.txt'):
	log ("Main: Warning: overwriting conditions file")
condition_file = open('conditions.txt', "w")
condition_file.write("action = stop_or\n")
condition_file.write("condition1 = {\n%s %s %s\n}\n" % (lambda_pf_name, lambda_pf_compar, str(lambda_pf_value)))
condition_file.write("condition2 = {\n%s %s %s\n}\n" % (lambda_m_name, lambda_m_compar, str(lambda_m_value)))
condition_file.close()

condition_file = open('conditions2.txt', "w")
condition_file.write("action = stop_or\n")
condition_file.write("condition1 = {\n%s %s %s\n}\n" % (lambda_f_name, lambda_f_compar, str(lambda_f_value)))
condition_file.write("condition2 = {\n%s %s %s\n}\n" % (lambda_m_name, lambda_m_compar, str(lambda_m_value)))
condition_file.close()

# base command line; all features that need to be in the input file
# must be specified here
base_command = [executable, input, 'print_energy_every=1e5', 'print_conf_every=1e6','no_stdout_energy=1','refresh_velocity=0','restart_step_counter=1']
base_command_string = ''.join (str(w) + ' ' for w in base_command)
log("Main: COMMAND: " + base_command_string)

if not os.path.exists(input):
	log ("the input file provided (%s) does not exist. Aborting" % (input))
	sys.exit(-3)

if not os.path.exists(executable):
	log ("the executable file provided (%s) does not exist. Aborting" % (executable))
	sys.exit(-3)

rnd.seed (initial_seed)

log ("Main: STARTING new shooting for %d" % desired_success_count)
log ("Main: number of initial confs: %d" % len(starting_confs))
desired_success_count += initial_success_count

# this function does the work of running the simulation, identifying a
# success or a failure, and taking appropriate actions
def f(idx):
	log ("Worker %d started" % idx)
	global success_count, success_from, success_pruning_count
	while success_count.value < desired_success_count:
		# choose a starting configuration
		#log ("Worker %d started" % idx)
		conf_index = rnd.choice (range(1, len(starting_confs))) 
		conf_file = starting_confs[conf_index]
		global attempt_from, attempt_count, time_succeeding 
		attempt_count.value += 1
		attempt_from[conf_index] += 1
		
		mytime = time.time()	
		
		# the seed is the index + initial seed, and the last_conf has an index as well
		seed = initial_seed + attempt_count.value
		last_conf = 'last_conf' + str(idx)
		
		# edit the command to be launched
		command = base_command + ['conf_file=%s' % (conf_file), 'lastconf_file=%s' % (last_conf), 'seed=%d' % (seed), 'ffs_file=conditions.txt']
		
		# open a file to handle the output
		output = tf.TemporaryFile ('r+')
		
		# here we run the command
		# print command
		r = sp.call (command, stdout=output, stderr=sp.STDOUT)
		if r != 0:
			print >> sys.stderr, "Error running program"
			print >> sys.stderr, "command line:"
			txt = ''
			for c in command:
				txt += c + ' '
			print >> sys.stderr, txt
			print >> sys.stderr, 'output:'
			output.seek(0)
			for l in output.readlines():
				print >> sys.stderr, l,
			output.close()
			sys.exit(-2)
		
		# now we process the output to find out wether the run was a success
		# (interface lambda_m reached) or a failure (interface lamda_f reached)
		output.seek(0)
		for line in output.readlines():
			words = line.split()
			if len(words) > 1:
				if words[1] == 'FFS' and words[2] == 'final':
					data = [w for w in words[4:]]
		op_names = data[::2]
		op_value = data[1::2]
		op_values = {}
		for ii, name in enumerate(op_names):
			op_values[name[:-1]] = float(op_value[ii][:-1])
		
		# now op_values is a dictionary representing the status of the final
		# configuration. We now need to find out whether it is a success or
		# a failure.
		# print op_values, 'op_values["%s"] %s %s' % (lambda_m_name, lambda_m_compar, str(lambda_m_value)), 'op_values["%s"] %s %s' % (lambda_f_name, lambda_f_compar, str(lambda_f_value))
		success = eval ('op_values["%s"] %s %s' % (lambda_m_name, lambda_m_compar, str(lambda_m_value)))
		#failure = eval ('op_values["%s"] %s %s' % (lambda_f_name, lambda_f_compar, str(lambda_f_value)))
		failure = eval ('op_values["%s"] %s %s' % (lambda_pf_name, lambda_pf_compar, str(lambda_pf_value)))
		
		if success and not failure:
			with success_lock:
				success_count.value += 1
				success_from[conf_index] += 1
				shutil.copy (last_conf, success_pattern + str(success_count.value))
				time_succeeding.value += (time.time () - mytime)
			os.remove(last_conf)
			log ("SUCCESS: worker %d: starting from conf_index %d and seed %d" % (idx, conf_index, seed))
		elif not success and failure:
			# do else
			log ("FAILURE (AT PRUNING INTERFACE): worker %d: starting from conf_index %d and seed %d" % (idx, conf_index, seed))
			
			# PRUNING
			if rnd.random() > prune_prob:
				log ("PRUNING: worker %d: starting from conf_index %d and seed %d" % (idx, conf_index, seed))
				# edit the command to be launched
				command = base_command + ['conf_file=%s' % (last_conf), 'lastconf_file=%s' % (last_conf), 'seed=%d' % (seed), 'ffs_file=conditions2.txt']
				
				# open a file to handle the output
				output = tf.TemporaryFile ('r+')
				
				# here we run the command
				# print command
				r = sp.call (command, stdout=output, stderr=sp.STDOUT)
				if r != 0:
					print >> sys.stderr, "Error running program"
					print >> sys.stderr, "command line:"
					txt = ''
					for c in command:
						txt += c + ' '
					print >> sys.stderr, txt
					print >> sys.stderr, 'output:'
					output.seek(0)
					for l in output.readlines():
						print >> sys.stderr, l,
					output.close()
					sys.exit(-2)
				
				# now we process the output to find out wether the run was a success
				# (interface lambda_m reached) or a failure (interface lamda_f reached)
				output.seek(0)
				for line in output.readlines():
					words = line.split()
					if len(words) > 1:
						if words[1] == 'FFS' and words[2] == 'final':
							data = [w for w in words[4:]]
				op_names = data[::2]
				op_value = data[1::2]
				op_values = {}
				for ii, name in enumerate(op_names):
					op_values[name[:-1]] = float(op_value[ii][:-1])
				
				# now op_values is a dictionary representing the status of the final
				# configuration. We now need to find out wether it is a success or
				# a failure.
				# print op_values, 'op_values["%s"] %s %s' % (lambda_m_name, lambda_m_compar, str(lambda_m_value)), 'op_values["%s"] %s %s' % (lambda_f_name, lambda_f_compar, str(lambda_f_value))
				success = eval ('op_values["%s"] %s %s' % (lambda_m_name, lambda_m_compar, str(lambda_m_value)))
				#failure = eval ('op_values["%s"] %s %s' % (lambda_f_name, lambda_f_compar, str(lambda_f_value)))
				failure = eval ('op_values["%s"] %s %s' % (lambda_f_name, lambda_f_compar, str(lambda_pf_value)))
				if success:
					with success_lock:
						success_count.value += 1
						success_from[conf_index] += 1
						success_pruning_count.value += 1
						shutil.copy (last_conf, success_pattern + str(success_count.value))
						time_succeeding.value += (time.time() - mytime)
					log ("SUCCESS from PRUNING: worker %d: starting from conf_index %d and seed %d" % (idx, conf_index, seed))
				else:
					log ("FAILURE from PRUNING: worker %d: starting from conf_index %d and seed %d" % (idx, conf_index, seed))
				os.remove (last_conf)
		else:
			# do undetermined
			txt = ''	
			output.seek(0)
			for l in output.readlines():
				txt += l
			log ("UNDETERMINED: worker %d: starting from conf_index %d and seed %d\n%s" % (idx, conf_index, seed, txt))
			with undetermined_lock:
				undetermined_count.value += 1
				if (keep_undetermined):
					shutil.copy (last_conf, undetermined_pattern + str (undetermined_count.value))
				else:
					os.remove (last_conf)
	log ("Enough processes are started. Worker %d returning" % idx)

# timer function: it spits out things
def timer ():
	log ("Timer started at %s" % (time.asctime(time.localtime())))
	itime = time.time()
	while True:
		time.sleep (10)
		now = time.time()
		with success_lock:
			ns = success_count.value - initial_success_count
			na = attempt_count.value 
			if (ns > 0):
				log ("Timer: at %s: successes: %d (%d from pruning), attemps: %d: prob: %g time per success: %g (%g sec) (%d%% succeeding)" % (time.asctime(time.localtime()), ns, success_pruning_count.value, na, float(ns)/na, (now-itime)/float(ns), now - itime, int(100 * time_succeeding.value / ((now - itime) * float(ncpus)))))
			else:
				log ("Timer: at %s, no successes yet over %d attempts" % (time.asctime(time.localtime()), na))

if __name__ == '__main__':
	processes = []
	for i in range (ncpus):
		p = Process(target=f, args=(i,))
		processes.append(p)
	
	tp = Process(target=timer)
	tp.start()
	log ("starting processes...")
	for p in processes:
		p.start()
	
	log ("waiting for processes to finish")
	for p in processes:
		p.join()
	
	log ("Terminating timer")
	tp.terminate() # terminate timer

	nsuccesses = success_count.value - initial_success_count
	# print >> sys.stderr, "nstarted: %d, nsuccesses: %d success_prob: %g" % (nstarted, nsuccesses, nsuccesses/float(nstarted))
	log ("## log of successes probabilities from each starting conf")
	log ("conf_index nsuccesses nattempts prob")
	for k, v in enumerate(success_from):
		txt = "%3d    %3d    %3d   " % (k, v, attempt_from[k])
		if attempt_from[k] > 0:
			txt += '%g' % (float(v)/float(attempt_from[k]))
		else:
			txt += 'NA'
		log(txt)
	log ("# SUMMARY")
	log ("# nsuccesses: %d (of which %d from pruning) nattempts: %d success_prob: %g (%g accounting for pruning), undetermined: %d" % (nsuccesses, success_pruning_count.value, attempt_count.value, nsuccesses / float(attempt_count.value), ((nsuccesses - success_pruning_count.value + success_pruning_count.value / (1. - prune_prob)) /float(attempt_count.value)), undetermined_count.value))
