#!/usr/bin/env python

'''
Forward Flux sampling: Flux generator a-la-Tom

Flavio
'''

####################################
# edit here
####################################
# number of successes, can be changed via command line
desired_success_count = 100

# executable set-up
#precommand = 'mosrun -J 12'
executable = '../oxDNA'
input = 'input'
logfilename = 'ffs.log'
starting_conf = 'flux_initial.dat'
success_pattern = './success_'

# interfaces 
# interface lambda_{-1}
lambda_f_name = 'native'
lambda_f_value = 5

# interface lambda_{0}
lambda_0_name = 'native'
lambda_0_compar = '<='
lambda_0_value = 3

# interface lambda_{success}
lambda_s_name = 'dist'
lambda_s_compar = '>'
lambda_s_value = 5.
####################################

import os, sys, getopt
import subprocess as sp
import time, random as rnd, tempfile as tf
import shutil, glob
from multiprocessing import Process, Lock, JoinableQueue, Value, Array

def usage():
	print('usage: %s %s' % (sys.argv[0], '[-n <num_sucesses>] [-s <seed>] [-c <ncpus>] [-k <success_count>]'), file=sys.stderr) 

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
			initial_success_count = int(v) - 1
		else:
			print("Warning: option %s not recognized" % (k), file=sys.stderr)
except:
	print("Error parsing options", file=sys.stderr)
	sys.exit(3)

log_lock = Lock()
log_file = open (logfilename, 'w', 1)
def log(text):
	log_lock.acquire()
	log_file.write(text + '\n')
	if Verbose:
		print(text, file=sys.stdout)
	log_lock.release()

# check that we can write to the success pattern
try:
	checkfile = open (success_pattern + '0', 'w')
	checkfile.close()
	os.remove (success_pattern + '0')
except:
	print("could not write to success_pattern", success_pattern, file=sys.stderr)
	sys.exit(3)
	
success_lock = Lock()
success_count = Value ('i', initial_success_count)

# write the condition file
#condition_file = open('close.txt', 'w')
#condition_file.write("condition1 = {\n%s %s %s\n}\n" % (lambda_0_name, lambda_0_compar, str(lambda_0_value)))
#condition_file.close()

#condition_file = open('apart-bw.txt', 'w')
#condition_file.write("condition1 = {\n%s > %s\n}\n" % (lambda_f_name, str(lambda_f_value)))
#condition_file.close()
condition_file = open('apart-or-success.txt', 'w')
condition_file.write("condition1 = {\n%s > %s\n}\n" % (lambda_f_name, str(lambda_f_value)))
condition_file.write("condition2 = {\n%s %s %s\n}\n" % (lambda_s_name, lambda_s_compar, str(lambda_s_value)))
condition_file.close()

condition_file = open('apart-fw.txt', 'w')
condition_file.write("condition1 = {\n%s <= %s\n}\n" % (lambda_f_name, str(lambda_f_value)))
condition_file.close()

condition_file = open('both.txt', 'w')
condition_file.write("condition1 = {\n%s %s %s\n}\n" % (lambda_0_name, lambda_0_compar, str(lambda_0_value)))
condition_file.write("condition2 = {\n%s > %s\n}\n" % (lambda_f_name, str(lambda_f_value)))
condition_file.close()

# base command line; all features that need to be in the input file
# must be specified here
try:
	base_command = precommand.split()
except:
	base_command = []
base_command += [executable, input, 'print_energy_every=1e5', 'print_conf_every=1e6','no_stdout_energy=1','refresh_vel=0','restart_step_counter=0']
base_command_string = ''.join (str(w) + ' ' for w in base_command)
log("Main: COMMAND: " + base_command_string)

if not os.path.exists(starting_conf):
	log ("the input file provided (%s) does not exist. Aborting" % (starting_conf))
	sys.exit(-3)

if not os.path.exists(input):
	log ("the input file provided (%s) does not exist. Aborting" % (input))
	sys.exit(-3)
# check of the input file. If it contains an entry for the log file, spit out an error
inf = open (input, 'r')
log_found = False
for line in inf.readlines():
	words = line.split ("=")
	if len(words) >= 1:
		if words[0].lstrip().startswith("log_file"):
			log_found = True
if (log_found):
	print("\nERROR: This script does not work if \"log_file\" is set in the input file. Remove it! :)\n", file=sys.stderr)
	sys.exit (-2)
inf.close()

if not os.path.exists(executable):
	log ("the executable file provided (%s) does not exist. Aborting" % (executable))
	sys.exit(-3)

#rnd.seed (initial_seed)

log ("Main: STARTING new shooting for %d" % desired_success_count)
desired_success_count += initial_success_count

# this function does the work of running the simulation, identifying a
# success or a failure, and taking appropriate actions
def f(lidx):
	idx = lidx
	
	# the seed is the index + initial seed, and the last_conf has an index as well
	seed = initial_seed + idx
	myrng = rnd.Random()
	myrng.seed (seed)
	myrng.jumpahead(1)
	
	global success_count
	while success_count.value < desired_success_count:
		# we copy the initial configuration here
		my_conf = 'conf' + str(idx)
		shutil.copy(starting_conf, my_conf) 
		
		# edit the command to be launched
		my_base_command = base_command + ['conf_file=%s' % (my_conf), 'lastconf_file=%s' % (my_conf)]
		
		# initial equilibration
		command = my_base_command + ['seed=%d' % myrng.randint(1,50000), 'sim_type=MD', 'steps=1e5', 'refresh_vel=1']
		# open a file to handle the output
		output = tf.TemporaryFile ('r+', suffix=str(idx))
		
		log ("Worker %d: equilibration started " % idx)
		r = sp.call (command, stdout=output, stderr=sp.STDOUT)
		assert (r == 0)
		log ("Worker %d: equilibrated " % idx)

		# edit the command; we set to 0 the timer ONLY for the first time
		#command = my_base_command + ['ffs_file=apart-bw.txt', 'restart_step_counter=1', 'seed=%d' % myrng.randint(1,50000)]
		command = my_base_command + ['ffs_file=apart-or-success.txt', 'restart_step_counter=1', 'seed=%d' % myrng.randint(1,50000)]
	
		output.seek (0)
		
		# here we run the command
		# print command
		r = sp.call (command, stdout=output, stderr=sp.STDOUT)
		if r != 0:
			print("Error running program", file=sys.stderr)
			print("command line:", file=sys.stderr)
			txt = ''
			for c in command:
				txt += c + ' '
			print(txt, file=sys.stderr)
			print('output:', file=sys.stderr)
			output.seek(0)
			for l in output.readlines():
				print(l, end=' ', file=sys.stderr)
			output.close()
			sys.exit(-2)
		# now we process the output to find out wether the run was a complete success
		# (interface lambda_s reached) or a complete failure (interface lamda_f reached)
		output.seek(0)
		for line in output.readlines():
			words = line.split()
			if len(words) > 1:
				if words[1] == 'FFS' and words[2] == 'final':
					#print line,
					data = [w for w in words[4:]]
		op_names = data[::2]
		op_value = data[1::2]
		op_values = {}
		for ii, name in enumerate(op_names):
			op_values[name[:-1]] = float(op_value[ii][:-1])
		complete_failure = eval ('op_values["%s"] %s %s' % (lambda_f_name, '>', str(lambda_f_value)))
		complete_success = eval ('op_values["%s"] %s %s' % (lambda_s_name, lambda_s_compar, str(lambda_s_value)))
		if (complete_success):
			log ("Worker %d has reached a complete success: returning with the tail in between my legs");
			continue
		
		log ("Worker %d: reached Q_{-2}..." % idx)
		# now the system is far apart;

		while success_count.value < desired_success_count:
			# cross lamnda_{-1} going forwards
			output.seek (0)
			command = my_base_command + ['ffs_file=apart-fw.txt', 'seed=%d' % myrng.randint(1,50000)]
			r = sp.call (command, stdout=output, stderr=sp.STDOUT)
			assert (r == 0)
			log ("Worker %d: reached lambda_{-1} going forwards" % idx)
			
			# we hope to get to success
			output.seek(0)
			command = my_base_command + ['ffs_file=both.txt', 'seed=%d' % myrng.randint(1,50000)]
			r = sp.call(command, stdout=output, stderr=sp.STDOUT)
			assert (r == 0)
		
			# now we process the output to find out wether the run was a success
			# (interface lambda_m reached) or a failure (interface lamda_f reached)
			output.seek(0)
			for line in output.readlines():
				words = line.split()
				if len(words) > 1:
					if words[1] == 'FFS' and words[2] == 'final':
						#print line,
						data = [w for w in words[4:]]
			op_names = data[::2]
			op_value = data[1::2]
			op_values = {}
			for ii, name in enumerate(op_names):
				op_values[name[:-1]] = float(op_value[ii][:-1])
			
			# now op_values is a dictionary representing the status of the final
			# configuration.
			# print op_values, 'op_values["%s"] %s %s' % (lambda_m_name, lambda_m_compar, str(lambda_m_value)), 'op_values["%s"] %s %s' % (lambda_f_name, lambda_f_compar, str(lambda_f_value))
			success = eval ('op_values["%s"] %s %s' % (lambda_0_name, lambda_0_compar, str(lambda_0_value)))
			failure = eval ('op_values["%s"] %s %s' % (lambda_f_name, '>', str(lambda_f_value)))
			
			#print "EEE", op_values, success, failure #, 'op_values["%s"] %s %s' % (lambda_0_name, lambda_0_compar, str(lambda_0_value)), 'op_values["%s"] %s %s' % (lambda_f_name, '<', str(lambda_f_value))
			
			if success and not failure:
				with success_lock:
					success_count.value += 1
					shutil.copy (my_conf, success_pattern + str(success_count.value))
				log ("Worker %d: crossed interface lambda_{0} going forwards: SUCCESS" % idx)
				output.seek(0)
				#command = my_base_command + ['ffs_file=apart-bw.txt', 'restart_step_counter=1', 'seed=%d' % myrng.randint(1,50000)]
				command = my_base_command + ['ffs_file=apart-or-success.txt', 'restart_step_counter=1', 'seed=%d' % myrng.randint(1,50000)]
				r = sp.call(command, stdout=output, stderr=sp.STDOUT)
				assert (r == 0)
				output.seek(0)
				for line in output.readlines():
					words = line.split()
					if len(words) > 1:
						if words[1] == 'FFS' and words[2] == 'final':
							#print line,
							data = [w for w in words[4:]]
				op_names = data[::2]
				op_value = data[1::2]
				op_values = {}
				for ii, name in enumerate(op_names):
					op_values[name[:-1]] = float(op_value[ii][:-1])
				complete_failure = eval ('op_values["%s"] %s %s' % (lambda_f_name, '>', str(lambda_f_value)))
				complete_success = eval ('op_values["%s"] %s %s' % (lambda_s_name, lambda_s_compar, str(lambda_s_value)))
				if (complete_success):
					shutil.copy (my_conf, "full_success" + str(success_count.value))
					log ("Worker %d has reached a complete success: restarting from equilibration" % idx)
					break # this breakes the innermost while cycle
				else:
					log ("Worker %d: crossed interface lambda_{-1} going backwards after success" % idx)
			elif failure and not success:
				log ("Worker %d: crossed interface lambda_{-1} going backwards" % (idx))
			else:
				output.seek(0)
				#for l in output.readlines():
				#	print l,
				print(op_values) #, 'op_values["%s"] %s %s' % (lambda_0_name, lambda_0_compar, str(lambda_0_value)), 'op_values["%s"] %s %s' % (lambda_f_name, '<', str(lambda_f_value))
				log ("Worker %d: UNDETERMINED" % (idx))
				#sys.exit()

	os.remove(my_conf)


# timer function: it spits out things
def timer ():
	log ("Timer started at %s" % (time.asctime(time.localtime())))
	itime = time.time()
	while True:
		time.sleep (10)
		now = time.time()
		with success_lock:
			ns = success_count.value - initial_success_count
			if (ns > 1):
				log ("Timer: at %s: successes: %d, time per success: %g (%g sec)" % (time.asctime(time.localtime()), ns, (now-itime)/float(ns), now - itime))
			else:
				log ("Timer: at %s: no successes yet (at %d)" % (time.asctime(time.localtime()), success_count.value))

if __name__ == '__main__':
	processes = []
	for i in range (ncpus):
		p = Process(target=f, args=(i,))
		processes.append(p)
	
	tp = Process(target=timer)
	tp.start()
	log ("Main: Starting processes...")
	for p in processes:
		p.start()
	
	log ("Main: waiting for processes to finish")
	for p in processes:
		p.join()
	
	log ("Main: Terminating timer")
	tp.terminate() # terminate timer

	nsuccesses = success_count.value - initial_success_count
	# print >> sys.stderr, "nstarted: %d, nsuccesses: %d success_prob: %g" % (nstarted, nsuccesses, nsuccesses/float(nstarted))
	print("terminating processes", file=sys.stderr)
	
	log ("Main: nsuccesses: %d in this run" % (success_count.value - initial_success_count))

	# final coputation of the flux
	stime = 0
	confs = glob.glob(success_pattern + '*')
	for conf in confs:
		with open (conf, 'r') as f:
			t = int(f.readline().split('=')[1])
			stime += t
	log ("Main: average number of timesteps taken to reach a success (including possibly previous runs with the same pattern) (aka 1./flux): %g"  % (float(stime)/len(confs)))
	log ("Main: initial flux (includes previous runs if they were there): %g"  % (len(confs)/float(stime)))


