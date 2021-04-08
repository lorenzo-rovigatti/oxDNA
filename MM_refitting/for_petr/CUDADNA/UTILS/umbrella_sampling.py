import base, readers
import order_parameters as op
import copy, math, random, os, shutil, sys
import subprocess as sp
import time

try:
	import numpy as np
except:
	import mynumpy as np

# i assume that weights behave as follows:
# w = Weights ( dict ), where dict has the form {'name', pointer_to_object}.
# object must be able to return get_value (self, system, ops)
# same thing for OrderParameters, but now we want a dictionary
# with pointers to objects that can return get_state

class Timer:
    def __init__ (self):
        self.tot = 0.
        self.start = 0.
        self.end = 0.
        self.ison = False

    def reset (self):
        self.tot = 0.

    def on (self):
        if self.ison:
            raise ValueError ('timer is already on')
        self.start = time.time ()
        self.ison = True

    def off (self):
        if not self.ison:
            raise ValueError ('timer is already off')
        self.end = time.time ()
        self.tot = self.tot + self.end - self.start
        self.ison = False

    def __str__ (self):
        return "%g" % (self.tot)


def pycat (fsrc, fdst):
    dst = open (fdst, 'ab')
    src = open (fsrc, 'rb')
    shutil.copyfileobj (src, dst)
    src.close()
    dst.close()

class UmbrellaSampling:
    
    initError = "Could not create UmbrellaSampling object: " 

    def __init__ (self, dir='dir', exe='./oxDNA', input='input',
                  ordparams=None, weights=None, debug = False, seed = None,
                  output_interval = 1):
        
        if os.path.isdir (dir):
            self.dir = dir
        else:
            print >> sys.stderr, UmbrellaSampling.initError, dir, "does not exist. Aborting now"
            sys.exit (-1)
        
        self.exe = exe
        self.input = input
        self.ordparams = ordparams
        self.weights = weights
        self.debug = debug
        self.histo = op.Histo () # some histogram object
        self.seed = seed
        if self.debug and self.seed is None:
            print >> sys.stderr, "debugging: using seed = 1 for US"
            self.seed = 1 # to make it all deterministic
        random.seed (self.seed)
        
        self.total = 1.
        self.accepted = 1.
        self.exchanged = 1.
        self.output_interval = output_interval
        self.simtimer = Timer ()
        self.auxtimer = Timer ()
    
    def print_histo (self, outfile=sys.stdout, mode='w'):
        self.histo.print_n (outfile, mode)
    
    def print_info (self, fptr = sys.stdout):
        acc_ratio = self.accepted / float (self.total)
        exc_ratio = (self.exchanged + (self.total - self.accepted)) / (self.total)
        sim_timep = 100. * self.simtimer.tot / (self.simtimer.tot + self.auxtimer.tot)
        info = "%i %lf %lf %4.1f%s ##; it, acc_rate, exc_rate, sim_time_perc" % (int(self.total), acc_ratio, exc_ratio, sim_timep, '%')
        if isinstance (fptr, file):
            print >> fptr, info
        elif isinstance (fptr, str):
            try:
                out = open (fptr, 'a')
                print >> out, info
                out.close ()
            except:
                print >> sys.stderr, "could not write info to file ", fptr
        else:
            print >> sys.stderr, "Should not get here... Hoping for the best"   

    def run_test (self):
        # run a test run to see if the program behaves as expected
        if self.debug:
            print >> sys.stderr, "US: running test simulation..."

        tempdir = ''.join(random.choice("ABCDEFGHILMNOPQRSTUVZ0123456789") for x in range(20))
        shutil.copytree (self.dir, tempdir)
        
        ret = self.run (1)
        
        shutil.rmtree (self.dir)
        shutil.move (tempdir, self.dir)

        if ret != 0:
            print >> sys.stderr, "test failed. Aborting"
            sys.exit (-1)
        
        return
     
    def run (self, nsteps = 1):
        self.auxtimer.on ()
        # run the small simulation in dir
        os.chdir (self.dir)
        stored = 'stored' + ''.join (random.choice ("ABCDEFGHILMNOPQRSTUVZ0123456789") for x in range(5))
        # get relevant parameters from input
        try:
            inp = open (self.input, 'r')
        except:
            print >> sys.stderr, "Could not open input file", self.input, "in directory", self.dir, "Aborting now!"
            os.chdir ('..')
            return 1
        
        conf_file = ''
        topology = ''
        energy_file = 'energy.dat'
        last_conf_file = 'last_conf.dat'
        trajectory_file = 'trajectory.dat'
        for line in inp.readlines ():
            words = line.strip().split()
            if not words:
                continue
            if words[0] == 'conf_file':
                conf_file = words[2]
            elif words[0] == 'last_conf_file':
                last_conf_file = words[2]
            elif words[0] == 'trajectory_file':
                trajectory_file = words[2]
            elif words[0] == 'energy_file':
                energy_file = words[2]
            elif words[0] == 'topology':
                topology = words[2]
        inp.close ()
        
        if self.debug:
            print "Found input file with", conf_file, topology, energy_file, last_conf_file, trajectory_file
        try:
            shutil.copy (conf_file, stored)
        except:
            print >> sys.stderr, "Could not copy", conf_file, "in", stored
            os.chdir ('..')
            return 1
        
        #print self.input, conf_file, topology
        system = op.load_system (self.input, conf_file, topology)
        #we = self.weights.get_value (None, system, self.ordparams)
        statenow = self.ordparams.get_all_states (system)
        we = self.weights.get_all_factors (system, self.ordparams)

        for it in xrange (1, nsteps + 1):
            # copio last_conf in stored;
            shutil.copy (stored, conf_file)
            try:
                os.remove (trajectory_file)
                os.remove (energy_file)
            except:
                pass
            
            self.auxtimer.off ()
            self.simtimer.on ()
            
            myinput = sp.Popen ([self.exe, self.input], stdout=sp.PIPE, stderr=sp.PIPE)
            
            # check the return value. if different from 0, die badly
            if myinput.wait () != 0:
                print >> sys.stderr, "execute program died unexpectedely. Not following it all the way down to hell, but expect more of these messages..."
                os.chdir  ('..')
                break
            
            self.simtimer.off ()
            self.auxtimer.on ()
            
            # calcolo del peso sul nuovo sistema
            we_old = we
            newsystem = op.load_system (self.input, last_conf_file, topology)
            we = self.weights.get_value (None, newsystem, self.ordparams)
            
            # Metropolis:
            if random.random () < (float(we) / float(we_old)):
                if self.debug:
                    print we_old, '--->', we, "accepting"
                # if accepted, we cat energy and trajectory to the current
                # directory
                shutil.copy (last_conf_file, stored)
                system = newsystem
                self.accepted += 1.
            else:
                if self.debug:
                    print we_old, '--->', we, "refusing"
                we = we_old
            
            if self.output_interval is not None:
                if it % self.output_interval == 0:
                    pycat (stored, '../' + trajectory_file)
            
            # change the seed in the input file
            inp = open (self.input, 'r')
            ret = "seed = %d\n" % (random.randrange (0,50000))
            for line in inp.readlines ():
                if line.find ('seed') == -1:
                    ret += line
            
            inp.close ()
            inp = open (self.input, 'w')
            inp.write (ret)
            inp.close ()
            
            self.total += 1
            
            # save old state to see if it has changed
            oldstate = {}
            for key, val in statenow.items():
                oldstate[key] = val
            
            statenow = self.ordparams.get_all_states (system)
            print statenow, self.weights.get_all_factors (system, self.ordparams)
            
            if not statenow == oldstate:
                self.exchanged += 1
            
            # update histogram
            self.histo.update (statenow, we)
            # end of cycle
        try:
            os.remove (stored)
        except:
            pass
        
        self.auxtimer.off ()
        os.chdir ('..')
        return 0

'''
TODO:
    c) make it check for the presence of energy and trajectory
    e) (MAYBE) make it build the directory instead of having it. We do need conf
    and topology... or MAYBE build the input
    f) make projection of histograms along given order parameters
    m) run offline
    j) intercept signals maybe?
    b) clean up code as a whole
    h) timings? We may want to check not to have the python thing slow down
       the whole program too much
   #
'''


