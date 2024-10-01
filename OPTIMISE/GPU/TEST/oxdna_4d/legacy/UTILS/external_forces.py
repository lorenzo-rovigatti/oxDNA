#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys

class Trap:
    types = ['string','trap']
    def __init__ (self, type, particle = None, pos0 = np.array([0,0,0]), rate = 0., stiff = 0., dir=np.array([0,0,1]), F0=0.):
        
        if type not in Trap.types:
            print >> sys.stderr, "ERROR: cannot create trap with type=%s" % (type)
            raise ValueError
            return None
        if particle is None:
            print >> sys.stderr, "ERROR: cannot create trap with particle=None" % (type)
            raise ValueError
            return None
        self.type = type
        self.particle = particle;
        self.pos0 = pos0;
        self.F0 = F0;
        self.rate = rate;
        self.stiff = stiff;
        self.dir = dir;

    def get_pot (self, pos, time = 0.):
        if self.type == 'trap':
            dr = pos - self.pos0
            return 1./2. * np.dot (dr, dr) * self.stiff 
        elif self.type == 'string':
            return - (self.F0 + self.rate * time) * np.dot (self.dir, pos)
        else:
            pass
    
    def get_force (self, pos, time = 0.):
        if self.type == 'trap':
            dr = pos - self.pos0
            return - self.stiff * dr
        elif self.type == 'string':
            return (self.F0 + self.rate * time) * self.dir
        else:
            pass

    def get_force_abs (self, pos):
        force = self.get_force (pos)
        return np.sqrt (np.dot (force, force))
    
    def __str__ (self):
        if self.type == 'trap':
            return "TRAP: part=%i, stiff=%lf, pos0=(%g,%g,%g)" % (self.particle, self.stiff, self.pos0[0], self.pos0[1], self.pos0[2])
        elif self.type =='string':
            return "STRING: part=%i, F0=%lf, rate=%g, dir=(%g,%g,%g)" % (self.particle, self.F0, self.rate, self.dir[0], self.dir[1], self.dir[2])
        else:
            pass


def parse_traps (inp = 'external.conf'):
    traps = []
    # parsing di external.conf
    try:
        external = open (inp)
    except:
        raise ValueError
        #%base.Logger.die ("Could not open %s. Dying now" % sys.argv[3])

    lines = external.readlines ()
    stripped = ""
    for line in lines:
        # remove trailing spaces
        line = line.lstrip ()
        if line.startswith ('#'):
            continue
        # remove comments
        i = line.find ('#')
        if i > 0:
            line = line[0:i] + '\n'
        # remove empty lines
        if line == '\n':
            continue
        
        stripped += line

    external.close ()

    strtraps = []
    while True:
        i1 = stripped.find('{')
        i2 = stripped.find('}')

        if i1 == -1 or i2 == -1:
            break
        
        strtraps.append (stripped[i1 + 1:i2])
        
        stripped = stripped[i2+1:]
    
    myrate = myF0 = mydir = mypos0 = mytype = mystiff = None
    for s in strtraps:
        lines = s.split('\n')
        for l in lines:
            words = [s.strip() for s in l.split('=')]
            if len(words) > 1:
                if words[0] == 'rate':
                    myrate = float(words[1])
                elif words[0] == 'particle':
                    mypart = int (words[1])
                elif words[0] == 'pos0':
                    mypos0 = np.array([float(i) for i in words[1].split(',')])
                elif words[0] == 'dir':
                    mydir = np.array([float(i) for i in words[1].split(',')])
                elif words[0] == 'stiff':
                    mystiff = float (words[1])
                elif words[0] == 'type':
                    mytype = words[1]
                elif words[0] == 'F0':
                    myF0 = float (words[1])
                else:
                    pass
         
        mytrap = Trap (mytype, particle=mypart, pos0=mypos0, rate=myrate, stiff=mystiff, F0=myF0, dir=mydir)
        base.Logger.log("Adding trap %s..." % (str(mytrap)), base.Logger.INFO)
        traps.append (mytrap)
    
    return traps

if __name__ == '__main__':
    if len(sys.argv) > 1:
        traps = parse_traps (sys.argv[1])
    else:
        traps = parse_traps ()
    print "## traps found:" 
    for t in traps:
        print t

