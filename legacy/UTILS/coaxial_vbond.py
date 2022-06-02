#!/usr/bin/env python

import sys
try:
    import numpy as np
except:
    import mynumpy as np
import base
import utils
import readers
import energies

MAX_FOUND = 10000

def strand_dstrand_energy(s1, ds1, ds2, box):
    energy = energies.coaxial(s1._nucleotides[0], ds1._nucleotides[0], box)        
    energy += energies.coaxial(s1._nucleotides[0], ds1._nucleotides[-1], box)
    energy += energies.coaxial(s1._nucleotides[0], ds2._nucleotides[0], box)        
    energy += energies.coaxial(s1._nucleotides[0], ds2._nucleotides[-1], box)

    energy += energies.coaxial(s1._nucleotides[-1], ds1._nucleotides[0], box)        
    energy += energies.coaxial(s1._nucleotides[-1], ds1._nucleotides[-1], box)
    energy += energies.coaxial(s1._nucleotides[-1], ds2._nucleotides[0], box)        
    energy += energies.coaxial(s1._nucleotides[-1], ds2._nucleotides[-1], box)

    if energy < 0.:
        for p in s1._nucleotides:
            for q in ds1._nucleotides + ds2._nucleotides:
                energy += energies.excluded_volume(p, q, box)
                if energy > 0.: return energy

    return energy

def bring_in_dstrand(s1, s2, initial, radius):
    R = utils.get_random_rotation_matrix()
    s1.rotate(R)
    s2.rotate(R, s1.cm_pos)
    
    start_from = [s1, s2][np.random.random_integers(0, 1)]._nucleotides[np.random.random_integers(-1, 0)].cm_pos
    new_pos = initial + utils.get_random_vector_in_sphere(radius)
    
    disp = new_pos - start_from
    s1.cm_pos = s1.cm_pos + disp
    s2.cm_pos = s2.cm_pos + disp

def main():
    if len(sys.argv) < 3:
        base.Logger.log("Usage is %s configuration topology [N_skip]" % sys.argv[0], base.Logger.CRITICAL)
        exit()
        
    if len(sys.argv) > 3:
        N_skip = int(sys.argv[3])
    else: N_skip = 0
        
    output = sys.argv[1] + ".vbond"
    
    np.random.seed(39823)
    
    l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
    s = l.get_system(N_skip=N_skip)
    box = s._box
    
    # we keep fixed strand 0 and 1 and we insert randomly strands 2 and 3 while 
    # keeping their internal degrees of freedom fixed
    diff = s._strands[1].cm_pos - s._strands[0].cm_pos
    diff = box * np.rint(diff / box)
    s._strands[1].cm_pos = s._strands[1].cm_pos - diff
    initial = (s._strands[0]._nucleotides[0].cm_pos + s._strands[1]._nucleotides[-1].cm_pos) / 2.
    radius = s._strands[0]._nucleotides[0].cm_pos - s._strands[1]._nucleotides[-1].cm_pos
    radius = base.CXST_RCHIGH + 0.5 * np.sqrt(np.dot(radius, radius))
    
    f = open(output, "w")
    V4pi = 4 * np.pi * (4 * np.pi * radius**3 / 3.)
    f.write("%e\n" % V4pi)
    f.close()
    
    diff = s._strands[3].cm_pos - s._strands[2].cm_pos
    diff = box * np.rint(diff / box)
    s._strands[3].cm_pos = s._strands[3].cm_pos - diff
    found = counter = 0
    
    while found < MAX_FOUND:
        bring_in_dstrand(s._strands[2], s._strands[3], initial, radius)
        
        E = strand_dstrand_energy(s._strands[2], s._strands[0], s._strands[1], box)
        if E < 2.: E += strand_dstrand_energy(s._strands[3], s._strands[0], s._strands[1], box)
        
        if E < 0.:
            print E
            f = open(output, "a")
            f.write("%d %e\n" % (counter, E))
            f.close()
            found += 1
                
        counter += 1
                
    s.print_crepy_output("prova2.mgl", same_colors=True)
                
    # cycle over the trajectory
    '''
    counter = 0
    while s:
    base.Logger.log("Analyzed configuration", base.Logger.INFO)
    s = l.get_system(N_skip=N_skip)
    counter += 1
    '''
    
    base.Logger.log("Output printed on '%s'" % output, base.Logger.INFO)
    

if __name__ == '__main__': 
#    import cProfile, pstats
#    cProfile.run('main()', 'pyProf')
#
#    p = pstats.Stats('pyProf')
#    p.strip_dirs().sort_stats("cumulative").print_stats(15)
    main()

