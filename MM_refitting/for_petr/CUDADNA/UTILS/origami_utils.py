import base
import numpy as np
import sys
import subprocess
import pickle
import os
#import cgkit.cgtypes.mat3 as mat3

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")

def norm(vec):
    return vec / np.sqrt(np.dot(vec,vec))

def parse_scanfile(infile):
    try:
        f = open(infile, "r")
    except IOError:
        base.Logger.log("could not open file %s" % infile, base.Logger.CRITICAL)
        sys.exit()
    for line in f.readlines():
        if line.startswith('begin_vb'):
            begin_vb = int(line.split()[2])
        elif line.startswith('end_vb'):
            end_vb = int(line.split()[2])
        elif line.startswith('vh_list'):
            vh_list = [int(x) for x in (line.split()[2]).split(",")]
        elif line.startswith('trim'):
            trim = int(line.split()[2])
    try:
        return begin_vb, end_vb, vh_list, trim
    except NameError:
        base.Logger.log("error while reading scan file %s, dying now" % infile, base.Logger.CRITICAL)
        sys.exit()

def get_scaffold_index(system):
    # find the index of the scaffold strand in the system
    strand_lengths = [strand.get_length() for strand in system._strands]
    return strand_lengths.index(max(strand_lengths))

def get_bb_midpoint(system, strand, n_index, interaction_list):
    base.Logger.log("origami_utils.get_bb_midpoint: deprecated function, use origami_utils.Origami.get_bb_midpoint", base.Logger.WARNING)
    # get midpoint vector between 2 hybridised bases
    r1 = strand._nucleotides[n_index].get_pos_base()
    r2 = system._nucleotides[interaction_list[strand._nucleotides[n_index].index]].get_pos_base()
    vec = (r1+r2)/2
    return vec

def get_nucleotide(vhelix, vbase, vhelix_indices):
    # find the system nucleotide index of a nucleotide given a position on the origami
    if vhelix % 2 == 0:
        dir = 1
    else:
        dir = -1
    return vhelix_indices[vhelix] + vbase * dir

def parse_vh_data(filename, origami):
    # get vhelices data from file - format is either <auto,[number of vhelices]> or, if a region of an origami is to be analysed, <[region width],[list of starting nucleotide index for the region for each vhelix]>
    scaf_index = get_scaffold_index(origami._sys)
    vhelix_def_file = open(filename, "r")
    data = [x for x in vhelix_def_file.readline().replace("\n","").replace(" ","").split(",")]
    if data[0] == "auto":
        origami.num_vh = int(data[1])
        origami.width = origami._sys._strands[scaf_index].get_length() / origami.num_vh
        origami.vhelix_indices = []
        start_nuc_ind = -1
        for i in range(origami.num_vh):
            if i % 2 == 0:
                start_nuc_ind += 1
            else:
                start_nuc_ind += origami.width*2 - 1
            origami.vhelix_indices.append(start_nuc_ind)
    else:
        origami.width = int(data[0])
        origami.vhelix_indices = [int(x) for x in data[1:]]
        origami.num_vh = len(origami.vhelix_indices)
    base.Logger.log("using file data.vhd, %d virtual helices found" % origami.num_vh, base.Logger.INFO)

def print_arrow_debug_line(begin, end, file):
    if file:
        file.write("draw arrow {%f %f %f} {%f %f %f}\n" % (begin[0], begin[1], begin[2], end[0], end[1], end[2]))
    return 0

def open_arrow_debug_file(filename, type="w"):
    f_arrow = open(filename, type)
    f_arrow.write(
        "color Display Background white\n" +
        "set mName [mol new]\n" +
        "proc vmd_draw_arrow {mol start end} {\n" + 
        "# an arrow is made of a cylinder and a cone\n"
        "set middle [vecadd $start [vecscale 0.65 [vecsub $end $start]]]\n" + 
        "graphics $mol cylinder $start $middle radius 0.05\n" + 
        "graphics $mol cone $middle $end radius 0.15\n" + 
        "}\n")
    return f_arrow

def print_box_debug_line(vecs, file):
    if file:
        if len(vecs) == 4:
            file.write("draw box {%f %f %f} {%f %f %f} {%f %f %f} {%f %f %f}\n" % (vecs[0][0], vecs[0][1], vecs[0][2], vecs[1][0], vecs[1][1], vecs[1][2], vecs[2][0], vecs[2][1], vecs[2][2], vecs[3][0], vecs[3][1], vecs[3][2]))
        else:
            base.Logger.log("drawing boxes only works for 4 vertices at the moment", base.Logger.WARNING)

def open_box_debug_file(filename):
    f_boxes = open(filename, "w")
    f_boxes.write(
        "color Display Background white\n" +
        "set mName [mol new]\n" +
        "proc vmd_draw_box {mol vert1 vert2 vert3 vert4} {\n" + 
        "# a 'box' is a plane made of 2 triangles here\n" +
        "graphics $mol triangle $vert1 $vert2 $vert3\n" + 
        "graphics $mol triangle $vert1 $vert4 $vert3\n" + 
        "}\n")
    return f_boxes

class vhelix_vbase_to_nucleotide(object):
    # at the moment squares with skips in have entries in the dicts but with the nucleotide list empty (rather than having no entry) - I'm not sure whether or not this is desirable. It's probably ok
    def __init__(self):
        self._scaf = {}
        self._stap = {}
        self.nuc_count = 0 # record the nucleotide count, updated only after a whole strand is added
        self.strand_count = 0

    def add_scaf(self, vh, vb, strand, nuc):
        self._scaf[(vh, vb)] = (strand, nuc)
        
    def add_stap(self, vh, vb, strand, nuc):
        self._stap[(vh, vb)] = (strand, nuc)

    # these methods use a reference vhvb2n object to make the final vhvb2n object
    def add_scaf_strand(self, add_strand, reference, continue_join = False):
        count = 0
        size = len(self._scaf)
        for (vh, vb), [strand_ind, nuc] in reference._scaf.iteritems():
            if strand_ind == add_strand:
                self.add_scaf(vh, vb, self.strand_count, [x + self.nuc_count for x in nuc])
                count += len(nuc)
        self.nuc_count += count
        if len(self._scaf) == size:
            return 1
        else:
            if continue_join == False:
                self.strand_count += 1
            return 0

    def add_stap_strand(self, add_strand, reference, continue_join = False):
        count = 0
        size = len(self._stap)
        for (vh, vb), [strand_ind, nuc] in reference._stap.iteritems():
            if strand_ind == add_strand:
                self.add_stap(vh, vb, self.strand_count, [x + self.nuc_count for x in nuc])
                count += len(nuc)
        self.nuc_count += count
        if len(self._stap) == size:
            return 1
        else:
            if continue_join == False:
                self.strand_count += 1
            return 0

    def add_strand(self, add_strand, reference, continue_join = False):
        if self.add_scaf_strand(add_strand, reference, continue_join) and self.add_stap_strand(add_strand, reference, continue_join):
            base.Logger.log("while adding strand %s to vhelix_vbase_to_nucleotide object; either strand already present or strand not found in reference object" % add_strand, base.Logger.WARNING)
            return 1
        else:
            return 0
        
class Origami(object):
    
    def __init__(self, system=False, cad2cuda_file = False):
        self.width = 0
        self.num_vh = 0
        if system:
            self._sys = system
            self.interaction_list = [-1 for x in range(self._sys._N)]
        self.vhelix_indices = []
        self.vbase_indices = []
        self.vec_long = np.array([0., 0., 0.])
        self.vec_lat = np.array([0., 0., 0.])
        self._vhelix_pattern = False
        if cad2cuda_file:
            self.get_cad2cudadna(cad2cuda_file)
            # build list of complementary nucleotides (according to cadnano scheme)
            self.complementary_list = ["na" for x in range(self._sys._N)]
            for (vhelix, vbase), (strand1, nucs1) in self._cad2cudadna._scaf.iteritems():
                try:
                    (strand2, nucs2) = self._cad2cudadna._stap[vhelix, vbase]
                    for i in range(len(nucs1)):
                        # I'm pretty sure these should always line up for any possible insertion/deletion scheme
                        self.complementary_list[nucs1[i]] = nucs2[i]
                        self.complementary_list[nucs2[i]] = nucs1[i]
                except KeyError:
                    pass
                
            # build a list of the vhelix indices contained in the cadnano design
            for (vhelix, vbase) in iter(self._cad2cudadna._scaf):
                if vhelix not in self.vhelix_indices:
                    self.vhelix_indices.append(vhelix)
                if vbase not in self.vbase_indices:
                    self.vbase_indices.append(vbase)
            self.vhelix_indices.sort()
            self.vbase_indices.sort()
            
            # build a list of the occupied virtual bases in each virtual helix
            self.vh_vbase_indices = [[] for x in self.vhelix_indices]
            self.vh_vbase_indices_scaf = [[] for x in self.vhelix_indices]
            self.vh_vbase_indices_stap = [[] for x in self.vhelix_indices]
            self.vvib = [[] for x in self.vhelix_indices] # vhelix_vbase_indices_both
            # scaffold vbase occupation
            for (vh, vb) in iter(self._cad2cudadna._scaf):
                self.vh_vbase_indices_scaf[self.vhelix_indices.index(vh)].append(vb)
            for row in self.vh_vbase_indices_scaf:
                row.sort()
            # staples vbase occupation
            for (vh, vb) in iter(self._cad2cudadna._stap):
                self.vh_vbase_indices_stap[self.vhelix_indices.index(vh)].append(vb)
            for row in self.vh_vbase_indices_stap:
                row.sort()
            # occupation of both staple and scaffold strand
            for vhi in range(len(self.vh_vbase_indices_scaf)):
                for vb_scaf in self.vh_vbase_indices_scaf[vhi]:
                    if vb_scaf in self.vh_vbase_indices_stap[vhi]:
                        self.vvib[vhi].append(vb_scaf)
            for row in self.vvib:
                row.sort()
            # occupation of either staple or scaffold strand
            for vhi in range(len(self.vh_vbase_indices_scaf)):
                for vb_scaf in self.vh_vbase_indices_scaf[vhi]:
                    if vb_scaf not in self.vh_vbase_indices[vhi]:
                        self.vh_vbase_indices[vhi].append(vb_scaf)
                for vb_stap in self.vh_vbase_indices_stap[vhi]:
                    if vb_stap not in self.vh_vbase_indices[vhi]:
                        self.vh_vbase_indices[vhi].append(vb_stap)
            for row in self.vh_vbase_indices:
                row.sort()
            # nicer aliases
            self.vvi = self.vh_vbase_indices
            self.vvisc = self.vh_vbase_indices_scaf
            self.vvist = self.vh_vbase_indices_stap
            
            self.num_vh = len(self.vhelix_indices)
            self.scaf_index = get_scaffold_index(self._sys)
#            self.width = self._sys._strands[scaf_index].get_length() / self.num_vh
#            if self.width != len(self.vbase_indices):
#                pass #this warning got annoying base.Logger.log("not a rectangular origami!", base.Logger.WARNING)

        else:
            self._cad2cudadna = {}

    def update_system(self, system):
        self._sys = system
    
    def get_corners(self):
        if self._cad2cudadna == {}:
            base.Logger.log("get_corners: build cad2cudadna property first", base.Logger.CRITICAL)
            sys.exit()

        # make sure that neighbours in the list are neighbours in the origami
        a = self.get_nucleotides(self.vhelix_indices[0],self.vh_vbase_indices[0][0])[0]
        b = self.get_nucleotides(self.vhelix_indices[0],self.vh_vbase_indices[0][-1])[0]
        c = self.get_nucleotides(self.vhelix_indices[-1],self.vh_vbase_indices[-1][-1])[0]
        d = self.get_nucleotides(self.vhelix_indices[-1],self.vh_vbase_indices[-1][0])[0]
        return [a, b, c, d]
        
    def get_nucleotides(self, vhelix, vbase, type="default"):
        # tries to return scaffold strand nucleotide, failing that staple strand, failing that error
        if self._cad2cudadna:
            if type == "default":
                try:
                    strand, nucs = self._cad2cudadna._scaf[(vhelix, vbase)]
                except KeyError:
                    strand, nucs = self._cad2cudadna._stap[(vhelix, vbase)]
                return nucs
            elif type == "double":
                # return double strands
                strand, nucs1 = self._cad2cudadna._scaf[(vhelix, vbase)]
                strand, nucs2 = self._cad2cudadna._stap[(vhelix, vbase)]
                nucs = []
                nucs.extend(nucs1)
                nucs.extend(nucs2)
                return nucs
        else:
            base.Logger.log("no cadnano to cudadna file detected, using old and possibly wrong get_nucleotides function", base.Logger.WARNING)
            # find the system nucleotide index of a nucleotide given a position on the origami
            if vhelix % 2 == 0:
                dir = 1
            else:
                dir = -1
            return self.vhelix_indices[vhelix] + vbase * dir

    def get_vhelix_ds_length(self, vhi):
        # requires one continuous double strand - otherwise how is double strand length defined??
        nucleotide_count = 0
        for vb in self.vvib[vhi]:
            try:
                nucs = self.get_nucleotides(self.vhelix_indices[vhi], vb, type="double")
            except KeyError:
                continue
            for nuc in nucs:
                nucleotide_count += 1
        return nucleotide_count/2
            

    def get_flat_nucs(self, vh):
        vhi = self.vhelix_indices.index(vh)
        nucs_flat = []
        iterable = self.vh_vbase_indices[vhi]
        for vb in iterable:
            try:
                nucs = self.get_nucleotides(vh, vb)
            except:
                continue #pass
            nucs_flat.extend(nucs)

        return nucs_flat
    
    def vb2nuci(self, vhi, vb0):
        vh = self.vhelix_indices[vhi]
        nucs_count = []
        vbi = self.vvib[vhi].index(vb0)
        iterable = self.vvib[vhi][:vbi]
        for vb in iterable:
            try:
                nucs = self.get_nucleotides(vh,vb)
            except:
                continue
            nucs_count.extend(nucs)

        return len(nucs_count)

    def prepare_principal_axes_calc():
        # get the central nucleotides that will (may) be used as reference points for the principal axes calculations
        com = np.array([0.,0.,0.])
        for nuc in self._sys._nucleotides:
            com += nuc.cm_pos
        com /= self._sys._N

        displ = range(self._sys._N)
        for ii in range(self._sys._N):
            rr = self._sys._nucleotides[ii].cm_pos - com
            displ[ii] = np.dot(rr,rr)

        minnuc = displ.index(min(displ))

        self._sys.map_nucleotides_to_strands()
        if minnuc+1 < self._sys._N and self._sys._nucleotide_to_strand[minnuc] == self._sys._nucleotide_to_strand[minnuc+1]:
            minnucs = (minnuc, minnuc+1)
        elif minnuc > 0 and self._sys._nucleotide_to_strand[minnuc] == self._sys._nucleotide_to_strand[minnuc-1]:
            minnucs = (minnuc, minnuc-1)
        else:
            base.Logger.die("something went wrong when trying to find the reference nucleotides for the principal axes calculation")
        return minnucs
        
    def get_principal_axes(self, approxaxes, minnucs):
        print "origami_utils.py: Origami.get_principal_axes: unsupported function, dying"
        sys.exit()
        # find moment of inertia tensor I, then eigenvectors are the principal axes
        # first get centre of mass
        com = np.array([0.,0.,0.])
        for nuc in self._sys._nucleotides:
            com += nuc.cm_pos
        com /= self._sys._N

        # get global rotation, the rotation to ensure that a particular vector always lies on [1,0,0] for every configuration
        vecref = norm(self._sys._nucleotides[minnucs[0]].cm_pos - self._sys._nucleotides[minnucs[1]].cm_pos)
        globalrotcg = 0#mat3().fromToRotation([vecref[0], vecref[1], vecref[2]], [1.,0.,0.])
        globalrot = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
        # convert to numpy array
        for ii in range(3):
            for jj in range(3):
                globalrot[ii][jj] = globalrotcg[ii][jj]
        # find I wrt centre of mass
        I = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        for nuc in self._sys._nucleotides:
            # rotate system so that vecref along x axis
            rk = np.dot(globalrot,(nuc.cm_pos - com))
            I[0,0] += rk[1]*rk[1] + rk[2]*rk[2]
            I[1,1] += rk[0]*rk[0] + rk[2]*rk[2]
            I[2,2] += rk[0]*rk[0] + rk[1]*rk[1]
            I[0,1] -= rk[0]*rk[1]
            I[0,2] -= rk[0]*rk[2]
            I[1,2] -= rk[1]*rk[2]
            I[1,0] -= rk[0]*rk[1]
            I[2,0] -= rk[0]*rk[2]
            I[2,1] -= rk[1]*rk[2]

        eigvals, eigvecs = np.linalg.eig(I)

        # order eigenvectors by size of eigenvalue
        i1 = np.where(eigvals == max(eigvals))[0][0]
        i3 = np.where(eigvals == min(eigvals))[0][0]
        i2 = 3 - i1 - i3
        v1 = eigvecs[:,i1]
        v2 = eigvecs[:,i2]
        v3 = eigvecs[:,i3]

        # order eigenvectors by how much they overlap with the approximate axes
        eigvecsl = [eigvecs[:,i] for i in range(len(eigvecs))] # make a list of 1d arrays from a 2d array
        eigvecs_ordered = []
        for ii in range(len(approxaxes)):
            res = range(len(eigvecsl))
            for jj in range(len(eigvecsl)):
                res[jj] = abs(np.dot(approxaxes[ii], norm(eigvecsl[jj]))) # to allow antiparallel vectors to count as aligned - not sure whether this is the best method

            eigvecs_ordered.append(eigvecsl.pop(res.index(max(res))))

        return np.array([eigvecs_ordered[0], eigvecs_ordered[1], eigvecs_ordered[2]])
    
    def get_1d_vecs(self, discard_unbonded = True):
        # returns orthonormal vectors based on the orientation of the strand
        # "invariant" co-ordinates x', y', z'. x': along double helix; yy: average backbone-base vector for the 0th strand; z': x' cross yy; y': z' cross x'

        trim = 3 # ignore ends of the double strand; for very short strands this may not be desirable
        
        vh = self.vhelix_indices[0]
        nucs_flat = self.get_flat_nucs(vh)
        
        # x'
        if len(nucs_flat) > trim * 2:
            n1 = nucs_flat[trim]
            n2 = nucs_flat[-trim]
        else:
            n1 = nucs_flat[0]
            n2 = nucs_flat[0]
            
        # check they exist # surely we already know they do??? 13/06/12
        bbm1 = self.get_bb_midpoint(n1)
        bbm2 = self.get_bb_midpoint(n2)
        if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
            xprime = bbm2 - bbm1

        # yy
        yy = np.zeros(3)
        if len(nucs_flat) > trim * 2:
            nuciter = nucs_flat[trim:-trim]
        else:
            nuciter = nucs_flat
        for nuci in nuciter:
            nuc = self._sys._nucleotides[nuci]
            yy += nuc._a1

        zprime = np.cross(xprime, yy)
        yprime = np.cross(zprime, xprime)

        xprime = norm(xprime)
        yprime = norm(yprime)
        zprime = norm(zprime)
        m = np.array( [xprime,yprime,zprime] )
        
        return m
        
    def get_3d_vecs(self, vhelix_extent = False, discard_unbonded = True):
        # returns orthonormal vectors based on the orientation of the origami.
        # z: average vector along helices; y: average vec between vertical helices; x: y cross z
        # vhelix_extent currently unused.... uses extent of double helices (assuming there is one duplex per virtual helix)
        # assumes the skips are lined up across vhelices - in particular it will presumably break if there is a vb with a skip adjacent (in vhelix dimension) to a vb that represents the end of a double helix
        
        if not self._vhelix_pattern:
            base.Logger.log("origami_utils: Origami.get_3d_vecs(): no virtual helix pattern found; either missing virt2nuc file in this directory, or virt2nuc file is out of date; Aborting now", base.Logger.CRITICAL)
            sys.exit(1)

        y = np.zeros(3, dtype = "float64")
        z = np.zeros(3, dtype = "float64")
        for vhi in range(len(self.vhelix_indices)):
            vh = self.vhelix_indices[vhi]
            # z: vec along helices
            vb1 = self.vvib[vhi][0]
            vb2 = self.vvib[vhi][-1]
            n1 = self.get_nucleotides(vh, vb1)[0]
            n2 = self.get_nucleotides(vh, vb2)[-1]
            if (self.interaction_list[n1] == -1 or self.interaction_list[n2] == -1) and discard_unbonded:                
                base.Logger.log("unbonded base when calculating interaction between %d,%d and %d,%d; interaction will not be used" % (vh, 0, vh, self.width - 1), base.Logger.WARNING)
            else:
                bbm1 = self.get_bb_midpoint(n1)
                bbm2 = self.get_bb_midpoint(n2)
                if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
                    z += bbm2 - bbm1

            # y: vec across virtual helices that are arranged vertically on the honeycomb lattice
            row = self._vhelix_pattern[vh][0]
            col = self._vhelix_pattern[vh][1]
            # check for a virtual helix above the current one
            result = [vh1 for vh1, (row1,col1) in self._vhelix_pattern.iteritems() if row1 == row - 1 and col1 == col]
            if len(result) > 0:
                # a virtual helix above the current one was found
                vhin = result[0]
                vb1 = max(self.vvib[vhi][0], self.vvib[vhin][0])
                vb2 = min(self.vvib[vhi][-1], self.vvib[vhin][-1])
                for vb in (vb1, vb2):
                    lnucs = len(self.get_nucleotides(vh, vb))
                    if lnucs > 0:
                        if vb in self.vvib[vhin]:
                            n1 = self.get_nucleotides(self.vhelix_indices[vhi],vb)[0]
                            n2 = self.get_nucleotides(self.vhelix_indices[vhin],vb)[0]
                            bbm1 = self.get_bb_midpoint(n1)
                            bbm2 = self.get_bb_midpoint(n2)
                            if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
                                y += bbm2 - bbm1

        x = np.cross(y,z)
        # redefine y to make sure it's perpendicular to z
        y = np.cross(z,x)

        x /= np.sqrt(np.dot(x,x))
        y /= np.sqrt(np.dot(y,y))
        z /= np.sqrt(np.dot(z,z))
        m = np.array([x,y,z])
        
        return m
        
    def get_plane_vecs(self, discard_unbonded = True):
        # first find average bbm-bbm vector along helices
        av_long = np.zeros(3, dtype = "float64")
        for vhi in range(len(self.vhelix_indices)):
            vh = self.vhelix_indices[vhi]
            vb1 = self.vvib[vhi][0]
            vb2 = self.vvib[vhi][-1]
            n1 = self.get_nucleotides(vh, vb1)[0]
            n2 = self.get_nucleotides(vh, vb2)[-1]
            '''
            for vbi in range(len(self.vvib[vhi][:-1])):
                vb = self.vvib[vhi][vbi]
                nucs = self.get_nucleotides(vh, vb)
                for nuci in range(len(nucs)):
                    # find the next nucleotide along
                    n1 = nucs[nuci]
                    try:
                        n2 = nucs[nuci+1]
                    except IndexError:
                        ii = vbi
                        while ii < len(self.vvib[vhi]):
                            try:
                                n2 = self.get_nucleotides(vh, vb+1)[0]
                                break
                            except:
                                ii += 1

            '''
            if (self.interaction_list[n1] == -1 or self.interaction_list[n2] == -1) and discard_unbonded:
                base.Logger.log("unbonded base when calculating interaction between %d,%d and %d,%d; interaction will not be used" % (vh, vb1, vh, vb2), base.Logger.WARNING)
            else:
                bbm1 = self.get_bb_midpoint(n1)
                bbm2 = self.get_bb_midpoint(n2)
                if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
                    av_long += bbm2 - bbm1

        # next find average bbm-bbm vector across helices
        av_lat = np.zeros(3, dtype = "float64")
        for vb in self.vbase_indices:
            # get a list of all vhelices with that vbase occupied
            vhs = []
            for x in iter(self._cad2cudadna._scaf):
                if x[1] == vb:
                    vhs.append(x[0])
            vhs.sort()
            if len(vhs) < 2:
                continue
            # if there is a pattern of skips and loops that is not the same across vhelices, this code tries to cope with it but it's not the most rigorous way...
            nucs = self.get_nucleotides(vhs[0], vb)
            for nuci in range(len(nucs)):
                n1 = nucs[nuci]
                skip = False
                for jj in range(len(vhs))[::-1]:
                    try:
                        n2 = self.get_nucleotides(vhs[jj], vb)[nuci]
                        break
                    except:
                        if jj in [0,1]:
                            skip = True
                            break
                        else:
                            pass
                if skip:
                    continue
                if (self.interaction_list[n1] == -1 or self.interaction_list[n2] == -1) and discard_unbonded:
                    base.Logger.log("unbonded base when calculating interactions between %d,%d and %d,%d; interaction will not be used" % (0, vb, self.num_vh - 1, vb), base.Logger.WARNING)
                else:
                    bbm1 = self.get_bb_midpoint(n1)
                    bbm2 = self.get_bb_midpoint(n2)
                    if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
                        av_lat += bbm2 - bbm1

        av_long /= np.sqrt(np.dot(av_long, av_long))
        av_lat /= np.sqrt(np.dot(av_lat, av_lat))

        self.vec_long = av_long
        self.vec_lat = av_lat
        
    def get_bb_midpoint(self, n_index):
        # get midpoint vector between 2 bases that are hybridised according to cadnano scheme
        if self.complementary_list[n_index] == "na":
            return False
        else:
            r1 = self._sys._nucleotides[n_index].get_pos_base()
            r2 = self._sys._nucleotides[self.complementary_list[n_index]].get_pos_base()
            vec = (r1+r2)/2
            return vec
        
    def get_backback_midpoint(self, n_index):
        # get midpoint vector between 2 bases that are hybridised according to cadnano scheme
        if self.complementary_list[n_index] == "na":
            return False
        else:
            r1 = self._sys._nucleotides[n_index].get_pos_back()
            r2 = self._sys._nucleotides[self.complementary_list[n_index]].get_pos_back()
            vec = (r1+r2)/2
            return vec
        
    def get_bb_vec(self, n_index):
        # get midpoint vector between 2 bases that are hybridised according to cadnano scheme
        if self.complementary_list[n_index] == "na":
            return False
        else:
            r1 = self._sys._nucleotides[n_index].get_pos_base()
            r2 = self._sys._nucleotides[self.complementary_list[n_index]].get_pos_base()
            vec = r1-r2
            return vec
        
    def get_backback_vec(self, n_index):
        # get midpoint vector between 2 bases that are hybridised according to cadnano scheme
        if self.complementary_list[n_index] == "na":
            return False
        else:
            r1 = self._sys._nucleotides[n_index].get_pos_back()
            r2 = self._sys._nucleotides[self.complementary_list[n_index]].get_pos_back()
            vec = r1-r2
            return vec
        
    def get_flat_units(self, discard_unbonded):
        # find width, height units of the origami if it had been laid out flat
        flat_width = 0
        flat_height = 0
        for vh in range(self.num_vh):
            for vb in range(self.width - 1):
                n1 = self.get_nucleotide(vh, vb)
                n2 = self.get_nucleotide(vh, vb + 1)
                if (self.interaction_list[n1] != -1 and self.interaction_list[n2] != -1) or not discard_unbonded:
                    dist = self.get_bb_midpoint(n1) - self.get_bb_midpoint(n2)
                dist = np.sqrt(np.dot(dist, dist))
                flat_width += dist
        flat_width /= self.num_vh * (self.width - 1)
        
        for vh in range(self.num_vh - 1):
            for vb in range(self.width):
                n1 = self.get_nucleotide(vh, vb)
                n2 = self.get_nucleotide(vh + 1, vb)
                if (self.interaction_list[n1] != -1 and self.interaction_list[n2] != -1) or not discard_unbonded:
                    dist = self.get_bb_midpoint(n1) - self.get_bb_midpoint(n2)
                dist = np.sqrt(np.dot(dist, dist))
                flat_height += dist
        flat_height /= (self.num_vh - 1) * self.width
        return flat_width, flat_height

    def get_com(self):
        # unsupported/may no longer work - e.g. take care in using self.width and get_bb_midpoint
        com = np.zeros(3)
        for vh in range(self.num_vh):
            for vb in range(self.width):
                com += self.get_bb_midpoint(self.get_nucleotide(vh, vb))
        com /= self.num_vh * self.width
        return com

    def get_h_bond_list(self, infile, conffile, conf_num):
        system = self._sys
        system.map_nucleotides_to_strands()
        try:
            open(infile)
        except:
            base.Logger.log("unable to find file %s, exit" % infile, base.Logger.CRITICAL)
            sys.exit()
        launchargs = [PROCESSDIR + 'output_bonds',infile,conffile,str(conf_num)]
        myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        system.read_H_bonds(myinput.stdout.readlines())

        self.interaction_list = [[] for x in range(self._sys._N)]
        for nucleotide in system._nucleotides:
            for i in nucleotide.interactions:
                self.interaction_list[nucleotide.index].append(i)

    def get_cad2cudadna(self, infile):
        f = open(infile, "r")
        data = pickle.load(f)
        if type(data) is tuple:
            self._cad2cudadna = data[0]
            self._vhelix_pattern = data[1]
        else:
            self._cad2cudadna = data
        f.close()

    def update_system(self, system):
        self._sys = system

    def map_base_to_vbase(self, vhelices_len, begin_bases, end_bases):
        # not used right now
        # returns virtual base AND the index of the nucleotide at that virtual base (ie to deal with the possiblity of a loop)
        self._base_to_vbase = [[] for x in range(vhelices_len)]
        for vhelix in range(vhelices_len):
            for mybase in range(begin_bases[vhelix], end_bases[vhelix]+1):
                current_base = -1
                current_vbase = -1
                while current_base < mybase:
                    current_vbase += 1
                    try:
                        strand, nucs = self._cad2cudadna._stap[(vhelix, current_vbase)]
                    except KeyError:
                        nucs = []
                    for i in range(len(nucs)):
                        current_base += 1
                        if current_base == mybase:
                            self._base_to_vbase[vhelix].append((current_vbase, i))
                            break
                    if current_vbase > (end_bases[vhelix]+1) * 2:
                        base.Logger.die("Unable to find virtual base in map_base_to_vbase, dying now")

    def base_to_vbase(self, vhelix, base):
        # not used right now
        return self._base_to_vbase[(vhelix, base)]
    
    def vbase_to_base(self, vhelix, vbase):
        # not used right now
        # return base corresponding to first base at a particular virtual base i.e. taking into account skips and loops (counting from the left of a cadnano diagram)
        base = 0
        for i in range(vbase - 1):
            try:
                strand, nucs = self._cad2cudadna._stap[(vhelix, i)]
            except KeyError:
                nucs = []
            base += len(nucs)
        base += 1
        return base

    def get_holliday_junctions(self, single=False):
        # assumes no insertions/deletions when using cad2cudadna
        # set single=True to get a list of single crossovers rather than holliday junctions ie pairs of crossovers
        # if somehow there were two holliday junctions between the same pairs of vbases (one for scaffold and one for staple strand), the 2nd would not be detected here
        # returns list of hjs in format (vh, vb, vh_neighbour, vb_next)
        vh_pattern = self._vhelix_pattern
        vh_neighbours_below = [[] for x in range(max(vh_pattern.keys())+1)]
        for vh, (row,col) in vh_pattern.iteritems():
            neighbours = [vh1 for vh1, (row1,col1) in vh_pattern.iteritems() if (row1 == row and (col1 - col) == 1) or ((row1 - row) == 1 and col1 == col)]
            for vh1 in neighbours:
                vh_neighbours_below[vh].append(vh1)
                
        hjs = []
        skip = False
        lone_crossover = 0
        for vhi in range(len(self.vhelix_indices)):
            vh = self.vhelix_indices[vhi]
            for ii in range(len(self.vvib[vhi])):
                if skip:
                    skip = False
                    continue
                vb = self.vvib[vhi][ii]
                try:
                    vhn = vh_neighbours_below[vh][0]
                except IndexError:
                    continue
                try:
                    # for a crossover we have no insertions/deletions - here we make sure that this is the case, because if there are 0 or > 1 nucleotides on this square, a ValueError will be thrown
                    (cstap, [cnuc]) = self._cad2cudadna._stap[(vh,vb)]
                    (nstap, [nnuc]) = self._cad2cudadna._stap[(vhn,vb)]
                except (KeyError, ValueError):
                    continue
                # if there is a crossover between the squares (vh,vb) and (vhn,vb), the strand will be the same and the nucleotides will be consecutive
                if cstap == nstap and abs(cnuc - nnuc) == 1:
                    if single:
                        hjs.append([vh, vhn, vb])
                    else:
                        this_hj = [vh, vb, vhn]
                        # we expect the staple crossovers to be in pairs
                        vb += 1
                        try:
                            (cstap, [cnuc]) = self._cad2cudadna._stap[(vh,vb)]
                            (nstap, [nnuc]) = self._cad2cudadna._stap[(vhn,vb)]
                        except (KeyError, ValueError):
                            lone_crossover += 1
                            continue
                        if cstap == nstap and abs(cnuc - nnuc) == 1:
                            this_hj.append(vb)
                            hjs.append(this_hj)
                            skip = True
                        else:
                            lone_crossover += 1
                            continue
                else:
                    # if no staple crossover, we check for a scaffold crossover
                    try:
                        (cscaf, [cnuc]) = self._cad2cudadna._scaf[(vh,vb)]
                        (nscaf, [nnuc]) = self._cad2cudadna._scaf[(vhn,vb)]
                    except (KeyError, ValueError):
                        continue
                    if cscaf == nscaf and abs(cnuc - nnuc) == 1:
                        if single:
                            hjs.append([vh, vhn, vb])
                        else:
                            this_hj = [vh, vb, vhn]
                            # we expect the scaffold crossovers to be in pairs
                            vb += 1
                            try:
                                (cscaf, [cnuc]) = self._cad2cudadna._scaf[(vh,vb)]
                                (nscaf, [nnuc]) = self._cad2cudadna._scaf[(vhn,vb)]
                            except (KeyError, ValueError):
                                lone_crossover += 1
                                continue
                            if cscaf == nscaf and abs(cnuc - nnuc) == 1:
                                this_hj.append(vb)
                                hjs.append(this_hj)
                                skip = True
                            else:
                                lone_crossover += 1
                                continue
                    
        if not single and lone_crossover > 0:
            base.Logger.log("%d lone crossovers found, they will not be used in analysis" % lone_crossover, base.Logger.INFO)
        return hjs

    def get_vhelix_neighbours(self, type):
        # checked for square lattice, not yet double checked for honeycomb lattice
        if type not in ("he", "sq"):
            base.Logger.log("origami_utils: get_vhelix_neighbours(): error while building neighbour list; unknown lattice type %s, dying now" % type, base.Logger.CRITICAL)
            sys.exit()
        vh_pattern = self._vhelix_pattern
        vh_neighbour_list = [[] for x in range(max(vh_pattern.keys())+1)]
        for vh, (row,col) in vh_pattern.iteritems():
            if type == "sq":
                neighbours = [vh1 for vh1, (row1,col1) in vh_pattern.iteritems() if (row1 == row and abs(col1 - col) == 1) or (abs(row1 - row) == 1 and col1 == col)]
            elif type == "he":
                if vh % 2 == 0:
                    neighbours = [vh1 for vh1, (row1,col1) in vh_pattern.iteritems() if (row1 == row and abs(col1 - col) == 1) or (row1 == row - 1 and col1 == col)]
                if vh % 2 == 1:
                    neighbours = [vh1 for vh1, (row1,col1) in vh_pattern.iteritems() if (row1 == row and abs(col1 - col) == 1) or (row1 == row + 1 and col1 == col)]
            for vh1 in neighbours:
                vh_neighbour_list[vh].append(vh1)

        return vh_neighbour_list

    def get_local_twist(self, vh_id, vvib_id, type, conf):
        # find local twist between given vhelix, vbase and the next one along
        # assumes no skip/loop / insertion/deletion
        ss = self._sys
        n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
        n2 = self.complementary_list[n1]
        try:
            n3 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id+1])[0]
        except KeyError:
            base.Logger.die("origami_utils.Origami.get_local_twist: KeyError - probably the virtual base index argument was too high")
        n4 = self.complementary_list[n3]
        # get the 3 normalised vectors we need: the base-base (or back-back) vectors and the vector from one base-base (or back-back) midpoint to the other
        if type == "base":
            r1 = ss._nucleotides[n1].get_pos_base()
            r2 = ss._nucleotides[n2].get_pos_base()
            r3 = ss._nucleotides[n3].get_pos_base()
            r4 = ss._nucleotides[n4].get_pos_base()
        elif type == "back":
            r1 = ss._nucleotides[n1].get_pos_back()
            r2 = ss._nucleotides[n2].get_pos_back()
            r3 = ss._nucleotides[n3].get_pos_back()
            r4 = ss._nucleotides[n4].get_pos_back()
        else:
            base.Logger.die("origami_utils.Origami.get_local_twist: unknown type %s; use either base or back" % type)

        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        n = mid34 - mid12
        v12 = r2 - r1
        v34 = r4 - r3
        n = norm(n)
        v12 = norm(v12)
        v34 = norm(v34)
        # find the component of the base-base (or back-back) vectors that is in the plane normal to n
        v12prime = v12 - np.dot(v12, n) * n
        v34prime = v34 - np.dot(v34, n) * n
        v12prime = norm(v12prime)
        v34prime = norm(v34prime)

        
        bp_twist = np.arccos(np.dot(v12prime, v34prime))
        bp_twist *= 180./np.pi
        xp = np.cross(v12prime, v34prime)
        if np.dot(xp,n) < 0:
            # this is interesting to note: happens fairly often in edge bp (of course we should be discarding those ones anyway)
            pass #print "not aligned", bp_twist, vh_id, vvib_id, conf
        return bp_twist

    def get_radius(self, vh_id, vvib_id):
        # find radius of a base pair (i.e. to find radius of helix)
        # assumes no skip/loop
        ss = self._sys
        n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
        n2 = self.complementary_list[n1]
        # use backbone centres
        r1 = ss._nucleotides[n1].get_pos_back()
        r2 = ss._nucleotides[n2].get_pos_back()
        d = r1-r2
        d = np.sqrt(np.dot(d,d))#/2 # radius
#        d += base.EXCL_S1/2

        return d
    
    def get_rise(self, vh_id, vvib_id):
        # find rise for a base pair
        # assumes no skip/loop
        ss = self._sys
        n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
        n2 = self.complementary_list[n1]
        try:
            n3 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id+1])[0]
        except KeyError:
            base.Logger.die("origami_utils.Origami.get_rise: KeyError - probably the virtual base index argument was too high")
        n4 = self.complementary_list[n3]

        r1 = ss._nucleotides[n1].get_pos_base()
        r2 = ss._nucleotides[n2].get_pos_base()
        r3 = ss._nucleotides[n3].get_pos_base()
        r4 = ss._nucleotides[n4].get_pos_base()
        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        d = mid12 - mid34
        d = np.sqrt(np.dot(d,d))

        return d
    
    def get_weave(self, vh_id, vvib_id):
        # find weave for a pair of base pairs located at (vh_id, vvib_id) and (vh_id + 1, vvib_id)
        # assumes no skip/loop
        ss = self._sys
        n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
        n2 = self.complementary_list[n1]
        try:
            n3 = self.get_nucleotides(self.vhelix_indices[vh_id+1], self.vvib[vh_id+1][vvib_id])[0]
        except KeyError:
            base.Logger.die("origami_utils.Origami.get_weave: KeyError - probably the virtual base index argument was too high")
        n4 = self.complementary_list[n3]

        r1 = ss._nucleotides[n1].get_pos_base()
        r2 = ss._nucleotides[n2].get_pos_base()
        r3 = ss._nucleotides[n3].get_pos_base()
        r4 = ss._nucleotides[n4].get_pos_base()
        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        d = mid12 - mid34
        d = np.sqrt(np.dot(d,d))

        return d

    def get_alignment(self, trim, period):
        # find alignment of a pair of nucleotides, ie to see if ones that are supposed to line up do.
        # assumes no skip/loop / insertion/deletion
        pair_count = len(self.vvib[0]) - trim * 2 - period
        twists = range(pair_count)
        for ii in range(pair_count):
            ss = self._sys
            n1 = self.get_nucleotides(self.vhelix_indices[0], self.vvib[0][trim+ii])[0]
            n2 = self.complementary_list[n1]
            n3 = self.get_nucleotides(self.vhelix_indices[0], self.vvib[0][trim+period+ii])[0]
            n4 = self.complementary_list[n3]
            r1 = ss._nucleotides[n1].get_pos_base()
            r2 = ss._nucleotides[n2].get_pos_base()
            r3 = ss._nucleotides[n3].get_pos_base()
            r4 = ss._nucleotides[n4].get_pos_base()

            ##
            mid12 = (r1+r2)/2
            mid34 = (r3+r4)/2
            n = mid34 - mid12
            a11 = ss._nucleotides[n1]._a1
            a13 = ss._nucleotides[n3]._a1
            n = norm(n)
            a11 = norm(a11)
            a13 = norm(a13)
            # find the component of the a1 vectors that is in the plane with normal n
            a11prime = a11 - np.dot(a11, n) * n
            a13prime = a13 - np.dot(a13, n) * n
            a11prime = norm(a11prime)
            a13prime = norm(a13prime)

            twist = np.arccos(np.dot(a11prime, a13prime))
            twist *= 180./np.pi
            xp = np.cross(a11prime, a13prime)
            if np.dot(xp,n) < 0:
                twist *= -1

            twists[ii] = twist

        return twists

    def get_hj_alignment(self, start, period):
        # find alignment of a pair of nucleotides, ie to see if ones that are supposed to line up do.
        # assumes exactly 2 virtual helices
        # assumes no skip/loop / insertion/deletion
        # start should be virtual base index of the first crossover (reading left to right)
        ss = self._sys
        twists = range(2)
        for ii in (0,1):
            # when ii = 0
            #                period
            #               <------->
            # vhelix 1 ====2===  ===4=====
            #              ||       ||
            # vhelix 0 ====1========3=====
            #
            # or when ii = 1
            # vhelix 1 =====2===  ===4====
            #              ||       ||
            # vhelix 0 =====1========3====
            n1 = self.get_nucleotides(self.vhelix_indices[0], start+ii)[0]
            n2 = self.get_nucleotides(self.vhelix_indices[1], start+ii)[0]
            n3 = self.get_nucleotides(self.vhelix_indices[0], start+period+ii)[0]
            n4 = self.get_nucleotides(self.vhelix_indices[1], start+period+ii)[0]

            bbm1 = self.get_bb_midpoint(n1)
            bbm2 = self.get_bb_midpoint(n2)
            bbm3 = self.get_bb_midpoint(n3)
            bbm4 = self.get_bb_midpoint(n4)

            ##
            n = bbm3 - bbm1
            align12 = bbm2 - bbm1
            align34 = bbm4 - bbm3
            n = norm(n)
            align12 = norm(align12)
            align34 = norm(align34)
            # find the component of the align vectors that is in the plane with normal n
            align12prime = align12 - np.dot(align12, n) * n
            align34prime = align34 - np.dot(align34, n) * n
            align12prime = norm(align12prime)
            align34prime = norm(align34prime)


            twist = np.arccos(np.dot(align12prime, align34prime))
            twist *= 180./np.pi
            xp = np.cross(align12prime, align34prime)
            if np.dot(xp,n) < 0:
                twist *= -1
            twists[ii] = twist
            
        return twists

    def get_arm_vec(self, vhc,vbc,vbn):
        system = self._sys
        nuc0scid = self.get_nucleotides(vhc,vbc)[0]
        nuc0stid = self.complementary_list[nuc0scid]
        nuc1scid = self.get_nucleotides(vhc,vbn)[0]
        nuc1stid = self.complementary_list[nuc1scid]

        nuc0sc = system._nucleotides[nuc0scid]
        nuc0st = system._nucleotides[nuc0stid]
        nuc1sc = system._nucleotides[nuc1scid]
        nuc1st = system._nucleotides[nuc1stid]

        v0 = nuc0sc.get_pos_base() + nuc0st.get_pos_base()
        v1 = nuc1sc.get_pos_base() + nuc1st.get_pos_base()

        v = norm(v1-v0)
        return v

    def get_junction_normal(self, hj, armlen=3):
        # hj should come from hjs[ii] where hjs=origami.get_holliday_junctions()
        # assume double crossover, no skip/loop (not sure if matters)
        vh1 = hj[0]
        vb1 = hj[1]
        vh2 = hj[2]
        vb2 = hj[3]
        # arm vectors - in cadnano representation:
        # A====>C
        #    ||
        # D<====B

        vA = self.get_arm_vec(vh1,vb1,vb1-armlen)
        vB = self.get_arm_vec(vh2,vb2,vb2+armlen)
        vC = self.get_arm_vec(vh1,vb2,vb2+armlen)
        vD = self.get_arm_vec(vh2,vb1,vb1-armlen)

        n = np.cross(vA,vB) + np.cross(vB,vC) + np.cross(vC,vD) + np.cross(vD,vA)
        n = norm(n)

        return n

    def get_corrug_weave(self, vh1, vh2, vb, n, plane_norm):
        # find distance between adjacent helices and decompose into weave and corrugation
        # for a pair of base pairs located at (vh1, vb) and (vh2, vb)
        # assumes no skip/loop
        ss = self._sys
        n1 = self.get_nucleotides(vh1, vb)[0]
        n2 = self.complementary_list[n1]
        n3 = self.get_nucleotides(vh2, vb)[0]
        n4 = self.complementary_list[n3]

        r1 = ss._nucleotides[n1].get_pos_base()
        r2 = ss._nucleotides[n2].get_pos_base()
        r3 = ss._nucleotides[n3].get_pos_base()
        r4 = ss._nucleotides[n4].get_pos_base()
        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        d = mid12 - mid34

        # caution - wd is scalar, d and cd are vectors here
        wd = np.dot(n, d) # weave_distance, corrugation_distance
        cd = d - wd*n

        # corrugation defined negative if opposite dir to plane normal
        if np.dot(cd, plane_norm) > 0:
            fac = 1
        else:
            fac = -1
            
        
        cd = np.sqrt(np.dot(cd,cd))

        return wd, fac * cd

    def get_weave_in_plane(self, vh_id, vvib_id):
        # find weave for a pair of base pairs located at (vh_id, vvib_id) and (vh_id + 1, vvib_id)
        # assumes no skip/loop
        ss = self._sys
        n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
        n2 = self.complementary_list[n1]
        try:
            n3 = self.get_nucleotides(self.vhelix_indices[vh_id+1], self.vvib[vh_id+1][vvib_id])[0]
        except KeyError:
            base.Logger.die("origami_utils.Origami.get_weave_in_plane: KeyError - probably the virtual base index argument was too high")
        n4 = self.complementary_list[n3]

        r1 = ss._nucleotides[n1].get_pos_base()
        r2 = ss._nucleotides[n2].get_pos_base()
        r3 = ss._nucleotides[n3].get_pos_base()
        r4 = ss._nucleotides[n4].get_pos_base()
        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        d = mid12 - mid34
        d = np.sqrt(np.dot(d,d))

        return d

