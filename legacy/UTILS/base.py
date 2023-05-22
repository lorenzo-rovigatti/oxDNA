"""
Utility functions for oxDNA.
base.py includes the classes: System, Strand, Nucleotide
    - Make initial configurations (generate.py)
    - Generate PDB, Chimera, VMD files from trajectory/config (traj2pdb.py, traj2chimera.py, traj2tcl.py)
    - Get detailed energy information (process_data/)
    - If you want to use it with oxRNA, you have to set environment variable OXRNA to 1  (export OXRNA=1) 
"""
import sys, os
try:
    import numpy as np
except:
    import mynumpy as np

def partition(s, d):
    if d in s:
        sp = s.split(d, 1)
        return sp[0], d, sp[1]
    else:
        return s, "", ""

# every defined macro in model.h must be imported in this module
def import_model_constants():
    PI = np.pi
    model = os.path.join(os.path.dirname(__file__), "../../src/model.h")
    f = open(model)
    for line in f.readlines():
        # line = line.strip().partition("//")[0].strip()
        line = (partition (line.strip (), "//")[0]).strip ()
        #macro = line.partition("#define ")[2].strip().split(" ", 1)
        macro = (partition (line, "#define ")[2]).strip().split(" ", 1)
        if len(macro) > 1:
            key, val = [x.strip() for x in macro]
            # the f needed by c to interpret numbers as floats must be removed
            # this could be a source of bugs
            val = val.replace("f", "")
            # this awful exec is needed in order to get the right results out of macro definitions
            exec "tmp = %s" % (val)
            globals()[key] = tmp
    f.close()
    
    
# if you are using the oxRNA model, then environmental variable RNA has to be set to 1. base.py then uses the constants defined in rna_model.h
def import_rna_model_constants():
    PI = np.pi
    model = os.path.join(os.path.dirname(__file__), "../../src/Interactions/rna_model.h")
    f = open(model)
    for line in f.readlines():
        # line = line.strip().partition("//")[0].strip()
        line = (partition (line.strip (), "//")[0]).strip ()
        #macro = line.partition("#define ")[2].strip().split(" ", 1)
        if  ('=' in line) and (not "float" in line ) and (not  "KEY_FOUND" in line ) and (not 'NULL' in line) and (not 'fopen' in line): 
			macro = line.strip().split('=')
				#macro = (partition (line, "=")[2]).strip().split(" ", 1)
			if len(macro) > 1:
			     key, val = [x.strip() for x in macro]
			     # the f needed by c to interpret numbers as floats must be removed
					 # this could be a source of bugs
			     val = val.replace("f", "")
			     val = val.replace(";","")
			     # this awful exec is needed in order to get the right results out of macro definitions
			     #print 'Importing ',key,val
			     exec "tmp = %s" % (val)
			     globals()[key] = tmp
    f.close()


import_model_constants()	

number_to_base = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T'}


base_to_number = {'A' : 0, 'a' : 0, 'G' : 1, 'g' : 1,
                  'C' : 2, 'c' : 2, 'T' : 3, 't' : 3,
                  'U' : 3, 'u' : 3, 'D' : 4}
                  
                  

RNA_ENV_VAR = "OXRNA"
if os.environ.get(RNA_ENV_VAR) == '1':
    RNA=True
else:
    RNA=False                  

if RNA:
	import_rna_model_constants()
	number_to_base = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'U'}
	

try:
    FLT_EPSILON = np.finfo(np.float).eps
except:
    FLT_EPSILON = 2.2204460492503131e-16

CM_CENTER_DS = POS_BASE + 0.2

OUT_TOM = 0
OUT_LORENZO = 1
OUT_VMD = 2
OUT_CREPY = 3
OUT_VMD_XYZ = 4
OUT_TEP_VMD_XYZ = 5

LENGTH_FACT = 8.518
BASE_BASE = 0.3897628551303122

RC2_BACK = EXCL_RC1**2
RC2_BASE = EXCL_RC2**2
RC2_BACK_BASE = EXCL_RC3**2

CREPY_COLOR_TABLE = ['red', 'blue', '0,0.502,0', '1,0.8,0', '0.2,0.8,1']
VMD_ELEMENT_TABLE = ['C', 'S', 'N', 'P', 'K', 'F', 'Si']

#INT_POT_TOTAL = 0
INT_HYDR = 4
INT_STACK = 2
INT_CROSS_STACK = 5
INT_COAX_STACK = 6
INT_FENE = 0
INT_EXC_BONDED = 1
INT_EXC_NONBONDED = 3


#H_CUTOFF = -0.5
H_CUTOFF = -0.1

GROOVE_ENV_VAR = "OXDNA_GROOVE"
if os.environ.get(GROOVE_ENV_VAR) == '1':
    MM_GROOVING=True
else:
    MM_GROOVING=False

# static class
class Logger(object):
    debug_level = None
    DEBUG = 0
    INFO = 1
    WARNING = 2
    CRITICAL = 3

    messages = ("DEBUG", "INFO", "WARNING", "CRITICAL")

    @staticmethod
    def log(msg, level=None, additional=None):
        if level == None: 
            level = Logger.INFO
        if level < Logger.debug_level: 
            return

        if additional != None and Logger.debug_level == Logger.DEBUG:
            print >> sys.stderr, "%s: %s (additional info: '%s')" % (Logger.messages[level], msg, additional)
        else: print >> sys.stderr, "%s: %s" % (Logger.messages[level], msg)

    @staticmethod
    def die(msg):
        Logger.log(msg, Logger.CRITICAL)
        sys.exit()


class Printable(object):
    def __init__(self):
        self._output_callables = {OUT_TOM : self._get_tom_output,
                                  OUT_LORENZO : self._get_lorenzo_output,
                                  OUT_VMD : self._get_vmd_output,
                                  OUT_CREPY : self._get_crepy_output,
                                  OUT_VMD_XYZ : self._get_vmd_xyz_output,
                                  OUT_TEP_VMD_XYZ : self._get_TEP_vmd_xyz_output
                                  }

    def get_output(self, type):
        return self._output_callables[type]()

    def _get_tom_output(self):
        raise NotImplementedError

    def _get_lorenzo_output(self):
        raise NotImplementedError

    def _get_vmd_output(self):
        raise NotImplementedError

    def _get_vmd_xyz_output(self):
        raise NotImplementedError

    def _get_TEP_vmd_xyz_output(self):
        raise NotImplementedError

    def _get_crepy_output(self):
        raise NotImplementedError


class Nucleotide(Printable):
    """
    Nucleotides compose Strands

    cm_pos --- Center of mass position
        Ex: [0, 0, 0]

    a1 --- Unit vector indicating orientation of backbone with respect to base
        Ex: [1, 0, 0]

    a3 --- Unit vector indicating orientation (tilting )of base with respect to backbone
        Ex: [0, 0, 1]

    base --- Identity of base, which must be designated with either numbers or
        letters (this is called type in the c++ code). Confusingly enough, this
        is similar to Particle.btype in oxDNA.
        
        Number: {0,1,2,3} or any number in between (-inf,-7) and (10, inf)
        To use specific sequences (or an alphabet large than four) one should
        start from the complementary pair 10 and -7. Complementary pairs are
        such that base_1 + base_2 = 3;
        
        Letter: {A,G,T,C} (these will be translated to {0, 1, 3, 2}).
        
        These are set in the dictionaries: number_to_base, base_to_number
    
    btype--- Identity of base. Unused at the moment.

    """
    index = 0

    def __init__(self, cm_pos, a1, a3, base, btype=None, L=np.array([0., 0., 0.]), v=np.array([0., 0., 0.]), n3=-1):
        Printable.__init__(self)
        self.index = Nucleotide.index
        Nucleotide.index += 1
        self.cm_pos = np.array(cm_pos)
        self._a1 = np.array (a1)
        self._a3 = np.array (a3)
        #self._a2 = np.cross (a3, a1) # implemented as property
        # _base should be a number
        if isinstance(base, int):
            pass
        else:
            try:
                base = base_to_number[base]
            except KeyError:
                Logger.log("Invalid base (%s)" % base, level=Logger.WARNING)
        self._base = base
        if btype is None:
            self._btype = base
        else:
            self._btype = btype
        self._L = L
        self._v = v
        self.n3 = n3
        self.next = -1
        self.interactions = []  #what other nucleotide this nucleotide actually interacts with
        self.init_interactions()

    def get_pos_base (self):
        """
        Returns the position of the base centroid
        Note that cm_pos is the centrod of the backbone and base.
        """
        return self.cm_pos + self._a1 * POS_BASE

    pos_base = property(get_pos_base)

    def get_pos_stack (self):
        return self.cm_pos + self._a1 * POS_STACK

    pos_stack = property (get_pos_stack)

    def get_pos_back (self):
        """
        Returns the position of the backbone centroid
        Note that cm_pos is the centrod of the backbone and base.
        """
        if os.environ.get(GROOVE_ENV_VAR) == '1':
            return self.cm_pos + self._a1 * POS_MM_BACK1 + self._a2 * POS_MM_BACK2
        elif RNA:
			return self.cm_pos + self._a1 * RNA_POS_BACK_a1 + self._a2 * RNA_POS_BACK_a2 + self._a3 * RNA_POS_BACK_a3
        else:
            return self.cm_pos + self._a1 * POS_BACK

    pos_back = property (get_pos_back)

    def get_pos_back_rel (self):
        """
        Returns the position of the backbone centroid relative to the centre of mass
        i.e. it will be a vector pointing from the c.o.m. to the backbone
        """
        return self.get_pos_back() - self.cm_pos

    def get_a2 (self):
        return np.cross (self._a3, self._a1)

    _a2 = property (get_a2)

    def copy(self, disp=None, rot=None):
        copy = Nucleotide(self.cm_pos, self._a1, self._a3, self._base, self._btype, self._L, self._v, self.n3)
        if disp is not None:
            copy.translate(disp)
        if rot is not None:
            copy.rotate(rot)

        return copy

    def translate(self, disp):
        self.cm_pos += disp
        #self.pos_base += disp
        #self.pos_stack += disp
        #self.pos_back += disp
        try: self.cm_pos_box += disp
        except: pass

    def rotate(self, R, origin=None):
        if origin is None:
            origin = self.cm_pos

        self.cm_pos = np.dot(R, self.cm_pos - origin) + origin
        self._a1 = np.dot(R, self._a1)
        self._a3 = np.dot(R, self._a3)
        # the following have been removed
        #self._a2 = np.dot(R, self._a2)
        #self.pos_base = self.cm_pos + self._a1 * POS_BASE
        #self.pos_stack = self.cm_pos + self._a1 * POS_STACK
        #self.pos_back = self.cm_pos + self._a1 * POS_BACK

    def distance (self, other, PBC=True, box=None):
        if PBC and box is None:
            if not (isinstance (box, np.ndarray) and len(box) == 3):
                Logger.die ("distance between nucleotides: if PBC is True, box must be a numpy array of length 3");
        dr = other.pos_back - self.pos_back
        if PBC:
            dr -= box * np.rint (dr / box)
        return dr

    def get_base(self):
        """
        Returns a number containing base id
        >>> v1 = np.array([0.,0.,0.])
        >>> v2 = np.array([1.,0.,0.])
        >>> v3 = np.array([0.,0.,1.])
        >>> Nucleotide(v1, v2, v3, 'A').get_base()
        'A'
        >>> Nucleotide(v1, v2, v3, "C").get_base()
        'C'
        >>> Nucleotide(v1, v2, v3, "G").get_base()
        'G'
        >>> Nucleotide(v1, v2, v3, "T").get_base()
        'T'
        >>> Nucleotide(v1, v2, v3, 1).get_base()
        'G'
        >>> Nucleotide(v1, v2, v3, 103).get_base()
        '103'
        >>> Nucleotide(v1, v2, v3, -97).get_base()
        '-97'
        """
        if type(self._base) is not int:
            try:
                b = number_to_base[self._base] # b is unused, this is just a check
            except KeyError:
                Logger.log("Nucleotide.get_base(): nucleotide %d: unknown base type '%d', defaulting to 12 (A)" % (self.index, self._base), Logger.WARNING)
        
        if self._base in [0,1,2,3]:
            return number_to_base[self._base]
        else:
            return str(self._base)

    def _get_lorenzo_output(self):
        a = np.concatenate((self.cm_pos, self._a1, self._a3, self._v, self._L))
        return " ".join(str(x) for x in a)

    def _get_crepy_output(self):
        s1 = self.cm_pos_box + self.get_pos_back_rel()
        return "%lf %lf %lf" % tuple(s1)

    def _get_ribbon_output(self):
        s1 = self.cm_pos_box + self.get_pos_back_rel()
        return "%lf %lf %lf %lf %lf %lf %lf %lf %lf" % (tuple(s1) + tuple (self._a1) + tuple (self._a2))

    def _get_tcl_detailed_output(self):
        s1 = self.cm_pos_box + self.get_pos_back_rel()
        return "%lf %lf %lf" % tuple(s1)

    def _get_tcl_output(self):
        s1 = self.cm_pos_box + self.get_pos_back_rel()
        return "%lf %lf %lf" % tuple(s1)

    def _get_vmd_xyz_output(self):
        s1 = self.cm_pos_box + self.get_pos_back_rel()
        s2 = self.cm_pos_box + POS_BASE * self._a1
        res = "C %lf %lf %lf\n" % tuple(s1)
        res += "O %lf %lf %lf\n" % (s2[0], s2[1], s2[2])

        return res
		# This prints a sphere with a line around it, hopefully.
    def _get_TEP_vmd_xyz_output(self):
				s1 = self.cm_pos_box
				s2 = self.cm_pos_box + 0.3*self._a2
				s3 = self.cm_pos_box + 0.15*self._a1
				res = "C %lf %lf %lf\n" %tuple(s1)
				res += "H %lf %lf %lf\n" %tuple(s2)
				res += "He %lf %lf %lf\n" %tuple(s3)
				
				return res

    def get_pdb_output(self, strtype, strsubid):
        s1 = self.cm_pos_box + self.get_pos_back_rel()
        s2 = self.cm_pos_box + POS_BASE*self._a1
        res = "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (2 * self.index + 1,"C",strtype,'C',self.index,' ',s1[0],s1[1],s1[2],1,7.895)
        res += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (2 * self.index + 2,"O",strtype,'C',self.index,' ',s2[0],s2[1],s2[2],1,6.316)
        return res

    def get_pdb_output_chimera(self, strtype, strsubid):
        # the s# holds the position vector of each nucleotide element
        s1 = self.cm_pos_box + self.get_pos_back_rel()
        if os.environ.get(GROOVE_ENV_VAR) == '1':
            s2 = self.cm_pos_box
            index_jump = 3
        else:
            index_jump = 2
        s3 = self.cm_pos_box + (POS_BACK + 0.68)*self._a1
        
        if RNA:
			 #s1 = self_cm_pos_box + self.get_pos_back_rel()
			 s2 = self.cm_pos_box
			 s3 = self.cm_pos_box + (POS_BASE-0.1)*self._a1
			 index_jump = 3
			 
        # some magic to get the nice ellipse for the base particle
        #s1 = self.cm_pos_box + POS_BACK* self._a1
        #s2 = self.cm_pos_box + (POS_BACK + 0.68)*self._a1
        I_b = np.matrix( ((1.1,0,0),(0,1.5,0),(0,0,2.2)) )
        R = np.matrix( (self._a1,self._a2,self._a3) )
        I_l = R.T*I_b*R
        anis = np.multiply(-1,I_l)
        anis[0,0] = (I_l[1,1]+I_l[2,2]-I_l[0,0])/2.0
        anis[1,1] = (I_l[0,0]+I_l[2,2]-I_l[1,1])/2.0
        anis[2,2] = (I_l[0,0]+I_l[1,1]-I_l[2,2])/2.0

        # print the backbone site
        res = "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (index_jump * self.index + 1,"A",strtype,'A',self.index % 10000,' ',s1[0],s1[1],s1[2],1,7.895)
        res += "ANISOU%5d %4s %3s %c%4d%c %7i%7i%7i%7i%7i%7i\n" % (index_jump * self.index + 1,"A",strtype,'A',self.index % 10000,' ',1000,1000,1000,0,0,0)
        # print the centre of mass site (if grooving is on)
        if os.environ.get(GROOVE_ENV_VAR) == '1' or RNA:
            res += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (index_jump * self.index + 2,"B",strtype,'B',self.index % 10000,' ',s2[0],s2[1],s2[2],1,7.895)
            res += "ANISOU%5d %4s %3s %c%4d%c %7i%7i%7i%7i%7i%7i\n" % (index_jump * self.index + 2,"B",strtype,'B',self.index % 10000,' ',250,250,250,0,0,0)

        if self._base == 0:
            atomtype = 'O'
        elif self._base == 1:
            atomtype = 'S'
        elif self._base == 2:
            atomtype = 'K'
        elif self._base == 3:
            atomtype = 'N'
        else:
            print >> sys.stderr, "Should not happen..."
            atomtype = 'H'

        # print the base site
        res += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (index_jump * self.index + 3, atomtype, strtype, 'C', self.index % 10000,' ',s3[0],s3[1],s3[2],1,6.316)
        res += "ANISOU%5d %4s %3s %c%4d%c %7i%7i%7i%7i%7i%7i\n" % (index_jump * self.index + 3, atomtype, strtype, 'C', self.index % 10000, ' ' , anis[0,0]*1000, anis[1,1]*1000, anis[2,2]*1000, anis[0,1]*1000, anis[0,2]*1000, anis[1,2]*1000)

        return res

    def add_H_interaction (self,nucleotide):
        self.interactions.append(nucleotide)

    def check_H_interaction(self, nucleotide):
        #print self.index, self.interactions
        if (nucleotide in self.interactions):
            return True
        else:
            return False

    def check_interaction(self,interaction_type,nucleotide):
		if(nucleotide in self.all_interactions[interaction_type].keys()):
			return True
		else:
			return False

    def get_interaction(self,nucleotide,interaction_type):
		if(nucleotide in self.all_interactions[interaction_type].keys()):
			return self.all_interactions[interaction_type][nucleotide]
		else:
			return False

    def add_interaction(self,interaction_type,nucleotide,interaction_value):
		self.all_interactions[interaction_type][nucleotide] = interaction_value

    def init_interactions(self):
		self.all_interactions = {}
		for i in range(8):
			self.all_interactions[i] = {}

    def _get_cylinder_output(self):
        # assume up to 1 interaction (h bond) per nucleotide
        if self.interactions == []:
            return np.zeros(3)
        else:
            r1 = self.get_pos_base()
            r2 = self.interactions[0]


class Strand(Printable):
    """
    Strand composed of Nucleotides
    Strands can be contained in System
    """
    index = 0

    def __init__(self):
        Printable.__init__(self)
        self.index = Strand.index
        Strand.index += 1
        self._first = -1
        self._last = -1
        self._nucleotides = []
        self._cm_pos = np.array([0., 0., 0.])
        self._cm_pos_tot = np.array([0., 0., 0.])
        self._sequence = []
        self.visible = True
        self.H_interactions = {} #shows what strands it has H-bonds with
        self._circular = False #bool on circular DNA

    def get_length(self):
        return len(self._nucleotides)

    def get_sequence(self):
        return self._sequence

    def _prepare(self, si, ni):
        self.index = si
        self._first = ni

        for n in range(self.N):
            self._nucleotides[n].index = ni + n

        self._last = ni + n
        return ni + n + 1

    def copy(self):
        copy = Strand()
        for n in self._nucleotides:
            copy.add_nucleotide(n.copy())
        return copy

    def get_cm_pos(self):
        return self._cm_pos

    def set_cm_pos(self, new_pos):
        diff = new_pos - self._cm_pos
        for n in self._nucleotides: n.translate(diff)

        self._cm_pos = new_pos

    def translate(self, amount):
        new_pos = self._cm_pos + amount
        self.set_cm_pos(new_pos)

    def rotate(self, R, origin=None):
        if origin is None:
            origin = self.cm_pos

        for n in self._nucleotides: n.rotate(R, origin)

        self._cm_pos = np.dot(R, self._cm_pos - origin) + origin

    def append (self, other):
        if not isinstance (other, Strand):
            raise ValueError

        dr = self._nucleotides[-1].distance (other._nucleotides[0], PBC=False)
        if np.sqrt(np.dot (dr, dr)) > (0.7525 + 0.25):
            print >> sys.stderr, "WARNING: Strand.append(): strands seem too far apart. Assuming you know what you are doing."

        ret = Strand()

        for n in self._nucleotides:
            ret.add_nucleotide(n)

        for n in other._nucleotides:
            ret.add_nucleotide(n)

        return ret

    def get_slice(self, start=0, end=None):
        if end is None: end = len(self._nucleotides)
        ret = Strand()
        for i in range(start, end):
            ret.add_nucleotide(self._nucleotides[i].copy())
        return ret

    def set_sequence(self, seq):
        if isinstance (seq, str):
            seq = [base_to_number[x] for x in seq]
        if len(seq) != len(self._nucleotides):
            Logger.log ("Cannot change sequence: lengths don't match", Logger.WARNING)
            return
        i = 0
        for n in self._nucleotides:
            n._base = seq[i]
            i += 1
        self._sequence = seq

    def bring_in_box_nucleotides(self, box):
        diff = np.rint(self.cm_pos / box ) * box
        for n in self._nucleotides:
            n.cm_pos_box = n.cm_pos - diff

    def add_nucleotide(self, n):
        if len(self._nucleotides) == 0:
            self._first = n.index
        n.strand = self.index
        self._nucleotides.append(n)
        self._last = n.index
        self._cm_pos_tot += n.cm_pos
        self._cm_pos = self._cm_pos_tot / self.N
        self.sequence.append(n._base)

    def overlaps_with (self, other, box):
        import energies as en

        EXC_VOL_CUTOFF = 2.0

        for n1 in self._nucleotides:
            for n2 in other._nucleotides:
                if en.excluded_volume (n1, n2, box, nearest=False) > EXC_VOL_CUTOFF:
                    #print en.excluded_volume (n1, n2, box, nearest=False)
                    return True

        return False

    def _get_lorenzo_output(self):
        if not self.visible:
            return ""

        conf = "\n".join(n.get_output(OUT_LORENZO) for n in self._nucleotides) + "\n"

        top = ""
        for n in self._nucleotides:
            if self._circular:
                if n.index == self._first:
                    n3 = self._last
                else:
                    n3 = n.index - 1
                if n.index == self._last:
                    n5 = self._first
                else:
                    n5 = n.index + 1
            else:
                if n.index == self._first:
                    n3 = -1
                else:
                    n3 = n.index - 1
                if n.index == self._last:
                    n5 = -1
                else:
                    n5 = n.index + 1
            top += "%d %s %d %d\n" % (self.index+1, n.get_base(), n3, n5)

        return conf, top

    def _get_crepy_output(self):
        if not self.visible:
            return ""
        v = [n.get_output(OUT_CREPY) for n in self._nucleotides]
        v.insert(1, "@ 0.3 C[red] DNA")
        return " ".join(v)

    def _get_ribbon_output(self):
        if not self.visible:
            return ""
        v = [n._get_ribbon_output() for n in self._nucleotides]
        v.insert(0, "0. 0. 0. @ 0.3 C[red] RIB")
        return " ".join(v)

    def get_tcl_detailed_output (self, labels=True, index=None):
        if not self.visible:
            return ""
        v = [n._get_tcl_detailed_output() for n in self._nucleotides]

        ret = ""
        # label
        # ret += "graphics 0 color black\n"
        if labels:
            if index == None:
                index = self.index
            ret += 'graphics 0 text {%lf %lf %lf}' % tuple(self._nucleotides[0].cm_pos_box)
            ret += ' "%s" size 0.75\n' % (index)

        # backbone
        for i in range (0, self.N):
            ret += "graphics 0 sphere "
            ret += "{%s} " % v[i]
            ret += "radius 0.35 resolution 20\n"
            
        # backbone-backbone
        for i, n in enumerate(self._nucleotides):
            if n.n3 != -1:
                ret += "graphics 0 cylinder "
                ret += "{%s} " % v[i]
                j = (i - 1) % (len(self._nucleotides)) 
                ret += "{%s} radius 0.25 resolution 20 filled yes\n" % v[j]

        ret += "graphics 0 sphere {%s} radius 0.35 resolution 20\n" % v[self.N-1]
        bdir = self._nucleotides[-1].cm_pos - self._nucleotides[-2].cm_pos
        bdir /= np.sqrt (np.dot (bdir, bdir))
        bdir2 = self._nucleotides[-1].cm_pos_box + bdir / 2.
        start = self._nucleotides[-1].cm_pos_box + self._nucleotides[-1].get_pos_back_rel()
        end = self._nucleotides[-1].cm_pos_box + self._nucleotides[-1].get_pos_back_rel() + bdir * 2. / 3.
        ret += "graphics 0 cone {%lf %lf %lf} {%lf %lf %lf} radius 0.35 resolution 20\n" % (start[0], start[1], start[2], end[0], end[1], end[2])

        # bases
        for n in self._nucleotides:
            ii = self._nucleotides.index(n)
            rb = n.cm_pos_box + n._a1 * POS_BASE
            if os.environ.get(GROOVE_ENV_VAR) == '1':
                rcm = n.cm_pos_box
                ret += "graphics 0 cylinder {%s} {%lf %lf %lf} radius 0.15 resolution 20 filled yes\n" % (v[ii], rcm[0], rcm[1], rcm[2])
                ret += "graphics 0 sphere {%lf %lf %lf} radius 0.15 resolution 20\n" % (rcm[0], rcm[1], rcm[2])
                ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.15 resolution 20 filled yes\n" % (rcm[0], rcm[1], rcm[2], rb[0], rb[1], rb[2])
            else:
                # backbone-base
                rB = n.cm_pos_box + n._a1 * POS_BACK
                ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.15 resolution 20 filled yes\n" % (rB[0], rB[1], rB[2], rb[0], rb[1], rb[2])
            #base 
            ret += "graphics 0 sphere {%lf %lf %lf} radius 0.20 resolution 20\n" % (rb[0], rb[1], rb[2])

        #print ret
        return ret

    def get_tcl_output (self, labels=True, index=None):
        if not self.visible:
            return ""
        v = [n._get_tcl_output() for n in self._nucleotides]
        ret = ""

        # label
        # ret += "graphics 0 color black\n"
        if labels:
            if index == None:
                index = self.index
            ret += 'graphics 0 text {%lf %lf %lf}' % tuple(self._nucleotides[0].cm_pos_box)
            ret += ' "%s" size 0.75\n' % (index)
        # backbone
        for i in range (0, self.N - 1):
            if self._nucleotides[i].n3 != -1:
                ret += "graphics 0 cylinder "
                ret += "{%s} " % v[i]
                ret += "{%s} radius 0.25 resolution 20 filled yes\n" % v[(i - 1)%self.N]
            ret += "graphics 0 sphere "
            ret += "{%s} " % v[i]
            ret += "radius 0.30 resolution 20\n"

        ret += "graphics 0 sphere {%s} radius 0.35 resolution 20\n" % v[self.N-1]
        if self._nucleotides[self.N - 1].n3 != -1:
            ret += "graphics 0 cylinder "
            ret += "{%s} " % v[i]
            ret += "{%s} radius 0.25 resolution 20 filled yes\n" % v[i + 1]

        return ret

    def get_pdb_output(self, domain=[],strand=0):
        if not self.visible:
            return ""
        strtypes = ["ALA","GLY","CYS","TYR","ARG","PHE","LYS","SER","PRO","VAL","ASN","ASP","CYX","HSP","HSD","MET","LEU"]
        strid = strtypes[( (self.index + 1) % len(strtypes) )]
        atomoutput = ""
        nid = 0

        for nucleo in self._nucleotides:
            nid += 1
            atomoutput += nucleo.get_pdb_output(strid,nid)

        return atomoutput

    def get_pdb_output_chimera(self,domain=[],strand=0):
        if not self.visible:
            return ""
        if len(domain) >0:
            strtypes = ["ALA","GLY","CYS","TYR","ARG","PHE","LYS","SER","PRO","VAL","ASN","ASP","CYX","HSP","HSD","MET","LEU"]
            nid=0
            atomoutput = ""
            for nucleo in self._nucleotides:
                try:
                    domid = domain[strand][nid]
                    if domid<0:
                        domid = len(strtypes)+domid
                except:
                    domid=0
                domname = strtypes[(domid)%len(strtypes)]
                nid += 1
                atomoutput += nucleo.get_pdb_output_chimera(domname,nid)
            return atomoutput

        else:
            strtypes = ["ALA","GLY","CYS","TYR","ARG","PHE","LYS","SER","PRO","VAL","ASN","ASP","CYX","HSP","HSD","MET","LEU"]
            strid = strtypes[( (self.index + 1) % len(strtypes) )]
            atomoutput = ""
            nid = 0

            for nucleo in self._nucleotides:
                nid += 1
                atomoutput += nucleo.get_pdb_output_chimera(strid,nid)

            return atomoutput

    def _get_vmd_xyz_output(self):
        if not self.visible:
            return ""

        return "".join(n.get_output(OUT_VMD_XYZ) for n in self._nucleotides)

    def _get_TEP_vmd_xyz_output(self):
        if not self.visible:
            return ""

        return "".join(n.get_output(OUT_TEP_VMD_XYZ) for n in self._nucleotides)

    cm_pos = property(get_cm_pos, set_cm_pos)
    N = property(get_length)
    sequence = property(get_sequence)

    def add_H_interaction(self,other_strand):
        if other_strand in self.H_interactions.keys():
            self.H_interactions[other_strand] += 1
	else:
            self.H_interactions[other_strand] = 1

    def get_H_interactions(self):
        return self.H_interactions

    def make_circular(self, check_join_len=False):
        if check_join_len:
            dr = self._nucleotides[-1].distance (self._nucleotides[0], PBC=False)
            if np.sqrt(np.dot (dr, dr)) > (0.7525 + 0.25):
                Logger.log("Strand.make_circular(): ends of the strand seem too far apart. \
                            Assuming you know what you are doing.", level=Logger.WARNING)
        self._circular = True

    def make_noncircular(self):
        self._circular = False


def parse_visibility(path):
    try:
        inp = open (path, 'r')
    except:
        Logger.log ("Visibility file `" + path + "' not found. Assuming default visibility", Logger.WARNING)
        return []

    output = []
    for linea in inp.readlines():
        linea = linea.strip().lower()
        # remove everything that comes after '#'
        linea = linea.split('#')[0]
        if len(linea) > 0: output.append(linea)

    return output



class System(object):
    """
    Object representing an oxDNA system
    Contains strands

    Arguments:
    box -- the box size of the system
        Ex: box = [50, 50, 50]

    time --- Time of the system

    E_pot --- Potential energy

    E_kin --- Kinetic energy

    """
    chimera_count = 1

    def __init__(self, box, time=0, E_pot=0, E_kin=0):
        self._time = time
        self._ready = False
        self._box = np.array(box,np.float64)
        self._N = 0
        self._N_strands = 0
        self._strands = []
        self._nucleotide_to_strand = []
        self._N_cells = np.array(np.floor (self._box / 3.), np.int)
        for kk in [0, 1, 2]:
            if self._N_cells[kk] > 100:
                self._N_cells[kk] = 100
        self._cellsides = box / self._N_cells
        self._head = [False,] * int(self._N_cells[0] * self._N_cells[1] * self._N_cells[2])
        self.E_pot = E_pot
        self.E_kin = E_kin
        self.E_tot = E_pot + E_kin
        self.cells_done = False

    def get_sequences (self):
        return [x._sequence for x in self._strands]

    _sequences = property (get_sequences)

    def get_N_Nucleotides(self):
        return self._N

    def get_N_strands(self):
        return self._N_strands

    def _prepare(self, visibility):
        sind = 0
        nind = 0
        for sind in range(self._N_strands):
            nind = self._strands[sind]._prepare(sind, nind)

        if visibility != None: self.set_visibility(visibility)

        for s in self._strands:
            s.bring_in_box_nucleotides(self._box)

    def copy (self):
		copy = System (self._box)
		for s in self._strands:
			copy.add_strand (s.copy (), check_overlap=False)
		return copy

    def get_reduced(self, according_to, bbox=None, check_overlap=False):
        visibility_list = self.get_visibility(according_to)

        if bbox is None or bbox is True:
            bbox = self._box
        elif isinstance(bbox, list) and not isinstance(bbox, np.array):
            bbox = np.array(box)
        else:
            Logger.die("Cannot reduce system, bbox not correct")

        copy = System(bbox)
        for i in range(self._N_strands):
            if visibility_list[i]:
                copy.add_strand(self._strands[i].copy(), check_overlap)

        return copy

    def join(self, other, box=None):
        if box is None:
            box = np.array([0.,0.,0.])
            for i in range(3):
                if other._box[i] > self._box[i]:
                    box[i] = other._box[i]
                else:
                    box[i] = self._box[i]

        ret = System(np.array(box, np.float64))
        for s in self._strands:
            if s.visible:
                ret.add_strand(s.copy(), check_overlap=False)
        for s in other._strands:
            if s.visible:
                ret.add_strand(s.copy(), check_overlap=False)

        return ret

    def get_visibility(self, arg=None):
        actions = {'vis' : True, 'inv' : False}
        visibility_list = [True, ] * self._N_strands

        if isinstance (arg, str):
            Logger.log ("Setting visibility with method 'file'", Logger.INFO)
            lines = parse_visibility(arg)

            for line in lines:
                # [uno, due, tre] = [p.strip() for p in line.partition("=")]
                [uno, due, tre] = [p.strip() for p in partition (line, "=")]
                if due != "=" or uno not in ["inv", "vis", "default"]:
                    Logger.log ("Lines in visibility must begin with one of inv=, vis= and default=. Skipping this line: --" + line + "--", Logger.WARNING)
                    continue

                if uno == 'default':
                    if tre not in ['inv', 'vis']:
                        Logger.log ("Wrong default in visibility file. Assuming visible as default", Logger.WARNING)
                        tre = 'vis'
                    if tre == 'inv':
                        visibility_list = [False, ] * self._N_strands
                else:
                    # filter removes all the empty strings
                    arr = [a.strip() for a in filter(None, tre.split(','))]
                    for a in arr:
                        try:
                            ind = int(a)
                        except:
                            Logger.log ("Could not cast '%s' to int. Assuming 0" % a, Logger.WARNING)
                            ind = 0
                        try:
                            visibility_list[ind] = actions[uno]
                        except:
                            Logger.log ("Strand %i does not exist in system, cannot assign visibility. Ignoring" % ind, Logger.WARNING)

        elif isinstance (arg, list):
            Logger.log("Setting visibility with method 'list'", Logger.INFO)
            # first len(arg) elements will be overwritten
            visibility_list[0:len(arg)] = arg
        else:
            if arg is not None:
                Logger.log("Argument of System.set_visibility can be a string or a list. Skipping visibility settings, assuming all strands are visible", Logger.WARNING)
            return visibility_list

        return visibility_list

    def set_visibility(self, arg=None):
        visibility_list = self.get_visibility(arg)

        for i in range(self._N_strands):
            self._strands[i].visible = visibility_list[i]

        return

    def do_cells (self):
        self._N_cells = np.array(np.floor (self._box / 3.), np.int)
        for kk in [0,1,2]:
            if self._N_cells[kk] > 100:
                self._N_cells[kk] = 100
        self._cellsides = self._box / self._N_cells
        for n in self._nucleotides:
            n.next = -1
        self._head = [False,] * int(self._N_cells[0] * self._N_cells[1] * self._N_cells[2])
        for n in self._nucleotides:
            cs = np.array((np.floor((n.cm_pos/self._box - np.rint(n.cm_pos / self._box ) + 0.5) * (1. - FLT_EPSILON) * self._box / self._cellsides)), np.int)
            cella = cs[0] + self._N_cells[0] * cs[1] + self._N_cells[0] * self._N_cells[1] * cs[2]
	    n.next = self._head[cella]
            self._head[cella] = n
        self.cells_done = True
        return

    def is_overlapping_better(self, s):
        import energies as en
        # facciamolo con le celle
        # primo ciclo sulle particelle
        if not self.cells_done:
            self.do_cells()

        for n1 in s._nucleotides:
            # trovo la cella
            cs = np.array(np.floor((n1.cm_pos/self._box - np.rint(n1.cm_pos / self._box ) + 0.5) * (1.-FLT_EPSILON) * self._box / self._cellsides), np.int)
            cella = cs[0] + self._N_cells[0] * cs[1] + self._N_cells[0] * self._N_cells[1] * cs[2]
            px = cella % self._N_cells[0]
            py = (cella / self._N_cells[0]) % self._N_cells[1]
            pz = cella / self._N_cells[0] / self._N_cells[1]
            for k1 in (-1, 0, 1):
                for k2 in (-1, 0, 1):
                    for k3 in (-1, 0, 1):
                        nx = (px + k1 + self._N_cells[0]) % self._N_cells[0]
                        ny = (py + k2 + self._N_cells[1]) % self._N_cells[1]
                        nz = (pz + k3 + self._N_cells[2]) % self._N_cells[2]
                        cella2 = nx + ny * self._N_cells[0] + nz * self._N_cells[0]* self._N_cells[1]
                        # ora scorro le particelle della seconda cella
                        n2 = self._head[cella2]
                        while n2:
                            if en.excluded_volume (n1, n2, self._box, nearest=False) > 10.: return True
                            n2 = n2.next

        return False

    def is_overlapping(self, s):
        if not self.cells_done:
            self.do_cells();
        # facciamolo con le celle
        # primo ciclo sulle particelle
        res = False
        for n1 in s._nucleotides:
            # trovo la cella
            cs = np.array(np.floor((n1.cm_pos/self._box - np.rint(n1.cm_pos / self._box ) + 0.5) * (1.-FLT_EPSILON) * self._box / self._cellsides), np.int)
            cella = cs[0] + self._N_cells[0] * cs[1] + self._N_cells[0] * self._N_cells[1] * cs[2]
            px = cella % self._N_cells[0]
            py = (cella / self._N_cells[0]) % self._N_cells[1]
            pz = cella / self._N_cells[0] / self._N_cells[1]
            for k1 in (-1, 0, 1):
                for k2 in (-1, 0, 1):
                    for k3 in (-1, 0, 1):
                        nx = (px + k1 + self._N_cells[0]) % self._N_cells[0]
                        ny = (py + k2 + self._N_cells[1]) % self._N_cells[1]
                        nz = (pz + k3 + self._N_cells[2]) % self._N_cells[2]
                        cella2 = nx + ny * self._N_cells[0] + nz * self._N_cells[0]* self._N_cells[1]
                        # ora scorro le particelle della seconda cella
                        n2 = self._head[cella2]
                        while n2:
                            dr = n1.pos_back - n2.pos_back
                            dr -= self._box * np.rint (dr / self._box)
                            if np.dot(dr, dr) < RC2_BACK: res = True

                            dr = n1.pos_base - n2.pos_base
                            dr -= self._box * np.rint (dr / self._box)
                            if np.dot(dr, dr) < RC2_BASE: res = True

                            dr = n1.pos_back - n2.pos_base
                            dr -= self._box * np.rint (dr / self._box)
                            if np.dot(dr, dr) < RC2_BACK_BASE: res = True

                            dr = n1.pos_base - n2.pos_back
                            dr -= self._box * np.rint (dr / self._box)
                            if np.dot(dr, dr) < RC2_BACK_BASE: res = True

                            n2 = n2.next
        return res

    def contains_overlaps (self):
        res = False
        N = self._N_strands
        for i in range(0, N):
            for j in range(0, i):
                if self._strands[i].overlaps_with (self._strands[j], self._box):
                    return True

        return False

    def add_strand(self, s, check_overlap=True):
        """
        Add a Strand to the System

        Returns True if non-overlapping
        Returns False if there is overlap
        """
        if check_overlap and self.is_overlapping(s):
            Nucleotide.index -= s.N
            Strand.index -= 1
            return False
        '''
        # we now make cells off-line to save time when loading
        # configurations; interactions are computed with h_bonds.py
        # most of the time anyways
        for n in s._nucleotides:
            cs = np.array((np.floor((n.cm_pos/self._box - np.rint(n.cm_pos / self._box ) + 0.5) * (1. - FLT_EPSILON) * self._box / self._cellsides)), np.int)
            cella = cs[0] + self._N_cells[0] * cs[1] + self._N_cells[0] * self._N_cells[1] * cs[2]
            n.next = self._head[cella]
            self._head[cella] = n
        '''
        self._strands.append(s)
        self._N += s.N
        self._N_strands += 1
	self.cells_done = False
        return True

    def add_strands(self, ss, check_overlap=True):
        if isinstance(ss, tuple) or isinstance(ss, list):
		added = []
		for s in ss:
                    if self.add_strand(s, check_overlap):
                        added.append(s)
		if len(added) == len(ss):
			return True
		else:
			for s in added:
				Nucleotide.index -= s.N
				Strand.index -= 1
				self._strands.pop()
				self._N -= s.N
                                self._N_strands -= 1
				self._sequences.pop()
			return False

        elif not self.add_strand(ss, check_overlap): return False

        return True

    def get_unique_seq(self):
        # we need only the unique sequences of the system
        # see http://stackoverflow.com/questions/1143379/removing-duplicates-from-list-of-lists-in-python
        unique_seq = dict((str(x), x) for x in self._sequences).values()
        return unique_seq

    def rotate (self, amount, origin=None):
        for s in self._strands:
            s.rotate (amount, origin)

    def translate (self, amount):
        for s in self._strands:
            s.translate (amount)

    def print_tcl_detailed_output (self, outname="out.tcl", visibility=None):
        self._prepare(visibility)

        try:
            f = open (outname, 'w')
        except:
            Logger.die ("tcl_output: cannot open output file. Dying")

        f.write ("color Display Background white\n")
        f.write ("mol new\n")

        box_radius = 0.1
        f.write ("graphics 0 color 0\n")
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2., -self._box[2]/2., self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2., self._box[2]/2., self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., -self._box[2]/2., -self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., self._box[2]/2., -self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., +self._box[1]/2., self._box[2]/2., self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., +self._box[1]/2., -self._box[2]/2., self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., -self._box[2]/2., self._box[0]/2., +self._box[1]/2., -self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., self._box[2]/2., self._box[0]/2., +self._box[1]/2., self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2.,-self._box[2]/2.,-self._box[0]/2., -self._box[1]/2.,+self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., -self._box[1]/2.,-self._box[2]/2.,self._box[0]/2., -self._box[1]/2.,+self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., self._box[1]/2.,-self._box[2]/2.,self._box[0]/2., self._box[1]/2.,+self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., self._box[1]/2.,-self._box[2]/2.,-self._box[0]/2., self._box[1]/2.,+self._box[2]/2.))

        # evitiamo bianco (8) perche' come sfondo
        colorid = 0
        for s in self._strands:
            f.write ("graphics 0 color %i\n" % (colorid))
            f.write (s.get_tcl_detailed_output())
            colorid = (colorid + 1)
            if colorid == 8:
                colorid += 1
            colorid = colorid % 33
        f.close ()

    def print_tcl_output (self, outname="out.tcl", visibility=None):
        self._prepare(visibility)

        try:
            f = open (outname, 'w')
        except:
            Logger.die ("tcl_output: cannot open output file. Dying")

        f.write ("color Display Background white\n")
        f.write ("mol new\n")

        box_radius = 0.1
        f.write ("graphics 0 color 0\n")
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2., -self._box[2]/2., self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2., self._box[2]/2., self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., -self._box[2]/2., -self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., self._box[2]/2., -self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., +self._box[1]/2., self._box[2]/2., self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., +self._box[1]/2., -self._box[2]/2., self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., -self._box[2]/2., self._box[0]/2., +self._box[1]/2., -self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., self._box[2]/2., self._box[0]/2., +self._box[1]/2., self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2.,-self._box[2]/2.,-self._box[0]/2., -self._box[1]/2.,+self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., -self._box[1]/2.,-self._box[2]/2.,self._box[0]/2., -self._box[1]/2.,+self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., self._box[1]/2.,-self._box[2]/2.,self._box[0]/2., self._box[1]/2.,+self._box[2]/2.))
        f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., self._box[1]/2.,-self._box[2]/2.,-self._box[0]/2., self._box[1]/2.,+self._box[2]/2.))

        # evitiamo bianco (8) perche' come sfondo
        colorid = 0
        for s in self._strands:
            f.write ("graphics 0 color %i\n" % (colorid))
            f.write (s.get_tcl_output())
            colorid = (colorid + 1)
            if colorid == 8:
                colorid += 1
            colorid = colorid % 33
        f.close ()

    def print_ribbon_output(self, name, same_colors=False, visibility=None, constr_size=None):
        self._prepare(visibility)
        if constr_size != None:
            for i in range(self.N_strands/constr_size):
                cind = constr_size * i
                s1 = self._strands[cind]
                s2 = self._strands[cind+1]
                s3 = self._strands[cind+2]
                s4 = self._strands[cind+3]

                diff1 = np.rint(s1.cm_pos / self._box ) * self._box
                diff2 = np.rint(s2.cm_pos / self._box ) * self._box
                diff3 = np.rint(s3.cm_pos / self._box ) * self._box
                diff4 = np.rint(s4.cm_pos / self._box ) * self._box

                s2.translate(np.rint((s1.cm_pos - s2.cm_pos - diff1 + diff2) / self._box) * self._box)
                s3.translate(np.rint((s1.cm_pos - s3.cm_pos - diff1 + diff3) / self._box) * self._box)
                s4.translate(np.rint((s1.cm_pos - s4.cm_pos - diff1 + diff4) / self._box) * self._box)

        unique_seq = list(self.get_unique_seq())

        if same_colors:
            n = len(unique_seq)
            colors = CREPY_COLOR_TABLE
            while len(colors) < n: colors *= 2

        f = open(name, "w")
        f.write(".Box:%lf,%lf,%lf\n" % tuple(self._box))
        for s in self._strands:
            out = s._get_ribbon_output() + "\n"
            if same_colors:
                color = colors[unique_seq.index(s.sequence)]
                out = out.replace("C[red]", "C[%s]" % color)
            f.write(out)
        f.close()

    def print_crepy_output(self, cpy_name, same_colors=False, visibility=None):
        self._prepare(visibility)

        #unique_seq = tuple(self.get_unique_seq())
        unique_seq = list(self.get_unique_seq())

        if same_colors:
            n = len(unique_seq)
            colors = CREPY_COLOR_TABLE
            while len(colors) < n: colors *= 2

        f = open(cpy_name, "w")
        f.write(".Box:%lf,%lf,%lf\n" % tuple(self._box))
        for s in self._strands:
            out = s.get_output(OUT_CREPY) + "\n"
            if same_colors:
                color = colors[unique_seq.index(s.sequence)]
                out = out.replace("C[red]", "C[%s]" % color)
            f.write(out)
        f.close()

    def print_vmd_xyz_output(self, xyz_name="out.xyz", append=False, same_colors=False, visibility=None):
        self._prepare(visibility)
        unique_seq = self.get_unique_seq()

        if same_colors:
            n = len(unique_seq)
            types = VMD_ELEMENT_TABLE
            while len(types) < n: types *= 2

        if append: flag = 'a'
        else: flag = 'w'

        # get the number of bases
        visible_nucleotides = 0
        for s in self._strands:
            if s.visible:
                visible_nucleotides += s.N

        f = open(xyz_name, flag)
        f.write("%d\n#%lf %lf %lf\n" % (2 * visible_nucleotides, self._box[0], self._box[1], self._box[2]))
        #f.write("%d\n\n" % (2*self._N, ))
        for s in self._strands:
             out = s.get_output(OUT_VMD_XYZ)
             if same_colors:
                 type = types[unique_seq.index(s.sequence)]
                 out = out.replace("C", type)
             f.write(out)
        f.close()

    def print_TEP_vmd_xyz_output(self, xyz_name="out.xyz", append=False, same_colors=False, visibility=None):
        self._prepare(visibility)
        unique_seq = self.get_unique_seq()

        if same_colors:
            n = len(unique_seq)
            types = VMD_ELEMENT_TABLE
            while len(types) < n: types *= 2

        if append: flag = 'a'
        else: flag = 'w'

        # get the number of bases
        visible_nucleotides = 0
        for s in self._strands:
            if s.visible:
                visible_nucleotides += s.N

        f = open(xyz_name, flag)
        f.write("%d\n#%lf %lf %lf\n" % (3 * visible_nucleotides, self._box[0], self._box[1], self._box[2]))
        #f.write("%d\n\n" % (2*self._N, ))
        for s in self._strands:
             out = s.get_output(OUT_TEP_VMD_XYZ)
             if same_colors:
                 type = types[unique_seq.index(s.sequence)]
                 out = out.replace("C", type)
             f.write(out)
        f.close()

    def print_lorenzo_output(self, conf_name, top_name, visibility=None):
        self._prepare(visibility)
        #print self._time, self._box[0], self._box[1], self._box[2], self.E_tot, self.E_pot, self.E_kin
        conf = "t = %lu\nb = %f %f %f\nE = %lf %lf %lf\n" % (int(self._time), self._box[0], self._box[1], self._box[2], self.E_tot, self.E_pot, self.E_kin)

        visible_strands = 0
        visible_nucleotides = 0
        for s in self._strands:
            if s.visible:
                visible_strands += 1
                visible_nucleotides += s.N

        topology = "%d %d\n" % (visible_nucleotides, visible_strands)
        for s in self._strands:
            sc, st = s.get_output(OUT_LORENZO)
            topology += st
            conf += sc

        f = open(conf_name, "w")
        f.write(conf)
        f.close()

        f = open(top_name, "w")
        f.write(topology)
        f.close()

    def print_pdb_output(self,filename,append=False, visibility=None):
        self._prepare(visibility)

        if append:
            flag = 'a'
        else:
            flag = 'w'

        result = "HEADER    frame t= " +str(self._time)+ " \nMODEL        0 \nREMARK ## 0,0\n"

        for ss in self._strands:
            result += ss.get_pdb_output()
        result += "REMARK  ######### \n\nTER \nENDMDL \n "

        f = open(filename, flag)
        f.write(result)
        f.close()

    def print_pdb_output_chimera(self,filename,append=False, visibility=None, domain=[],colour_by_seq=False):
        """
        Outputs to chimera format using custon procedure
        Note that this is the only output procedure which supports incrementing MODEL numbers
        This allowed MODELs to be read as a trajectory or set submodels when loaded in Chimera
        Also supports circular DNA

        """


        self._prepare(visibility)

        if self._N > 9999:
            Logger.log("More than 9999 nucleotides in system; recycling nucleotide indices in pdb file", Logger.WARNING)

        if append:
            flag = 'a'
            System.chimera_count += 1
        else:
            flag = 'w'
            System.chimera_count = 1

        # Increment the MODEL number to allow loading as submodels or trajectories in chimera and pymol
        #result = "HEADER    frame t= " +str(self._time)+ " \nMODEL        0 \nREMARK ## 0,0\n"
        result = "HEADER    frame t= " +str(self._time)+ " \nMODEL     " + str(System.chimera_count) + " \nREMARK ## 0,0\n"
        strand_tally=0
        for ss in self._strands:
            result += ss.get_pdb_output_chimera(domain,strand_tally)
            strand_tally +=1

        result += "REMARK  ######### \n\nTER \nENDMDL \n "
        f = open(filename, flag)
        f.write(result)
        f.close()
	
        # we now generate a command file for chimera that can be loaded into
        # chimera to make sure that all the bonds are correct and that things
        # are colored the way we want
        # http://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/framecommand.html
        chimera_colors = ["sandy brown", "blue", "red", "green", "purple","dark gray","yellow","orange","deep pink","magenta","sienna","goldenrod","gray","plum","olive drab","dark red","steel blue","sandy brown"]
        strtypes = ["ALA","GLY","CYS","TYR","ARG","PHE","LYS","SER","PRO","VAL","ASN","ASP","CYX","HSP","HSD","MET","LEU"]
        mycolors = [chimera_colors[i % len(chimera_colors)] for i in range ( len (strtypes))]
        commands = []
        # cambiamo il colore dello sfondo
        commands.append ("set bg_color white")
        # distruggiamo i legami e li ricreiamo uno per uno per essere
        # sicuri che siano giusti
        commands.append ("~bond #0")
        i = 0
        # make the bonds within each nucleotide
        for s in self._strands:
            for n in s._nucleotides:
                if os.environ.get(GROOVE_ENV_VAR) == '1' or RNA:
                    commands.append ("bond #0:%i.A:%i.B" % (i, i))
                    commands.append ("bond #0:%i.B:%i.C" % (i, i))
                else:
                    commands.append ("bond #0:%i" % (i))
                i += 1
                
        # make the bonds between nucleotide backbones
        oid = 0
        for s in self._strands:
            msid = 0
            for n in s._nucleotides:
                if n.n3 != -1:
                    commands.append ("bond #0:%i.A,%i.A" % (msid + oid, (msid-1)%(len(s._nucleotides)) + oid))
                msid += 1
            oid += len(s._nucleotides)
        # a questo punto vogliamo sistemare il nuovo comportamento su chimera-1.10.
        commands.append("preset apply int ribbons")
        commands.append("set bg_color white")
        
        # sistemiamo il colore degli strands (dopo quello delle basi)
        # DOPO aver ricreato i legami
        for i in range(len(strtypes)):
            commands.append ("color %s #0:%s" % (mycolors[i], strtypes[i]))
            #If colour_by_seq is True, colour bases separately
            if colour_by_seq==True:
                commands.append ("col cyan #0:%s@O" % (strtypes[i]))
                commands.append ("col coral #0:%s@S" % (strtypes[i]))
                commands.append ("col yellow #0:%s@K" % (strtypes[i]))
                commands.append ("col cornflower blue #0:%s@P" % (strtypes[i]))
            else:
                commands.append ("col deep sky blue #0:%s@N" % (strtypes[i]))
                commands.append ("col deep sky blue #0:%s@O" % (strtypes[i]))
                commands.append ("col deep sky blue #0:%s@S" % (strtypes[i]))
                commands.append ("col deep sky blue #0:%s@K" % (strtypes[i]))
            commands.append ("bondcolor %s #0:%s" % (mycolors[i], strtypes[i]))


        # facciamo gli ellissoidi, no? visto che ci siamo
        commands.append ("aniso scale 0.75 smoothing 4")

        # sistemiamo la dimensione dei legami
        commands.append ("setattr m stickScale 0.6 #0")

        # e per il momento ci accontentiamo
        f = open ("chimera.com", "w")
        for c in commands:
            #print c
            print >> f, c
        f.close ()

    N = property(get_N_Nucleotides)
    N_strands = property (get_N_strands)

    def get_nucleotide_list (self):
        ret = []
        for s in self._strands:
            ret += s._nucleotides
        return ret

    _nucleotides = property (get_nucleotide_list)

    def map_nucleotides_to_strands(self):
        #this function creates nucl_id -> strand_id array
        index = 0
        for i in range(len(self._strands)):
            for j in range(self._strands[i].get_length()):
                self._nucleotide_to_strand.append(i)

    def read_H_bonds(self, inputpipe):
        for line in inputpipe:
            if( len(line.split()) > 6 and line[0] != '#') :
                vals = line.split()
                nuclA = int(vals[0])
                nuclB = int(vals[1])
		#print 'ADDING',nuclA,nuclB, float(vals[INT_HYDR+2])
                self.add_H_interaction(nuclA,nuclB,float(vals[INT_HYDR+2]))

    def read_H_bonds_output_bonds(self, inputpipe):
        for line in inputpipe:
            if(line[0] != '#' and len(line.split()) > 6):
                vals = line.split()
                nuclA = int(vals[0])
                nuclB = int(vals[1])
                self.add_H_interaction(nuclA,nuclB,float(vals[3]))

    def read_all_interactions(self, inputpipe):
        for line in inputpipe:
            if(line[0] != '#' and len(line.split()) > 6):
                vals = line.split()
                nuclA = int(vals[0])
                nuclB = int(vals[1])
                for i in range(8):
                    self.add_interaction(nuclA,nuclB,i,float(vals[i+2]))
                self.add_H_interaction(nuclA,nuclB,float(vals[INT_HYDR+2]))

    def add_interaction(self,nuclA,nuclB,interaction_type,interaction_val):
            self._nucleotides[nuclA].add_interaction(interaction_type,nuclB,interaction_val)
  	    self._nucleotides[nuclB].add_interaction(interaction_type,nuclA,interaction_val)

    def get_interaction(self,nuclA,nuclB,interaction_type):
        return self._nucleotides[nuclA].get_interaction(nuclB,interaction_type)

    def add_H_interaction(self,nuclA,nuclB,interaction):
        strandA = self._nucleotide_to_strand[nuclA]
        strandB = self._nucleotide_to_strand[nuclB]
        if(interaction < H_CUTOFF):
            #print "Adding ",nuclA, " ", nuclB, " ",float(interaction)
            if strandA <= strandB:  #each interaction added just once
                self._strands[strandA].add_H_interaction(strandB)
            else:
                self._strands[strandB].add_H_interaction(standA)
	    self._nucleotides[nuclA].add_H_interaction(nuclB)
	    self._nucleotides[nuclB].add_H_interaction(nuclA)


    def show_H_interactions(self):
        sequences_with_sequences = {}
        print "# Strand1_id  Strand2_id   Number_of_H_bonds"
        for i in range(len(self._strands)):
            interactions = self._strands[i].get_H_interactions()
            for j in interactions.keys():
                print "%5d \t %5d \t %5d" % (i,j,interactions[j])

    def check_H_interaction(self,nucleotideA,nucleotideB):
        return self._nucleotides[nucleotideA].check_H_interaction(nucleotideB)

    def print_dot_bracket_output(self, filename):
        # assumes each nucleotide has at most 1 hydrogen bond, requires interactions already to be filled for nucleotide objects
        nupack_string = ""
        for n1 in range(self.get_N_Nucleotides()):
            interactions = self._nucleotides[n1].interactions
            if len(interactions) > 1:
                Logger.log ("more than 1 HB for a nucleotide", Logger.WARNING)
            if len(interactions) == 0:
                nupack_string += "."
            elif interactions[0] > n1:
                nupack_string += "("
            elif interactions[0] < n1:
                nupack_string += ")"
            else:
                Logger.log("unexpected interaction detected while building nupack string", base.Logger.CRITICAL)

        f = open(filename, "w")
        f.write(nupack_string)
        f.close()

    def print_tcl_cylinder_output(self, filename="out.tcl", show_box=True, show_labels=True, visibility=None):
        self._prepare(visibility)
        for nuc in self._nucleotides:
            nuc.printed_sphere = False
            nuc.printed_cylinder = False

        try:
            f = open(filename, "w")
        except:
            Logger.die("print_cylinder_output: cannot open file %s, dying now" % filename)

        f.write("color Display Background white\n")
        f.write ("mol new\n")

        if show_box:
            box_radius = 0.1
            f.write ("graphics 0 color 0\n")
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2., -self._box[2]/2., self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2., self._box[2]/2., self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., -self._box[2]/2., -self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., self._box[2]/2., -self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., +self._box[1]/2., self._box[2]/2., self._box[0]/2., -self._box[1]/2., self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., +self._box[1]/2., -self._box[2]/2., self._box[0]/2., -self._box[1]/2., -self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., -self._box[2]/2., self._box[0]/2., +self._box[1]/2., -self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., +self._box[1]/2., self._box[2]/2., self._box[0]/2., +self._box[1]/2., self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., -self._box[1]/2.,-self._box[2]/2.,-self._box[0]/2., -self._box[1]/2.,+self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., -self._box[1]/2.,-self._box[2]/2.,self._box[0]/2., -self._box[1]/2.,+self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (self._box[0]/2., self._box[1]/2.,-self._box[2]/2.,self._box[0]/2., self._box[1]/2.,+self._box[2]/2.))
            f.write ("graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n" % (-self._box[0]/2., self._box[1]/2.,-self._box[2]/2.,-self._box[0]/2., self._box[1]/2.,+self._box[2]/2.))

        colorid = 2
        for strand in self._strands:
            f.write ("graphics 0 color %i\n" % (colorid))
            f.write (self.get_tcl_cylinder(strand, show_labels))
        f.close ()

    def get_tcl_cylinder(self, mystrand, labels=True):
        self.map_nucleotides_to_strands()
        if not mystrand.visible:
            return ""

        # get the bb midpoints
        v = []
        for nuc1 in mystrand._nucleotides:
            if nuc1.interactions == []:
                v.append(np.zeros(3))
            else:
                nuc2 = self._nucleotides[nuc1.interactions[0]]
                r1 = nuc1.get_pos_base()
                r2 = nuc2.get_pos_base()
                v.append((r1 + r2) / 2)

        ret = ""
        # label
        if labels:
            index = mystrand.index
            ret += 'graphics 0 text {%lf %lf %lf}' % tuple(mystrand._nucleotides[0].cm_pos_box)
            ret += ' "%s" size 0.75\n' % (index)

        # cylinders
        for i in range (0, mystrand.N - 1):
            nuc1 = mystrand._nucleotides[i]
            if nuc1.printed_cylinder == False:
                nuc2 = mystrand._nucleotides[i+1]
                cylinder_type = self.tcl_cylinder_type(nuc1, nuc2, v[i], v[i+1])
                if cylinder_type == "big":
                    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 1.25 resolution 40 filled yes\n" % (v[i][0], v[i][1], v[i][2], v[i+1][0], v[i+1][1], v[i+1][2])
                    self._nucleotides[nuc2.interactions[0]].printed_cylinder = True
                elif cylinder_type == "small":
                    r1 = nuc1.cm_pos
                    r2 = nuc2.cm_pos
                    ret += "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.75 resolution 20 filled yes\n" % (r1[0], r1[1], r1[2], r2[0], r2[1], r2[2])

        # spheres
        for i in range(mystrand.N):
            nuc1 = mystrand._nucleotides[i]
            sphere_type = self.tcl_sphere_type(nuc1)
            if sphere_type == "big":
                ret += "graphics 0 sphere {%lf %lf %lf} radius 1.25 resolution 40\n" % (v[i][0], v[i][1], v[i][2])
                self._nucleotides[nuc1.interactions[0]].printed_cylinder = True
            elif sphere_type == "small":
                r1 = nuc1.cm_pos
                ret += "graphics 0 sphere {%lf %lf %lf} radius 0.75 resolution 20\n" % (r1[0], r1[1], r1[2])

        return ret

    def tcl_cylinder_type(self, nuc1, nuc2, v1, v2):
        bbm_sep = v1 - v2
        bbm_sep = np.sqrt(np.dot(bbm_sep, bbm_sep))
        if bbm_sep > 1:
            return "none"
        elif nuc1.interactions == [] or nuc2.interactions == []:
            return "small"
        else:
            return "big"

    def tcl_sphere_type(self, nuc1):
        if nuc1.interactions == []:
            return "small"
        else:
            return "big"
