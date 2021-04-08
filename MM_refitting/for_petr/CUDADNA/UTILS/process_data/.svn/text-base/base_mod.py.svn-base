import sys, os
import numpy as np

# every defined macro in model.h must be imported in this module
def import_model_constants():
    PI = np.pi
    model = os.path.join(os.path.dirname(__file__), "../../model.h")
    f = open(model)
    for line in f.readlines():
        line = line.strip().partition("//")[0].strip()
        macro = line.partition("#define ")[2].strip().split(" ", 1)
        if len(macro) > 1:
            key, val = [x.strip() for x in macro]
            # the f needed by c to interpret numbers as floats must be removed
            # this could be a source of bugs
            val = val.replace("f", "")
            # this awful exec is needed in order to get the right results out of macro definitions
            exec "tmp = %s" % (val)
            globals()[key] = tmp
    f.close()

import_model_constants()

number_to_base = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T'}
    
base_to_number = {'A' : 0, 'a' : 0, 'G' : 1, 'g' : 1,
                  'C' : 2, 'c' : 2, 'T' : 3, 't' : 3}

FLT_EPSILON = np.finfo(np.float).eps

CM_CENTER_DS = POS_BASE + 0.2

OUT_TOM = 0
OUT_LORENZO = 1
OUT_VMD = 2
OUT_CREPY = 3
OUT_VMD_XYZ = 4

LENGTH_FACT = 8.518
BASE_BASE = 0.3897628551303122

RC2_BACK = EXCL_RC1**2
RC2_BASE = EXCL_RC2**2
RC2_BACK_BASE = EXCL_RC3**2

CREPY_COLOR_TABLE = ['red', 'blue', '0,0.502,0', '1,0.8,0', '0.2,0.8,1']
VMD_ELEMENT_TABLE = ['C', 'S', 'N', 'P', 'K', 'F', 'Si']


INT_POT = 0
INT_HYDR = 1
INT_STACK = 2
INT_CROSS_STACK = 3
INT_COAX_STACK = 4
INT_EXC = 5
INT_FENE = 6

H_CUTOFF = -0.1

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
        if level == None: level = Logger.INFO
        if level < Logger.debug_level: return

        if additional != None and Logger.debug_level == Logger.DEBUG:
            print >> sys.stderr, "%s: %s (additional info: '%s')" % (Logger.messages[level], msg, additional)
        else: print >> sys.stderr, "%s: %s" % (Logger.messages[level], msg)

    @staticmethod
    def die(msg):
        Logger.log(msg, Logger.CRITICAL)
        exit()
        

class Printable(object):
    def __init__(self):
        self._output_callables = {OUT_TOM : self._get_tom_output,
                                  OUT_LORENZO : self._get_lorenzo_output,
                                  OUT_VMD : self._get_vmd_output,
                                  OUT_CREPY : self._get_crepy_output,
                                  OUT_VMD_XYZ : self._get_vmd_xyz_output
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
    
    def _get_crepy_output(self):
        raise NotImplementedError


class Nucleotide(Printable):
    index = 0
    
    def __init__(self, cm_pos, a1, a3, base, L=np.array([0., 0., 0.]), v=np.array([0., 0., 0.])):
        Printable.__init__(self)

        self.index = Nucleotide.index
        Nucleotide.index += 1
        self.cm_pos = np.array(cm_pos)
        self._a1 = np.array(a1)
        self._a3 = np.array(a3)
        self._a2 = np.cross(a3, a1)
        self._base = base
        self._L = L
        self._v = v
        self.pos_base = self.cm_pos + self._a1 * POS_BASE
        self.pos_stack = self.cm_pos + self._a1 * POS_STACK
        self.pos_back = self.cm_pos + self._a1 * POS_BACK
        self.next = -1
        self.interactions = {} #what other nucleotide this nucleotide actually interacts with
        
    def copy(self, disp=None, rot=None):
        copy = Nucleotide(self.cm_pos, self._a1, self._a3, self._base, self._L, self._v)
        if disp is not None:
            copy.translate(disp)
        if rot is not None:
            copy.rotate(rot)
        
        return copy

    def translate(self, disp):
        self.cm_pos += disp
        self.pos_base += disp
        self.pos_stack += disp
        self.pos_back += disp

    def rotate(self, R, origin=None):
        if origin == None: origin = self.cm_pos

        self.cm_pos = np.dot(R, self.cm_pos - origin) + origin
        self._a1 = np.dot(R, self._a1)
        self._a3 = np.dot(R, self._a3)
        self._a2 = np.dot(R, self._a2)
        self.pos_base = self.cm_pos + self._a1 * POS_BASE
        self.pos_stack = self.cm_pos + self._a1 * POS_STACK
        self.pos_back = self.cm_pos + self._a1 * POS_BACK
        
    def get_base(self):
        try:
            b = number_to_base[self._base]
        except KeyError:
            Logger.log("nucleotide %d: unknown base type '%d', defaulting to A", Logger.WARNING)
            b = 'A'
        return b
        
    def _get_lorenzo_output(self):
        a = np.concatenate((self.cm_pos, self._a1, self._a3, self._v, self._L))
        return " ".join(str(x) for x in a)
    
    def _get_crepy_output(self):
        s1 = self.cm_pos_box + POS_BACK*self._a1
        return "%lf %lf %lf" % tuple(s1)

    def _get_vmd_xyz_output(self):
        s1 = self.cm_pos + POS_BACK*self._a1
        s2 = self.cm_pos + POS_BASE*self._a1
        res = "C %lf %lf %lf\n" % tuple(s1)
        res += "O %lf %lf %lf\n" % (s2[0], s2[1], s2[2])
        
        return res
    
    def get_pdb_output(self, strtype, strsubid):
        s1 = self.cm_pos_box + POS_BACK*self._a1
        s2 = self.cm_pos_bo + POS_BASE*self._a1
        res = "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (2 * self.index + 1,"C",strtype,'C',self.index,' ',s1[0],s1[1],s1[2],1,7.895)
        res += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (2 * self.index + 2,"O",strtype,'C',self.index,' ',s2[0],s2[1],s2[2],1,6.316)
        # Might need the following if we ever want chimera output
    	#res = "ATOM  %d %s %s C %d "  % (self.index*2+1,"C",strtype,self.index)
    	#res += "%lf %lf %lf  1.00 7.895 \n" % tuple(s1)  
    	#res += "ANISOU  %d %s %s C %d "  % (self.index*2+1,"C",strtype,self.index)
    	#res += " 1000   1000   1000   0      0      0 \n"
        #res += "ATOM  %d %s %s C %d "  % (self.index*2+2,"O",strtype,self.index)
    	#res += "%lf %lf %lf  1.00 6.316 \n" % tuple(s2)  
    	#res += "ANISOU  %d %s %s C %d "  % (self.index*2+2,"C",strtype,self.index)
    	#res += " 1000   1000   1000   0      0      0 \n"  
        return res

class Strand(Printable):
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
        
    def rotate(self, R, origin=None):
        if origin == None: origin = self.cm_pos
        
        for n in self._nucleotides: n.rotate(R, origin)
            
        self._cm_pos = np.dot(R, self._cm_pos - origin) + origin

    def get_slice(self, start=0, end=None):
        if end is None: end = len(self._nucleotides)
        ret = Strand()
        for i in range(start, start + end):
            ret.add_nucleotide(self._nucleotides[i].copy())
        return ret

    def set_sequence(self, seq):
        if len(seq) != len(self._nucleotides):
            Logger.log ("Cannot change sequence: lengths don't match", Logger.WARNING)
            return
        
        i = 0
        for n in self._nucleotides:
            n._base = seq[i]
            i += 1
    
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
        
    def _get_lorenzo_output(self):
        if not self.visible:
            return ""

        conf = "\n".join(n.get_output(OUT_LORENZO) for n in self._nucleotides) + "\n"
        
        top = ""
        for n in self._nucleotides:
            n5 = n.index +1
            if n.index == self._first:
                n3 = -1
            else:
                n3 = n.index - 1
            if n.index == self._last:
                n5 = -1
            else:
                n5 = n.index + 1
            top += "%d %c %d %d\n" % (self.index+1, n.get_base(), n3, n5)
            
        return conf, top
    
    def _get_crepy_output(self):
        if not self.visible:
            return ""
        v = [n.get_output(OUT_CREPY) for n in self._nucleotides]
        v.insert(1, "@ 0.3 C[red] DNA")
        return " ".join(v)
    
    def get_pdb_output(self):
        if not self.visible:
            return ""
        
        strtypes = ["ALA","GLY","CYS","TYR","ARG","PHE","LYS","SER","PRO","VAL","TYR","ASN","ASP"]
        strid = strtypes[( (self.index + 1) % len(strtypes) )]
        atomoutput = ""
        nid = 0
        
        for nucleo in self._nucleotides:
            nid += 1
            atomoutput += nucleo.get_pdb_output(strid,nid)
        
        return atomoutput    
	 
    def _get_vmd_xyz_output(self):
        if not self.visible:
            return ""
        
        return "".join(n.get_output(OUT_VMD_XYZ) for n in self._nucleotides)

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
    def __init__(self, box, time=0, E_pot=0, E_kin=0):
        self._time = time
        self._ready = False
        self._box = np.array(box)
        self._N = 0
        self._N_strands = 0
        self._strands = []
        self._nucleotide_to_strand = []
        self._N_cells = np.array(np.floor (self._box / 3.), dtype = int)
        self._cellsides = box / self._N_cells
        self._head = [False,] * (self._N_cells[0] * self._N_cells[1] * self._N_cells[2])
        self._sequences = []
        self.E_pot = E_pot
        self.E_kin = E_kin
        self.E_tot = E_pot + E_kin

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

    def get_reduced(self, according_to, bbox=None, check_overlaps=False):
        visibility_list = self.get_visibility(according_to)

        if bbox == None or bbox == True:
            bbox = self._box
        elif isinstance(bbox, list) and not isinstance(bbox, np.array):
            bbox = np.array(box)
        else:
            Logger.die("Cannot reduce system, bbox not correct")

        copy = System(bbox)
        for i in range(self._N_strands):
            if visibility_list[i]:
                copy.add_strand(self._strands[i].copy(), check_overlaps)

        return copy

    def join(self, other):
        box = np.array([0,0,0])
        for i in range(3):
            if other.box[i] > self.box[i]:
                box[i] = other.box[i]
            else: box[i] = self.box[i]

        ret = System(np.array(box))
        for s in self._strands:
            if s.visible:
                ret.add_strand(s.copy(), check_overlaps=False)
        for s in other._strands:
            if s.visible:
                ret.add_strand(s.copy(), check_overlaps=False)

        return ret

    def get_visibility(self, arg=None):
        actions = {'vis' : True, 'inv' : False}
        visibility_list = [True, ] * self._N_strands

        if isinstance (arg, str):
            Logger.log ("Setting visibility with method 'file'", Logger.INFO)
            lines = parse_visibility(arg)

            for line in lines:
                [uno, due, tre] = [p.strip() for p in line.partition("=")]
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

    def is_overlapping(self, s):
        # facciamolo con le celle
        # primo ciclo sulle particelle
        res = False
        for n1 in s._nucleotides:
            # trovo la cella
            cs = np.array(np.floor((n1.cm_pos/self._box - np.rint(n1.cm_pos / self._box ) + 0.5) * (1.-FLT_EPSILON) * self._box / self._cellsides), dtype=int)
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
    
    def add_strand(self, s, check_overlap=True):
        if check_overlap and self.is_overlapping(s):
            Nucleotide.index -= s.N
            Strand.index -= 1
            return False

        for n in s._nucleotides:
            cs = np.array((np.floor((n.cm_pos/self._box - np.rint(n.cm_pos / self._box ) + 0.5) * (1. - FLT_EPSILON) * self._box / self._cellsides)), dtype = int)
            cella = cs[0] + self._N_cells[0] * cs[1] + self._N_cells[0] * self._N_cells[1] * cs[2]
            n.next = self._head[cella]
            self._head[cella] = n

        self._strands.append(s)
        self._N += s.N
        self._N_strands += 1
        self._sequences.append(s.sequence)

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
        
    def print_crepy_output(self, cpy_name, same_colors=False, visibility=None):
        self._prepare(visibility)

        unique_seq = tuple(self.get_unique_seq())

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
        
    def print_lorenzo_output(self, conf_name, top_name, visibility=None):
        self._prepare(visibility)
        conf = "t = 0\nb = %lf %lf %lf\nE = 0 0 0\n" % (tuple(self._box))
        
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
        result += "REMARK  ######### \n\nTER \nENDML \n "
	f = open(filename, flag) 
	f.write(result)
        f.close()

    N = property(get_N_Nucleotides)
    N_strands = property (get_N_strands)

    def map_nucleotides_to_strands(self): 
        #this function creates nucl_id -> strand_id array
        index = 0
        for i in range(len(self._strands)):
	  for j in range(self._strands[i].get_length()):
	     self._nucleotide_to_strand.append(i)
	 
    def read_H_bonds(self, inputpipe):
      for line in inputpipe:
        if(line[0] != '#' and len(line.split()) > 6):  
           vals = line.split()
           nuclA = int(vals[0])
           nuclB = int(vals[1])
           self.add_H_interaction(nuclA,nuclB,float(vals[INT_HYDR+2]))
      
           
    def add_H_interaction(self,nuclA,nuclB,interaction):
         strandA = self._nucleotide_to_strand[nuclA]
         strandB = self._nucleotide_to_strand[nuclB]
         if(interaction < H_CUTOFF):
	   #print "Adding ",nuclA, " ", nuclB, " ",float(interaction)
           if strandA <= strandB:  #each interaction added just once
               self._strands[strandA].add_H_interaction(strandB)
           else:
               self._strands[strandB].add_H_interaction(standA)
         
    def show_H_interactions(self):
        sequences_with_sequences = {}
        print "# Strand1_id  Strand2_id   Number_of_H_bonds"
        for i in range(len(self._strands)):
	  interactions = self._strands[i].get_H_interactions()
	  for j in interactions.keys():
	    print "%5d \t %5d \t %5d" % (i,j,interactions[j]) 
	    
	    
	    
	    
