import os.path
import sys
import numpy as np
from oxDNA_analysis_tools.UTILS import base

#The Python3 successor to the trusty LorenzoReader
#Reads in oxDNA trajectory files one configuration at a time to avoid memory overflow
#_get_system() returns a System object as defined in base.py
class LorenzoReader2:
    """
    Reads oxDNA trajectory files and creates System objects from individual configurations

    The System object is a descriptive representation of DNA which is made of strands which are made of nucleotides.  It has many built-in methods defined in base.py.
    """
    def __enter__(self):
        return self

    def __exit__(self,  exc_type, exc_val, exc_tb):
        self.__del__()

    def __init__(self, configuration, topology, check_overlap=False):
        # NOTE: on paths with symlinks kills the output stream 
        # so script continues to run ... probably  better to crash 
        if not os.path.isfile(configuration):
            print("Configuration file '{}' is not readable".format(configuration))
        if not os.path.isfile(topology):
            print("Topology file '{}' is not readable".format(topology))
        self._check_overlap = check_overlap
        self._conf = open(configuration, "r")
        with open(topology, "r") as file:
            file.readline()
            self._top_lines = file.readlines()
            self.particle_count = len(self._top_lines)

    def __del__(self):
        try: 
            if self._conf: 
                self._conf.close()
        except:
            print("ERROR: The reader could not load the provided configuration/topology pair")

    def __iter__(self):
        return self

    def __next__(self):
        s = self._get_system()
        if s:
            s.inbox_system()
            return s
        else:
            raise StopIteration
  
    def _read(self, only_strand_ends=False, skip=False):
        """
        Read a single configuration from the trajectory file

        Parameters:
            <optional> only_strand_ends (bool): If True, only the first and last nucleotide of a strand will be read.  Defaults to False.
            <optional> skip (bool): If True, the configuration will be skipped. This is how you track to a specific configuration.  Defaults to False.
        
        Returns:
            system (System): The System object representing the configuration
        """
        timeline = self._conf.readline()
        time = 0
        if  len(timeline) == 0:
            return False
        else:
            time = float(timeline.split()[2])
        
        box = np.array([float(x) for x in self._conf.readline().split()[2:]])
        [E_tot, E_pot, E_kin] = [float(x) for x in self._conf.readline().split()[2:5]]

        if skip: #this could be slightly optimized by moving it up, but then it wouldn't be same length as toplines
            for tl in self._top_lines:
                self._conf.readline()
            return False

        system = base.System(box, time=time, E_pot=E_pot, E_kin=E_kin)
        base.Nucleotide.index = 0
        base.Strand.index = 0

        s = False
        strandid_current = 0
        for tl in self._top_lines:
            tls = tl.split()
            n3 = int(tls[2])
            n5 = int(tls[3])
            strandid = int(tls[0])
            if (len (tls[1]) == 1):
                try:
                    if strandid > 0:
                        b = base.base_to_number[tls[1]]
                    if strandid < 0 :
                        b = base.aa_to_number[tls[1]]
                except:
                    b = 0
                bb = b
            else:
                try:
                    tmp = int (tls[1])
                except:
                    raise ValueError ("problems in topology file with specific base pairing")

                if tmp > 0:
                    b = tmp % 4
                else:
                    b = (3 - ((3 - tmp) % 4))
                bb = tmp

            if strandid != strandid_current:
                # check for circular strand
                if n3 != -1:
                    iscircular = True
                else:
                    iscircular = False

                if s:
                    system.add_strand(s, self._check_overlap)
                if strandid >= 0:
                    s = base.Strand()
                else:
                    s = base.Peptide()
                if iscircular:
                    s.make_circular()
                strandid_current = strandid

            ls = self._conf.readline().split()
            cm = [float(x) for x in ls[0:3]]
            a1 = [float(x) for x in ls[3:6]]
            a3 = [float(x) for x in ls[6:9]]
            v = [float(x) for x in ls[9:12]]
            L = [float(x) for x in ls[12:15]]
            if not only_strand_ends or n3 == -1 or n5 == -1:
                try:
                    s.add_nucleotide(base.Nucleotide(cm, a1, a3, b, bb, v, L, n3))
                except Exception as e:
                    print("ERROR: Reader died while reading configuration with t = {}.\n\nError message:\n{}".format(time, e))
                    return False

        system.add_strand(s, self._check_overlap)

        return system


    # if only_strand_ends == True then a strand will contain only the first and the last nucleotide
    # useful for some analysis like csd for coaxial interactions
    def _get_system(self, only_strand_ends=False, N_skip=0):
        for _ in range(N_skip):
            self._read(skip=True)

        return self._read(only_strand_ends=only_strand_ends, skip=False)
