from __future__ import annotations
from dataclasses import dataclass
from copy import deepcopy
from typing import List
from ctypes import c_ulong
import numpy as np

@dataclass
class Chunk:
    """
        Dataclass to hold information about a chunk of a trajectory
        
        Parameters:
            block (str) : The file contents of the chunk
            offset (int) : The start byte of the chunk
            is_last (bool) : Whether this is the last chunk of the trajectory
            file_size (int) : The total size of the file
    """
    block : str
    offset : int
    is_last : bool
    file_size : int

@dataclass
class ConfInfo:
    """
        Dataclass to hold metadata about a single configuration

        Parameters:
            offset (int) : The start byte of the configuration
            size (int) : The size of the configuration in bytes
            idx (int) : The index of the configuration in the trajectory
    """
    offset : int
    size : int
    id : int

@dataclass
class TrajInfo:
    """
        Dataclass to hold metadata about a trajectory file

        Parameters:
            path (str) : The path to the trajectory file
            nconfs (int) : The number of configurations in the trajectory
            idxs (List[ConfInfo]) : A list of ConfInfo objects locating configurations in the trajectory
            incl_v (bool) : Are the velocities included in the trajectory?
    """
    path : str
    nconfs : int
    idxs : List[ConfInfo]
    incl_v : bool

@dataclass
class Configuration:
    """
        Dataclass/numpy representation for oxDNA configurations

        Parameters:
            time (int) : The time step of the configuration
            box (numpy.ndarray) : The box size for the simulation
            energy (numpy.ndarray) : The potential, kinetic and total energy for the configuration
            positions (numpy.ndarray) : The positions of the nucleotides
            a1s (numpy.ndarray) : the orientations of the base pairing sites
            a3s (numpy.ndarray) : the orientations of the stacking sites
    """
    time : int
    box : np.array
    energy : np.array
    positions : np.array
    a1s : np.array
    a3s : np.array

@dataclass
class TopInfo:
    """
        Dataclass to hold metadata about a topology file

        Parameters:
            path (str) : The path to the topology file
            nbases (int) : The number of nucleotides in the topology
    """
    path : str
    nbases : int

class System:
    """
        Object hierarchy representation of an oxDNA configuration

        Parameters:
            strands (List[Strand]) : A list of Strand objects
    """
    __slots__ = ('strands')
    strands : list

    def __init__(self, strands:List = []):
        self.strands = strands

    def __getitem__(self, key):
        return self.strands[key]

    def __setitem__(self, key, value):
        self.strands[key] = value

    def __iter__(self):
        return (s for s in self.strands)

    def append(self, strand):
        """
            Append a strand to the system

            Parameters:
                strand (Strand) : The strand to append
            
            Modifies this system in-place
        """
        self.strands.append(strand)

class Strand:
    """
        Object hierarchy representation of an oxDNA strand

        Parameters:
            id (int) : The id of the strand
            monomers (List[Monomer]) : A list of Monomer objects
    """
    __slots__ = ('id', 'monomers')
    id : int
    monomers : list

    def __init__(self, id):
        self.id = id

    def __getitem__(self, key):
        return self.monomers[key]

    def __setitem__(self, key, value):
        self.monomers[key] = value

    def __iter__(self):
        return (m for m in self.monomers)

    def is_circular(self):
        return True if self.monomers[0].n3 == self.monomers[-1].id else False
    
    def get_length(self):
        return len(self.monomers)

#No, you cannot make n3, n5 and pair refs to other Monomers
#Scaffold strands are long enough that it would stack overflow while pickling for Pool processing
#Therefore you have to get the references from ids in the monomers array
@dataclass
class Monomer:
    """
        Object hierarchy representation of an oxDNA monomer

        Parameters:
            id (int) : The id of the monomer
            type (str) : The type of the monomer
            strand (Strand) : The strand the monomer belongs to
            n3 (int) : The id of the 3' neighbor of the monomer
            n5 (int) : The id of the 5' neighbor of the monomer
            pair (int) : The id of the pair of the monomer
    """
    id : int
    type : str
    strand : Strand
    n3 : int
    n5 : int
    pair : int