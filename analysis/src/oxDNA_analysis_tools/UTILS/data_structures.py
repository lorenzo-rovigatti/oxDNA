from __future__ import annotations
from dataclasses import dataclass
from typing import List
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
    box : np.ndarray
    energy : np.ndarray
    positions : np.ndarray
    a1s : np.ndarray
    a3s : np.ndarray

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
            top_file (str) : The path to the topology file creating the system
            strands (List[Strand]) : A list of Strand objects
    """
    __slots__ = ('top_file', 'strands')

    def __init__(self, top_file:str='', strands:List = []):
        self.top_file = top_file
        self.strands = strands

    def __getitem__(self, key):
        return self.strands[key]

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

    def __init__(self, id, *initial_data, **kwargs):
        self.id = id
        self.monomers = []
        self.type = ''
        self.circular = False
        for dictionary in initial_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __getitem__(self, key:int):
        return self.monomers[key]

    def __setitem__(self, key:int, value:Monomer):
        self.monomers[key] = value

    def __iter__(self):
        return (m for m in self.monomers)

    def is_circular(self) -> bool:
        return self.circular
    
    def get_length(self) -> int:
        return len(self.monomers)
    
    def get_sequence(self) -> str:
        return ''.join([m.btype for m in self])
    
    def set_sequence(self, new_seq:str) -> None:
        if len(new_seq) != self.get_length():
            raise RuntimeError("New sequence must be the same length as old sequence.")
        
        for m, s, in zip(self, new_seq):
            m.btype = s


#No, you cannot make n3, n5 and pair refs to other Monomers
#Scaffold strands are long enough that it would stack overflow while pickling for Pool processing
#Therefore you have to get the references from ids in the monomers array
@dataclass
class Monomer:
    """
        Object hierarchy representation of an oxDNA monomer

        Parameters:
            id (int) : The id of the monomer
            btype (str) : The type of the monomer
            strand (Strand | None) : The strand the monomer belongs to
            n3 (int) : The id of the 3' neighbor of the monomer
            n5 (int) : The id of the 5' neighbor of the monomer
            pair (int) : The id of the pair of the monomer
    """
    id : int
    btype : str
    strand : (Strand | None)
    n3 : (int | None)
    n5 : (int | None)
    pair : (int | None)