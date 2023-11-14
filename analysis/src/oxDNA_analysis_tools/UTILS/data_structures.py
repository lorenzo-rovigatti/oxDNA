from __future__ import annotations
from dataclasses import dataclass
from typing import List, Dict
import inspect
import numpy as np
from typing import Union


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

        Can be accessed and modified using Python list syntax (implements __getitem__, __setitem__, and __iter__) for accessing Strand objects.

        Parameters:
            top_file (str) : The path to the topology file creating the system
            strands (List[Strand]) : A list of Strand objects
    """

    def __init__(self, top_file:str='', strands:List = []):
        self.top_file = top_file
        self.strands = strands

    def __getitem__(self, key):
        return self.strands[key]

    def __iter__(self):
        return (s for s in self.strands)

    def append(self, strand:Strand):
        """
            Append a strand to the system.

            Modifies this system in-place.

            Parameters:
                strand (Strand) : The strand to append
        """
        self.strands.append(strand)

class Strand:
    """
        Object hierarchy representation of an oxDNA strand.

        Can be accessed and modified using Python list syntax (implements __getitem__, __setitem__, and __iter__) for accessing Monomer objects.

        Parameters:
            id (int) : The id of the strand
            *initial_data (list[dict]) : Set additional attributes of the strand object (currently unused)
            **kwargs (dict) : Set addtional attributes of the strand object from key:value pairs.

        Attributes:
            __from_old (bool) : Was this created from an old-style topology file? (default : False)
            id (int) : ID of this strand in the topology file
            monomers (list[Monomer]) : List of consitutent monomer objects (default : [])
            type (str) : Type of molecule this represents (default : DNA)
            circular (bool) : Is this a circular strand? (default : False)
    """

    def __init__(self, id, *initial_data, **kwargs):
        self.__from_old = False
        self.id = id
        self.monomers = []
        self.type = 'DNA'
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
    
    # Attributes which should not be included file write out should start with '__'
    def get_kwdata(self) -> Dict[str,str]:
        """
        Returns all attributes of this object which do not begin with '__'

        Used for writing out new-style topology files.
        """
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        attributes = [a for a in attributes if not(a[0].startswith('__') or a[0].endswith('__'))]
        d = {k:str(v) for k,v in attributes if k != 'monomers'}
        return d

    def is_old(self):
        """
        Returns whether this Strand came from an old-style topology file.
        """
        return self.__from_old
    
    def set_old(self, from_old):
        """
        Sets the __from_old attribute read by `is_old`.
        """
        self.__from_old = from_old

    def is_circular(self) -> bool:
        """
        Returns the `circular` attribute.
        """
        return self.circular
    
    def get_length(self) -> int:
        """
        Returns the number of monomers in the Strand.
        """
        return len(self.monomers)
    
    def get_sequence(self) -> str:
        """
        Returns the sequence of the Strand as a string.
        """
        return ''.join([m.btype for m in self])
    
    def set_sequence(self, new_seq:str) -> None:
        """
        Modify the sequence of the Strand.  Note that the new sequence must be the same length as the old sequence.
        """
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
        A Dataclass containing information about oxDNA monomers.
    """
    id : int
    btype : str
    strand : Union[Strand, None]
    n3 : Union[int, None]
    n5 : Union[int, None]
    pair : Union[int, None]
    
