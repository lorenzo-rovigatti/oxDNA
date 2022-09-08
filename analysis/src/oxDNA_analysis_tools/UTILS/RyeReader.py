from sys import stderr
import numpy as np
import pickle
from os.path import exists
from typing import List, Tuple
import os
from oxDNA_analysis_tools.UTILS.data_structures import *
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import get_chunk_size
from oxDNA_analysis_tools.UTILS.get_confs import cget_confs

####################################################################################
##########                             FILE READERS                       ##########
####################################################################################

def Chunker(file, fsize, size=1000000) -> Chunk:
    """
        Generator that yields chunks of a fixed number of bytes

        Parameters:
            file (file) : The file to read
            fsize (int) : The size of the file
            size (int) : The size of the chunks

        Returns:
            (Chunk) : The chunk
    """
    current_chunk = 0  
    while True:
        b = file.read(size)
        if not b: break
        yield Chunk(b,current_chunk*size, current_chunk * size + size > fsize, fsize)
        current_chunk+=1

def linear_read(traj_info:TrajInfo, top_info:TopInfo, chunk_size:int=None) -> List[Configuration]:
    """
        Read a trajecory without multiprocessing.  

        Produces an iterator that yields a list of <chunk_size> configurations.

        Parameters:
            traj_info (TrajInfo) : The trajectory info
            top_info (TopInfo) : The topology info
            ntopart (int) : The number of confs to read at a time

        Returns:
            list[Configuration] : list of configurations
    """
    if chunk_size is None:
        chunk_size = get_chunk_size()
    current_chunk = 0
    while True:
        print(f"INFO: processed {current_chunk*chunk_size} / {len(traj_info.idxs)} confs", end='\r', file=stderr)
        if current_chunk*chunk_size >= len(traj_info.idxs):
            break
        confs = get_confs(top_info, traj_info, current_chunk*chunk_size, chunk_size)
        yield confs
        current_chunk += 1

#calculates the length of a trajectory file
def _index(traj_file) -> List[ConfInfo]: 
    """
        Finds conf starts in a trajectory file.

        Parameters:
            traj_file (file) : The trajectory file

        Returns:
            List[ConfInfo] : The start byte and number of bytes in each conf.
    """
    def find_all(a_str, sub):
        #https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
        start = 0
        idxs = []
        while True:
            start = a_str.find(sub, start)
            if start == -1: return idxs
            idxs.append(start)
            start += len(sub) 
    # our trajectory files can be splitted by the occurence of t
    val = b"t" 
    conf_starts = []
    fsize = os.stat(traj_file).st_size 
    for chunk in Chunker(open(traj_file, 'rb'), fsize):
        idxs = np.array(find_all(chunk.block,val)) + chunk.offset # find all offsets of t in file
        conf_starts.extend(idxs)
    conf_starts = np.array(conf_starts)
    #generate a list of confs for all but the last one
    idxs = [ConfInfo(conf_starts[i], conf_starts[i+1] - conf_starts[i],i) 
                                            for i in range(len(conf_starts)-1)]
    #handle last offset
    idxs.append(ConfInfo(conf_starts[-1], fsize - conf_starts[-1], len(conf_starts)-1))
    return idxs

def get_confs(top_info:TopInfo, traj_info:TrajInfo, start_conf:int, n_confs:int) -> List[Configuration]:
    """
        Read a chunk of confs from a trajectory file.

        Parameters:
            top_info (TopInfo) : Contains the number of bases in the configuration
            traj_info (TrajInfo) : Contains metadata about the trajectory file
            start_conf (int) : The index of the first conf to read
            n_confs (int) : The number of confs to read

        Returns:
            List[Configuration] : list of configurations

    """
    indexes = traj_info.idxs
    traj_file = traj_info.path
    n_bases = top_info.nbases
    return cget_confs(indexes, traj_file, start_conf, n_confs, n_bases)

####################################################################################
##########                             FILE PARSERS                       ##########
####################################################################################

def get_top_info(top : str) -> TopInfo:
    """
        bare bones of topology info

        Parameters:
            top (str) : path to the topology file

        Returns:
            (TopInfo) : topology info
    """
    with open(top) as f:
        my_top_info = f.readline().split(' ')
        if len(my_top_info)  == 2:
            nbases = my_top_info[0]
        elif len(my_top_info) == 5:
            nbases, ndna, nres = (my_top_info[0], my_top_info[2], my_top_info[3])
        else:
            print("ERROR: malformed topology header, failed to read topology file", file=stderr)
            exit()
    return TopInfo(top, int(nbases))

def get_top_info_from_traj(traj : str) -> TopInfo:
    """
        Retrieve top and traj info without providing a topology. 

        Note that the resulting top_info will have 0 strands because that information cannot be found in the trajectory. 

        Parameters:
            traj (str) : path to the trajectory file

        Returns:
            (TopInfo, TrajInfo) : topology and trajectory info
    """
    with open(traj) as f:
        l = ''
        # dump the header
        for _ in range(4):
            l = f.readline()
        n_bases = 0
        while (l[0] != 't'):
            n_bases += 1
            l = f.readline()
            if l == '':
                break

    return TopInfo("", int(n_bases))

def get_traj_info(traj : str) -> TrajInfo:
    """
        Get the information of a trajectory file

        Parameters:
            traj (str) : path to the trajectory file

        Returns:
            (TrajInfo) : trajectory info

    """
    #if idxs is None: # handle case when we have no indexes provided
    if not(exists(traj+".pyidx")):
        idxs = _index(traj) # no index created yet
        with open(traj+".pyidx","wb") as file:
            file.write(pickle.dumps(idxs)) # save it
    else:
        #we can load the index file
        with open(traj+".pyidx","rb") as file:
            idxs = pickle.loads(file.read())
        
        #check if index file matches the trajectory, if not, regenerate.
        if idxs[-1].offset+idxs[-1].size != os.stat(traj).st_size:
            idxs = _index(traj)
            with open(traj+".pyidx","wb") as file:
                file.write(pickle.dumps(idxs))

    return TrajInfo(traj,len(idxs),idxs)

def describe(top : str, traj : str) -> Tuple[TopInfo, TrajInfo]:
    """
        retrieve top and traj info for a provided pair

        You can provide None as the topology and it will read the first conf of the traj to get the number of particles.
        Note that the TopInfo will be missing the path parameter if no topology is provided.

        Parameters:
            top (str or None) : path to the topology file
            traj (str) : path to the trajectory file

        Returns:
            (TopInfo, TrajInfo) : topology and trajectory info
    """
    if top is None:
        return (get_top_info_from_traj(traj), get_traj_info(traj))
    else:
        return (get_top_info(top), get_traj_info(traj))

def strand_describe(top) -> Tuple[System, list]:
    """
        Retrieve all information from topology file mapping nucleotides to strands.
        
        Parameters:
            top (str) : path to topology file

        Returns:
            system (System) : system object 
            monomers (list of Monomer) : list of monomers
    """
    get_neighbor = lambda x: monomers[x].id if x != -1 else None

    with open (top) as f:
        l = f.readline().split()
        nmonomers = int(l[0])
        nstrands = int(l[1])

        system = System([None] * nstrands)
        monomers = [Monomer(i, None, None, None, None, None) for i in range(nmonomers)]

        ls = f.readlines()

        l = ls[0].split()
        curr = int(l[0])
        mid = 0
        s_start = 0
        s = Strand(curr)
        monomers[mid].type = l[1]
        monomers[mid].strand = s
        monomers[mid].n3 = get_neighbor(int(l[2]))
        monomers[mid].n5 = get_neighbor(int(l[3]))
        l = ls[1].split()
        mid += 1
        while l:
            if int(l[0]) != curr:
                s.monomers = monomers[s_start:mid]
                system[curr-1] = s #this is going to do something weird with proteins
                curr = int(l[0])
                s = Strand(curr)
                s_start = mid
            
            monomers[mid].type = l[1]
            monomers[mid].strand = s
            monomers[mid].n3 = get_neighbor(int(l[2]))
            monomers[mid].n5 = get_neighbor(int(l[3]))

            mid += 1
            try:
                l = ls[mid].split()
            except IndexError:
                break  

    s.monomers = monomers[s_start:mid]
    system[curr-1] = s #this is going to do something weird with proteins

    return system, monomers

def get_input_parameter(input_file, parameter) -> str:
    """
    Gets the value of a parameter in an oxDNA input file
    
    Parameters:
        input_file (str): The path to the input file
        parameter (str): The parameter you want to get the value of

    Returns:
        value (str): The value of the parameter
    """
    fin = open(input_file)
    value = ''
    for line in fin:
        line = line.lstrip()
        if not line.startswith('#'):
            if parameter in line:
                value = line.split('=')[1].replace(' ','').replace('\n','')
    fin.close()
    if value == '':
        print("ERROR: Key {} not found in input file {}".format(parameter, input_file))
    return value

####################################################################################
##########                              CONF UTILS                        ##########
####################################################################################

def inbox(conf : Configuration, center=False) -> Configuration:
    """
        Modify the positions attribute such that all positions are inside the box.

        Parameters:
            conf (Configuration) : The configuration to inbox
            center (bool) : If True, center the configuration on the box

        Returns:
            (Configuration) : The inboxed configuration
    """
    def realMod (n, m):
        return(((n % m) + m) % m)
    def coord_in_box(p):
        p = realMod(p, conf.box)
        return(p)
    def calc_PBC_COM(conf):
        angle = (conf.positions * 2 * np.pi) / conf.box
        cm = np.array([[np.sum(np.cos(angle[:,0])), np.sum(np.sin(angle[:,0]))], 
        [np.sum(np.cos(angle[:,1])), np.sum(np.sin(angle[:,1]))], 
        [np.sum(np.cos(angle[:,2])), np.sum(np.sin(angle[:,2]))]]) / len(angle)
        return conf.box / (2 * np.pi) * (np.arctan2(-cm[:,1], -cm[:,0]) + np.pi)
    target = np.array([conf.box[0] / 2, conf.box[1] / 2, conf.box[2] / 2])
    cms = calc_PBC_COM(conf)
    positions = conf.positions + target - cms   
    new_poses = coord_in_box(positions)
    positions += (new_poses - conf.positions)
    if center:
        cms = np.mean(positions, axis=0)
        positions -= cms
    return Configuration(
        conf.time, conf.box, conf.energy,
        positions, conf.a1s, conf.a3s
    )

####################################################################################
##########                             FILE WRITERS                       ##########
####################################################################################

def write_conf(path : str, conf : Configuration, append=False) -> None:
    """
        write the conf to a file

        Parameters:
            path (str) : path to the file
            conf (Configuration) : the configuration to write
            append (bool) : if True, append to the file, if False, overwrite
    """
    out = []
    out.append('t = {}'.format(int(conf.time)))
    out.append('b = {}'.format(' '.join(conf.box.astype(str))))
    out.append('E = {}'.format(' '.join(conf.energy.astype(str))))
    for p, a1, a3 in zip(conf.positions, conf.a1s, conf.a3s):
        out.append('{} {} {} 0 0 0 0 0 0'.format(' '.join(p.astype(str)), ' '.join(a1.astype(str)), ' '.join(a3.astype(str))))
    
    mode = 'a' if append else 'w'
    with open(path,mode) as f:
        f.write("\n".join(out))

def conf_to_str(conf : Configuration) -> str:
    """
    Write configuration as a string

    Parameters:
        conf (Configuration) : The configuration to write

    Returns:
        (str) : The configuration as a string
    """
    # When writing a configuration to a file, the conversion from ndarray to string is the slowest part
    # This horrific list comp is the best solution we found
    header = f't = {int(conf.time)}\nb = {" ".join(conf.box.astype(str))}\nE = {" ".join(conf.energy.astype(str))}\n'
    return(''.join([header, ''.join([('{} {} {} 0 0 0 0 0 0\n'.format(' '.join(p.astype(str)), ' '.join(a1.astype(str)), ' '.join(a3.astype(str)))) for p, a1, a3 in zip(conf.positions, conf.a1s, conf.a3s)])]))

def get_top_string(system) -> str:
    """
        Write topology file from system object.

        Parameters:
            system (System) : system object

        Returns:
            (str) : string representation of the system in .top format

    """

    n_na = 0
    n_aa = 0
    na_strands = 0
    aa_strands = 0
    mid = -1


    header = []
    body = []

    # iterate through strands and assign sequential ids
    # this will break circular strands.
    for s in system.strands:
        #it's a nucleic acid strand
        if s.id > 0:
            na_strands += 1
            for i, m in enumerate(s.monomers):
                n_na += 1
                mid += 1
                body.append(f'{na_strands} {m.type} {-1 if i == 0 else mid-1} {-1 if i == len(s.monomers)-1 else mid+1}')

        # it's a peptide strand
        elif s.id < 0:
            aa_strands -= 1
            for i, m in enumerate(s.monomers):
                n_aa += 1
                mid += 1
                body.append(f'{aa_strands} {m.type} {-1 if i == 0 else mid-1} {-1 if i == len(s.monomers)-1 else mid+1}')

    header.append(f'{n_na+n_aa} {na_strands+aa_strands} {n_na if n_aa > 0 else ""} {n_aa if n_aa > 0 else ""}')

    out = header+body
    return '\n'.join(out)
