import cython
import numpy as np
cimport numpy as numpy
from sys import stderr
from cpython.bytes cimport PyBytes_Size 
from libc.stdio cimport fopen, fclose, fread, fseek, FILE
from libc.string cimport strtok, strcpy
from libc.stdlib cimport atoi, atof, atol, malloc, free
from oxDNA_analysis_tools.UTILS.data_structures import Configuration

"""
Cython function to extract a specified number of configurations from a text trajecory file.

Parameters:
    idxs (list) : The list of starting bytes for configurations in the trajectory
    traj_path (str) : The path to the trajectory
    start (int) : The ID of the first configuration to read in idxs
    nconfs (int) : How many confs to read
    nbases (int) : How many bases per conf
"""
@cython.wraparound(False)
@cython.boundscheck(False)
def cget_confs(list idxs, str traj_path, int start, int nconfs, int nbases, bint incl_vel=1):
    # Number of configurations to read
    cdef int conf_count = len(idxs)
    if (start+nconfs >= conf_count): #this handles the last chunk which may not have nconfs confs remaining.
        nconfs = conf_count - start

    # Configuration start/size markers within the chunk
    cdef int *sizes = <int *> malloc(nconfs * sizeof(int))
    cdef int *conf_starts = <int *> malloc(nconfs * sizeof(int))
    if not sizes or not conf_starts:
        raise MemoryError("ERROR: Could not allocate memory for the configuration sizes and starts", file=stderr)
    cdef int chunk_size = 0
    for i in range(nconfs):
        sizes[i] = idxs[start+i].size
        chunk_size += sizes[i]
        conf_starts[i] = idxs[start+i].offset - idxs[start].offset


    # Convert the path to something C can open
    cdef char *traj_path_c = <char *>malloc(len(traj_path)+1)
    strcpy(traj_path_c, traj_path.encode('utf-8'))
    traj_path_c[len(traj_path)] = b'\0'
    cdef FILE *traj_file = fopen(traj_path_c, "rb")
    if traj_file == NULL:
        print("ERROR: Could not open trajectory file %s" % traj_path, file=stderr)
        return
    fseek(traj_file, idxs[start].offset, 1)

    # Read in the current chunk
    cdef const char *chunk = <char *>malloc(chunk_size)
    fread(chunk, chunk_size, 1, traj_file)

    # Parse the chunk into Configurations
    cdef list confs = [None]*nconfs
    for i in range(nconfs):
        c = parse_conf(chunk, conf_starts[i], sizes[i], nbases, incl_vel)
        confs[i] = c

    fclose(traj_file)
    free(chunk)
    free(traj_path_c)
    free(sizes)
    free(conf_starts)

    return confs

@cython.wraparound(False)
@cython.boundscheck(False)
cdef parse_conf(char *chunk, int start_byte, int size, int nbases, bint incl_vel=1):
    cdef int THREE = 3
    cdef numpy.int64_t time #Windows and Unix use different precision for time. Using `long` means long trajectories can't be loaded on Windows systems.
    
    #allocate some memory for our configuration
    cdef cbox    = np.zeros(3, dtype = np.float64)
    cdef cenergy = np.zeros(3, dtype = np.float64)
    cdef cposes  = np.zeros(nbases * THREE, dtype = np.float64)
    cdef ca1s    = np.zeros(nbases * THREE, dtype = np.float64)
    cdef ca3s    = np.zeros(nbases * THREE, dtype = np.float64)

    cdef int j = 0
    cdef int i = 0

    # Get a pointer to the start of the configuration
    cdef const char *ptr = chunk + start_byte

    # Get the time
    ptr = strtok(ptr, 't =\n')
    time = np.int64(ptr)

    # Get the box
    ptr = strtok(NULL, 'b =')
    for j in range(THREE):
        cbox[j] = atof(ptr)
        ptr = strtok(NULL, ' ')
    ptr = strtok(NULL, '\nE =')
    
    cenergy[0] = atof(ptr)
    ptr = strtok(NULL, ' ')
    cenergy[1] = atof(ptr)
    ptr = strtok(NULL, ' \n')
    cenergy[2] = atof(ptr)

    # Parse the configuration itself
    for i in range(nbases):
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            cposes[i*THREE+j] = atof(ptr)
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            ca1s[i*THREE+j] = atof(ptr)
        if incl_vel:
            for j in range(THREE):
                ptr = strtok(NULL, ' ')
                ca3s[i*THREE+j] = atof(ptr)
            ptr = strtok(NULL, '\n')
        else:
            for j in range(2):
                ptr = strtok(NULL, ' ')
                ca3s[i*THREE+j] = atof(ptr)
            ptr = strtok(NULL, '\n')
            ca3s[i*THREE+2] = atof(ptr)

    return Configuration(
                time, 
                cbox, 
                cenergy, 
                cposes.reshape(nbases, THREE),
                ca1s.reshape(nbases, THREE), 
                ca3s.reshape(nbases, THREE)
            )