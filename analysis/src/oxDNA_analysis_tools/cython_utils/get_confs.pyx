import cython
import numpy as np
cimport numpy as numpy
from sys import stderr
from cpython.bytes cimport PyBytes_Size 
from libc.stdio cimport fopen, fclose, fread, fseek, FILE
from libc.string cimport strtok, strcpy, strlen
from libc.stdlib cimport atoi, atof, atol, malloc, free
from oxDNA_analysis_tools.UTILS.data_structures import Configuration

"""
Cython function to extract a specified number of configurations from a text trajecory file.

Parameters:
    idxs (list) : The list of starting bytes for configurations in the trajectory
    traj_path (str) : The path to the trajectory
    start (int) : The ID of the first configuration to read in idxs
    nconfs (int) : How many confs to read?
    nbases (int) : How many bases per conf?
    stride (int) : Return only every this many confs within the chunk.
    incl_vel (bool) : Are velocities included in the trajectory file?
"""
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def cget_confs(list idxs, str traj_path, int start, int nconfs, int nbases, bint incl_vel=1):
    # Number of configurations to read
    cdef int conf_count = len(idxs)
    cdef int cnconfs = nconfs
    if (start+cnconfs >= conf_count): #this handles the last chunk which may not have nconfs confs remaining.
        cnconfs = conf_count - start

    # Configuration start/size markers within the chunk
    cdef int *sizes = <int *> malloc(cnconfs * sizeof(int))
    cdef int *conf_starts = <int *> malloc(cnconfs * sizeof(int))
    if not sizes or not conf_starts:
        raise MemoryError("Could not allocate memory for the configuration sizes and starts", file=stderr)

    cdef int chunk_size = idxs[start+cnconfs-1].offset + idxs[start+cnconfs-1].size - idxs[start].offset
    for i in range(cnconfs):
        sizes[i] = idxs[start+i].size
        conf_starts[i] = idxs[start+i].offset - idxs[start].offset

    # Convert the path to something C can open
    cdef char *traj_path_c = <char *>malloc(len(traj_path)+1)
    strcpy(traj_path_c, traj_path.encode('utf-8'))
    traj_path_c[len(traj_path)] = b'\0'
    cdef FILE *traj_file = fopen(traj_path_c, "rb")
    if traj_file == NULL:
        raise OSError("Could not open trajectory file %s" % traj_path)

    fseek(traj_file, idxs[start].offset, 1)

    # Read in the current chunk
    cdef const char *chunk = <char *>malloc(chunk_size)
    fread(chunk, chunk_size, 1, traj_file)

    # Parse the chunk into Configurations
    cdef list confs = [None]*cnconfs
    for i in range(cnconfs):
        c = parse_conf(chunk, conf_starts[i], nbases, incl_vel)
        if c == 1:
            raise RuntimeError("Trajectory parsing failed on conf {} in chunk {}.  This likely means the previous conf was truncated.".format(i, start))
        else:
            confs[i] = c

    fclose(traj_file)
    free(chunk)
    free(traj_path_c)
    free(sizes)
    free(conf_starts)

    return confs

@cython.wraparound(False)
@cython.boundscheck(False)
cdef parse_conf(char *chunk, int start_byte, int nbases, bint incl_vel=1):
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
    if strlen(ptr) == 1:
        return 1

    # Get the time
    # Note that once strtok has been called, chunk is modified to have a \0 in place of t= and you can no longer get the size of chunk
    # The standard way around this is to make a copy of the target string, but we don't want to do that in case of large chunks.
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
            if not ptr:
                raise RuntimeError("Final configuration (t={}) ended earlier than expected.  It is probably truncated.".format(time))
            cposes[i*THREE+j] = atof(ptr)
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            if not ptr:
                raise RuntimeError("Final configuration (t={}) ended earlier than expected.  It is probably truncated.".format(time))
            ca1s[i*THREE+j] = atof(ptr)
        if incl_vel:
            for j in range(THREE):
                ptr = strtok(NULL, ' ')
                if not ptr:
                    raise RuntimeError("Final configuration (t={}) ended earlier than expected.  It is probably truncated.".format(time))
                ca3s[i*THREE+j] = atof(ptr)
            ptr = strtok(NULL, '\n')
        else:
            for j in range(2):
                ptr = strtok(NULL, ' ')
                if not ptr:
                    raise RuntimeError("Final configuration (t={}) ended earlier than expected.  It is probably truncated.".format(time))
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