from sys import stderr
from multiprocessing import Pool
from typing import Callable, NamedTuple

def get_chunk_size():
    from oxDNA_analysis_tools.UTILS.chunksize import CHUNKSIZE
    return CHUNKSIZE

# set_chunk_size is in config.py

def oat_multiprocesser(nconfs:int, ncpus:int, function:Callable, callback:Callable, ctx:NamedTuple):
    """
        Runs a function on a trajectory by distributing chunks of the trajectory to each processor. Accumulates the results with a callback function.

        Parameters:
            nconfs (int): The number of configurations to process at a time
            ncpus (int): The number of processors to use
            function (function): The function to run on each chunk
            callback (function): The function to call after each chunk is processed
            ctx (NamedTuple): A named tuple containing the arguments for the function

        The callback function must use the `nonlocal` keyword to update a variable in the main thread.
    """
    chunk_size = get_chunk_size()

    pool = Pool(ncpus)

    nchunks = int(nconfs / chunk_size +
                         (1 if nconfs % chunk_size else 0))

    print(f"INFO: Processing in blocks of {chunk_size} configurations", file=stderr)
    print(f"INFO: You can modify this number by running oat config -n <number>, which will be persistent between analyses.", file=stderr)

    ## Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {nchunks} chunks")
    responses = [pool.apply_async(function,(ctx,chunk_size,i)) for i in range(nchunks)]
    print("All spawned, waiting for results")

    pool.close()

    for i,r in enumerate(responses):
        callback(i, r.get())
        print(f"finished {i+1}/{nchunks}",end="\r")

    pool.join()