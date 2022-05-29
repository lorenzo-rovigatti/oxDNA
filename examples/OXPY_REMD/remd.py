#!/usr/bin/env python3

import sys
import math
import json

import numpy as np
from mpi4py import MPI
from random import random

import oxpy

try:
    from tqdm import tqdm
except ImportError:
    print("WARNING: tqdm not found, no progress meter will be shown", file=sys.stderr)
    def tqdm(iterable):
        return iterable


def get_order(i, rank, locations, mpi_nprocs):
    # figure out who is who ?
    im_rank = locations[rank]
    odd_pairs = (i % 2) == 0
    im_responsible = (im_rank % 2) == odd_pairs

    if im_responsible:
        talk_to = im_rank + 1       
    else:
        talk_to = im_rank - 1
    if talk_to >= 0 and talk_to < mpi_nprocs:
        talk_to = locations.index(talk_to)   
        return im_responsible, im_rank, talk_to 
    else:
        return None    


comm = MPI.COMM_WORLD
rank = comm.Get_rank()  # process id  
mpi_nprocs = comm.Get_size()  # number of processes running 
locations = list(range(mpi_nprocs))  # curent location of the exchange temperature
exchange = np.zeros(mpi_nprocs)  # keep track of exchanges
exchange_tries = np.zeros(mpi_nprocs)  # and attempts
rates = np.zeros(mpi_nprocs)

if rank == 0:
    history = []


def remd_log(pid, *args):
    if pid == rank:
        print("x", pid, ":", *args)


with oxpy.Context():
    input = oxpy.InputFile()  # load up conf and input
    input.init_from_filename(sys.argv[1])  # as 1st is the name of the script
    # have separate files for every running instance
    
    input["energy_file"] = f"energy_{rank}.dat"
    input["lastconf_file"] = f"last_conf_{rank}.dat"
    input["CUDA_device"] = str(rank)
    # figure out steps to run and our temperature
    pt_move_every = int(input["pt_move_every"])
    steps_to_run = int(int(input["steps"]) / pt_move_every) 
    # setup temperatures
    temperatures = input["pt_temp_list"].split(",")

    input["T"] = temperatures[rank]  # set process temperature
    input["trajectory_file"] = f"trajectory_{rank}.dat"

    # convert temperaures to ox units
    if "K" in input["T"]:
        temperatures = [oxpy.utils.Kelvin_to_oxDNA(t.replace("K", "")) for t in temperatures]   
    elif "C" in input["T"]:
        temperatures = [oxpy.utils.Celsius_to_oxDNA(t.replace("C", "")) for t in temperatures]  

    # setup simulation manager
    manager = oxpy.OxpyManager(input)
    
    # simulation loop
    for i in tqdm(range(steps_to_run)):
        manager.run(pt_move_every, True)

        # keep everyone informed where who is 
        locations = comm.bcast(locations, root=0)
        # let's wait for everyone to get started 
        swap = False
        resut = get_order(i, rank , locations, mpi_nprocs)
        if resut:
            # unpack result
            im_responsible, im_rank, talk_to = resut
            if im_responsible:
                irresp_energy = comm.recv(source=talk_to)

                irresp_temperature = temperatures[locations[talk_to]]
                # remd_log(rank,'t', i, talk_to, irresp_temperature)

                # compute my energy
                resp_energy = manager.system_energy() 
                
                fact = min(1,
                        math.exp(
                            (1 / temperatures[im_rank] - 1 / irresp_temperature)\
                                * (resp_energy - irresp_energy)
                        ))

                swap = random() < fact 
                # notify of exchange
                comm.send(swap, talk_to)
                if swap:
                    # accept configuration exchange
                    manager.update_temperature(irresp_temperature)
            else:
                
                resp_temperature = temperatures[locations[talk_to]]
                # remd_log(rank,'t', i, talk_to, resp_temperature)

                # compute the energy at resp_id T
                manager.update_temperature(resp_temperature)
                irresp_energy = manager.system_energy()        
                comm.send(irresp_energy, talk_to)
                # get decision 
                swap = comm.recv(source=talk_to) 
                if not swap:
                    manager.update_temperature(temperatures[im_rank]) 

        # handle updates to locations
        if rank == 0:
            swaped = [] 
            # handle 0
            if  resut and im_responsible:
                exchange_tries[0] += 1
                exchange_tries[talk_to] += 1 
                if swap:
                    locations[0] , locations[talk_to] = locations[talk_to] , locations[0]
                    swaped.append(0)
                    swaped.append(talk_to)    
                    exchange[0] += 1
                    exchange[talk_to] += 1
                    rates [0] = exchange[0] / exchange_tries[0]
                    rates [talk_to] = exchange[talk_to] / exchange_tries[talk_to]

            # receive from everyone but me the updates
            for j in range(1, mpi_nprocs):
                resut = comm.recv(source=j)
                if resut:
                    swap, im_responsible, im_rank, talk_to = resut
                    if im_responsible:
                        exchange_tries[j] += 1
                        exchange_tries[talk_to] += 1
                        if swap: 
                            locations[j] , locations[talk_to] = locations[talk_to] , locations[j]
                            swaped.append(j)
                            swaped.append(talk_to)
                            exchange[j] += 1
                            exchange[talk_to] += 1
                            rates [j] = exchange[0] / exchange_tries[j]
                            rates [talk_to] = exchange[talk_to] / exchange_tries[talk_to]

            history.append(locations)
            with open("history.json", "w") as file:
                file.write(
                    json.dumps(history)
                )
            remd_log(0, "rates:", rates)
        else:
            # notify 0 of update
            if resut:
                comm.send((swap, im_responsible, im_rank, talk_to), 0)
            else:
                comm.send(None, 0)  # nothing happened for this round

if rank == 0:
    print(locations)
    with open("history.json", "w") as file:
        file.write(
            json.dumps(history)
        )
