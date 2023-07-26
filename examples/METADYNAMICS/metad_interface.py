import numpy as np
import pandas as pd
import os, shutil
from copy import copy
import multiprocessing as mp
import traceback
import pickle as pkl
import glob
import oxpy


class oxDNARunner(mp.Process):

    def __init__(self, index, input_file, queue):
        mp.Process.__init__(self)

        self.index = index
        self.working_dir = f"{Estimator.RUN_BASEDIR}{index}"
        self.input_file = input_file
        self.queue = queue
        
        self._pconn, self._cconn = mp.Pipe()
        self._exception = None

    def run(self):
        os.chdir(self.working_dir)
        
        with oxpy.Context():
            input_file = oxpy.InputFile()
            input_file.init_from_filename(self.input_file)
            
            self.manager = oxpy.OxpyManager(input_file)
            
            while True:
                try:
                    steps, new_potential_grid = self.queue.get()
                    # if steps is None we have to break the loop and stop the walker
                    if steps is not None:
                        # update the lookup tables of all the metadynamics-related forces
                        for force in self.manager.config_info().forces:
                            if force.group_name == "metadynamics":
                                force.potential_grid = new_potential_grid
                        
                        self.manager.run(steps)
                    else:
                        break
                except Exception as e:
                    # if an exception is raised we save the traceback and send it (together with the exception)
                    # to this process' pipe, through which we can obtain and send it to the parent's process
                    tb = traceback.format_exc()
                    self._cconn.send((e, tb))
                finally:
                    self.queue.task_done()
                
    @property
    def exception(self):
        if self._pconn.poll():
            self._exception = self._pconn.recv()
        return self._exception
            

class Estimator():
    
    INPUT_FILE = "input-meta"
    EXT_FORCES_FILE = "ext-meta"
    LAST_CONF = "output/last_conf.dat"
    BIAS_DIR = "bias"
    RUN_BASEDIR = "run-meta_"

    def __init__(self, base_dir, dX=0.1, sigma=0.2, A=0.01, dT=10, Niter=200, meta=True, tau=int(1e5),
                    N_walkers=1, use_seq_GPUs=False, p_dict={},
                    dim=1, ratio=0, angle=0, xmin=0, xmax=20, conf_interval=int(1e3),
                    save_hills=1, continue_run=0, T=None):

        self.base_dir = base_dir
        self.dX = dX
        self.sigma = sigma
        self.A = A
        self.dT = dT
        self.Niter = Niter
        self.meta = meta
        self.tau = tau
        self.N_walkers = N_walkers
        self.p_dict = p_dict
        self.dim = dim
        self.ratio = ratio
        self.angle = angle
        self.save_hills = save_hills
        self.continue_run = continue_run

        self.xmin = xmin 
        self.xmax = xmax 

        self.N_grid = int ((self.xmax - self.xmin) / self.dX) + 1
        self.dX = (self.xmax - self.xmin) / (self.N_grid - 1) 

        # we update the spacing for integer division 

        if not self.continue_run:
            if self.dim == 1:
                self.potential_grid = np.zeros(self.N_grid)

            if self.dim == 2:
                self.potential_grid = np.zeros((self.N_grid, self.N_grid))
        else:  # load the most recent file
            center_fnames = glob.glob(f"./{Estimator.BIAS_DIR}/bias_*")
            self.max_index = max([int(i.split('_')[-1]) for i in center_fnames])
            max_fname = f"./{Estimator.BIAS_DIR}/bias_{self.max_index}" 
            print(f"restarting from bias file : {max_fname}")
            self.potential_grid = pkl.load(open(max_fname, 'rb'))

        self.r_cut_int = 6 * self.sigma / self.dX

        if self.ratio:
            if self.dim == 1:
                observable_string = f'''{{
    name = pos.dat
    print_every = 100
    col_1 = {{
        type = distance
        particle_1 = {self.p_dict['p1a']}
        particle_2 = {self.p_dict['p2a']}
    }}
    col_2 = {{
        type = distance
        particle_1 = {self.p_dict['p1b']}
        particle_2 = {self.p_dict['p2b']}
    }}
}}'''
    
        elif self.angle:
            if self.dim == 1:
                observable_string = f'''{{
    name = pos.dat
    print_every = 100
    col_1 = {{
        type = distance
        particle_1 = {self.p_dict['p1a']}
        particle_2 = {self.p_dict['p2a']}
    }}
    col_2 = {{
        type = distance
        particle_1 = {self.p_dict['p3a']}
        particle_2 = {self.p_dict['p2a']}
    }}
    col_3 = {{
        type = distance
        particle_1 = {self.p_dict['p1a']}
        particle_2 = {self.p_dict['p3a']}
    }}
}}'''

        else:
            if self.dim == 1:
                observable_string = f'''{{
    name = pos.dat
    print_every = 100
    col_1 = {{
        type = distance
        particle_1 = {self.p_dict['p1a']}
        particle_2 = {self.p_dict['p2a']}
    }}
}}'''
            elif self.dim == 2:
                observable_string = f'''{{
    name = pos.dat
    print_every = 100
    col_1 = {{
        type = distance
        particle_1 = {self.p_dict['p1a']}
        particle_2 = {self.p_dict['p2a']}
    }}
    col_2 = {{
        type = distance
        particle_1 = {self.p_dict['p1b']}
        particle_2 = {self.p_dict['p2b']}
    }}
}}'''
                
        self._init_walkers()

        # parse the user-provided input file and generate the input file that will be used to run metad simulations
        with oxpy.Context(print_coda=False):
            input_file = oxpy.InputFile()
            input_file.init_from_filename(os.path.join(self.base_dir, "input"))
            
            initial_conf = os.path.join(self.base_dir, input_file["conf_file"])
            top_file = os.path.join(self.base_dir, input_file["topology"])
            
            # remove some clutter from the output
            input_file["log_file"] = "oxDNA_log.txt"
            input_file["no_stdout_energy"] = "true"
            input_file["show_overwrite_warnings"] = "false"
            input_file["print_conf_interval"] = str(conf_interval)
            # we standardise the location of the last configuration, which is also the configuration we will start from
            input_file["lastconf_file"] = Estimator.LAST_CONF
            input_file["conf_file"] = Estimator.LAST_CONF
            
            if T == None:
                self.T = oxpy.get_temperature(input_file["T"])
            else:
                self.T = oxpy.get_temperature(T)
                input_file["T"] = T
                
            
            # look for the first available output stream index
            keep_searching = True
            i = 1
            while keep_searching:
                if f"data_output_{i}" in input_file:
                    i += 1
                else:
                    keep_searching = False
            input_file[f"data_output_{i}"] = observable_string
            
            self.other_forces = ""
            # if the external forces file is non empty we save its contents
            if "external_forces_file" in input_file:
                try:
                    ext_file_path = os.path.join(self.base_dir, input_file["external_forces_file"])
                    with open(ext_file_path) as f:
                        self.other_forces = f.read()
                except:
                    pass
            
            input_file["external_forces_file"] = Estimator.EXT_FORCES_FILE
            input_file["external_forces"] = "true"
            
            if not self.continue_run:
                # delete and recreate the bias directory
                shutil.rmtree(Estimator.BIAS_DIR, ignore_errors=True)
                os.mkdir(Estimator.BIAS_DIR)
                
                # initialise the directories where simulations will run
                for w in self.walkers:
                    shutil.rmtree(w.working_dir, ignore_errors=True)
                    
                    os.mkdir(w.working_dir)
                    os.mkdir(f"{w.working_dir}/positions")
                    os.mkdir(f"{w.working_dir}/output")
                    
                    shutil.copy(initial_conf, os.path.join(w.working_dir, Estimator.LAST_CONF))
                    shutil.copy(top_file, w.working_dir)
                    
                    # write the meta input
                    with open(os.path.join(w.working_dir, Estimator.INPUT_FILE), 'w+') as f:
                        if use_seq_GPUs:
                            input_file["CUDA_device"] = str(w.index)
                        f.write(str(input_file))

        # write the initial force files                
        for w in self.walkers:
            self.write_external_forces_file(w.working_dir)
            w.start()
                
    def get_new_data(self, dir_name):
        os.system(f"tail -n 3 ./{dir_name}/pos.dat > ./{dir_name}/pos_min.dat")
        data = np.array(pd.read_csv(f'./{dir_name}/pos_min.dat', sep='\s+',
                          header=None,
                          index_col=False))

        if self.ratio:
            if self.dim == 1:
                data[:, 0] = data[:, 1] / data[:, 0] 
                data = np.arctan(data[:, 0])

        elif self.angle:
            if self.dim == 1:
                A = data[:, 0]
                B = data[:, 1]
                C = data[:, 2]
                val = (A ** 2 + B ** 2 - C ** 2) / (2 * A * B)
                val[np.where(val >= 1)] = 1  # just in case numerical precision takes us outside the domain 
                data_new = np.arccos(val)
                data = data_new
        else:
            if self.dim == 1:
                data = data.flatten()
            else:
                pass

        return data

    def save_old_trajectory(self, dir_name, index):
        os.system("cp -r ./%s/output/trajectory.dat %s/traj_files/trajectory_%s.dat" % (dir_name, dir_name, index,))

    def write_external_forces_file(self, dir_name):
        # build the initial lookup table
        grid_string = ''
        oxDNA_potential_grid = self.potential_grid * self.T
        if self.dim == 1:
            for i in oxDNA_potential_grid:
                grid_string += f"{i},"
        elif self.dim == 2: 
            for i in oxDNA_potential_grid:
                for j in i:
                    grid_string += f"{j},"
                grid_string += '|'
        
        # use our p-dictionary! 
        common = "\n".join(f"{com_name}={particles}" for com_name, particles in self.p_dict.items())
        common += f'''
group_name = metadynamics
xmin = {self.xmin}
xmax = {self.xmax}
N_grid = {self.N_grid}
potential_grid = {grid_string}
PBC = false'''
        
        force_file_name = os.path.join(dir_name, Estimator.EXT_FORCES_FILE)
        
        with open(force_file_name, "w+") as f:
            f.write(self.other_forces)
            
            if self.ratio == 1:
                if self.dim == 1:
                    for mode in range(1, 5):
                        our_string = f"{{\ntype = meta_atan_com_trap\n{common}\nmode = {mode}\n}}\n"
                        f.write(our_string)

            elif self.angle == 1:
                if self.dim == 1:
                    for mode in range(1, 4):
                        our_string = f"{{\ntype = meta_com_angle_trap\n{common}\nmode = {mode}\n}}\n"
                        f.write(our_string)

            else:
                if self.dim == 1:
                    for mode in range(1, 3):
                        our_string = f"{{\ntype = meta_com_trap\n{common}\nmode = {mode}\n}}\n"
                        f.write(our_string)

                elif self.dim == 2:
                    for mode in range(1, 5):
                        our_string = f"{{\ntype = meta_2D_com_trap\n{common}\nymin = {self.xmin}\nymax = {self.xmax}\nmode = {mode}\n}}\n"
                        f.write(our_string)

    def save_potential_grid(self, index):
        with open(f"{Estimator.BIAS_DIR}/bias_{index}", 'wb+') as f:
            pkl.dump(self.potential_grid, f)

    def save_positions(self, new_data, index, walker_index):
        with open(f"{Estimator.RUN_BASEDIR}{walker_index}/positions/all-pos", 'w+') as f:
            pkl.dump(new_data, f)

    def interpolatePotential1D(self, x, potential_grid):
        x_left = dX * np.floor(x / self.dX)
        x_right = x_left + self.dX
        ix_left = np.floor((x - self.xmin) / self.dX)
        ix_right = ix_left + 1
        f1 = potential_grid[int(ix_left)]
        f2 = potential_grid[int(ix_right)]
        fx = (x_right - x) / self.dX * f1 + (x - x_left) / self.dX * f2
        return fx

    def interpolatePotential2D(self, x, y, potential_grid):
        dX = self.dX
        xmin = self.xmin
        x_left = dX * np.floor(x / dX)
        y_left = dX * np.floor(y / dX)
        x_right = x_left + dX
        y_right = y_left + dX

        ix_left = int((x - xmin) / dX)
        iy_left = int((y - xmin) / dX)
        ix_right = ix_left + 1
        iy_right = iy_left + 1

        f11 = potential_grid[ix_left, iy_left]
        f12 = potential_grid[ix_left, iy_right]
        f21 = potential_grid[ix_right, iy_left]
        f22 = potential_grid[ix_right, iy_right]

        fxy1 = (x_right - x) / dX * f11 + (x - x_left) / dX * f21 
        fxy2 = (x_right - x) / dX * f12 + (x - x_left) / dX * f22
        local_potential = (y_right - y) / dX * fxy1 + (y - y_left) / dX * fxy2
        return local_potential

    def get_metad_potential(self, x, cx):
        U = self.A * np.exp(-((x - cx) ** 2) / (2 * self.sigma ** 2))
        return U

    def get_metad_potential2D(self, x, cx, y, cy):
        U = self.A * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * self.sigma ** 2))
        return U

    def _init_walkers(self):
        self.queue = mp.JoinableQueue()
        self.walkers = []
        for walker_index in range(self.N_walkers):
            self.walkers.append(oxDNARunner(walker_index, "input-meta", self.queue))
            
    def stop_runners(self):
        runner_args = [None, 0]
        for w in self.walkers:
            self.queue.put(runner_args)
        self.queue.join()

    def do_metadynamics_iteration(self, index):
        print("iteration %s" % (index,), flush=True)
        
        new_potential_grid = self.potential_grid * self.T
            
        # run parallel computation
        runner_args = [tau, new_potential_grid]
        for w in self.walkers: 
            self.write_external_forces_file(w.working_dir)
            self.queue.put(runner_args)
        self.queue.join()
        
        # check whether exceptions were raised by the walkers during the previous run
        for w in self.walkers:
            if w.exception is not None:
                error, traceback = w.exception
                print(f"The following error was raised during the execution of one of the child processes, aborting run:")
                print(error)
                print(traceback)
                self.stop_runners()
                exit(1)

        if index % self.save_hills == 0:
            self.save_potential_grid(index)
        '''
        for dir_name in self.dir_names:                     
            self.save_old_trajectory(dir_name,index)
        '''
       
        # load new data 
        new_collective_data = []
        for w in self.walkers:
            new_data = self.get_new_data(w.working_dir)
            # self.save_positions(new_data, index, walker_index)
            new_collective_data.append(copy(new_data))
           
        # update forces
        if self.meta:
            for data_selection in new_collective_data:
                x = data_selection[-1]
                if self.dim == 1:
                    local_potential = self.interpolatePotential1D(x, self.potential_grid)

                elif self.dim == 2:
                    local_potential = self.interpolatePotential2D(x[0], x[1], self.potential_grid)

                prefactor = np.exp(-local_potential / self.dT)

                if self.dim == 1:
                    ix_left = int((x - self.xmin) / self.dX);
                    ix_start = max(ix_left - self.r_cut_int, 0)
                    ix_end = min(ix_left + self.r_cut_int + 1, self.N_grid)
                    for ix in range(int(ix_start), int(ix_end)):
                        x_val = self.xmin + self.dX * ix
                        new_potential = self.get_metad_potential(x, x_val)
                        self.potential_grid[ix] += prefactor * new_potential

                elif self.dim == 2:
                    ix_left = int((x[0] - self.xmin) / self.dX);
                    ix_start = max(ix_left - self.r_cut_int, 0)
                    ix_end = min(ix_left + self.r_cut_int + 1, self.N_grid)

                    iy_left = int((x[1] - self.xmin) / self.dX);
                    iy_start = max(iy_left - self.r_cut_int, 0)
                    iy_end = min(iy_left + self.r_cut_int + 1, self.N_grid)

                    for ix in range(int(ix_start), int(ix_end)):
                        for iy in range(int(iy_start), int(iy_end)):
                            x_val = self.xmin + self.dX * ix
                            y_val = self.xmin + self.dX * iy
                            new_potential = self.get_metad_potential2D(x[0], x_val, x[1], y_val)
                            self.potential_grid[ix, iy] += prefactor * new_potential

    def do_run(self):
        if self.continue_run:
            for index in range(self.max_index + 1, self.max_index + 1 + self.Niter):
                self.do_metadynamics_iteration(index)
        else:
            for index in range(self.Niter):
                self.do_metadynamics_iteration(index)
                
        self.stop_runners()

        
if __name__ == '__main__': 

    import argparse
    
    parser = argparse.ArgumentParser(description="Metadynamics interface for oxDNA. Requires oxpy (oxDNA's Python bindings).")
    parser.add_argument("base_directory", help="The directory storing the base simulation files")
    parser.add_argument("--A", default=0.1, help="Initial bias-height increment")
    parser.add_argument("--sigma", default=0.05, help="Width of the deposited Gaussian")
    parser.add_argument("--tau", default=10000, help="Length of a metadynamics iteration (in time steps)")
    parser.add_argument("--dX", default=0.001, help="The spacing of the grid used to approximate the potential")
    parser.add_argument("--N_walkers", default=6, help="Number of parallel processes to be launched")
    parser.add_argument("--dT", default=3, help="Strength of the tempering (dT -> 0 leads to conventional simulations, dt -> infinity leads to conventional, non-well-tempered metadynamics)")
    parser.add_argument("--dim", default=1, help="Dimension of the order parameter")
    parser.add_argument("--p_fname", default="locs.meta", help="File storing the indexes of the particles whose coordinates are used to build the order parameters")
    parser.add_argument("--ratio", default=0, help="Use the angle defined from the ratio of the distances between centres of mass as the order parameter")
    parser.add_argument("--angle", default=0, help="Use the angle defined from three centres of mass as the order parameter")
    parser.add_argument("--Niter", default=10000, help="Number of metadynamics iterations")
    parser.add_argument("--xmin", default=0, help="The lower boundary of the potential grid")
    parser.add_argument("--xmax", default=30, help="The upper boundary of the potential grid")
    parser.add_argument("--conf_interval", default=int(1e3), help="Frequency with which configurations should be saved (in time steps)")
    parser.add_argument("--save_hills", default=1, help="Frequency with which the potential grid is sampled")
    parser.add_argument("--continue_run", default=0, help="Whether the simulation should continue or start anew")
    parser.add_argument("--T", default=None, help="The temperature at which the simulations will be run. If not set, the temperature in the initial input file will be used")
    parser.add_argument("--use_sequential_GPUs", action="store_true", help="Each walker will try to use a dedicated GPU: the first one will attempt to use GPU 1, the second one GPU 2, etc.")
    
    args = parser.parse_args()
    
    base_dir = args.base_directory
    A = float(args.A)
    sigma = float(args.sigma)
    tau = int(args.tau)
    dX = float(args.dX)
    N_walkers = int(args.N_walkers)
    dT = float(args.dT)
    dim = int(args.dim)
    p_fname = str(args.p_fname)
    ratio = bool(int(args.ratio))
    Niter = int(args.Niter)
    angle = bool(int(args.angle))
    xmin = float(args.xmin)
    xmax = float(args.xmax)
    conf_interval = int(args.conf_interval)
    save_hills = int(args.save_hills)
    continue_run = bool(int(args.continue_run))
    T = args.T
    use_seq_GPUs = args.use_sequential_GPUs
    
    # load in the pfile here.
    with open(p_fname) as f:
        p_dict = {}
        for line in f.readlines():
            com_name, particles = [l.strip() for l in line.split(':')]
            p_dict[com_name] = particles
    
    estimator = Estimator(base_dir, Niter=Niter, meta=True, dT=dT,
                    sigma=sigma, dX=dX, A=A, tau=tau,
                    N_walkers=N_walkers, use_seq_GPUs=use_seq_GPUs,
                    p_dict=p_dict, dim=dim, ratio=ratio, angle=angle, xmin=xmin, xmax=xmax,
                    conf_interval=conf_interval, save_hills=save_hills, continue_run=continue_run, T=T)
    
    estimator.do_run()
						      
