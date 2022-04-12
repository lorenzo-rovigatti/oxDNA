import numpy as np
import pandas as pd
import os
from copy import copy
import subprocess
import multiprocessing as mp
from pickle import dump
import pickle as pkl
import glob
import oxpy
import os

class oxDNARunner(mp.Process):
    def __init__(self, working_dir, input_file, conf_interval, queue):
        mp.Process.__init__(self)

        self.working_dir = working_dir
        self.input_file = input_file
        self.conf_interval = conf_interval
        self.queue = queue

    def initialise(self):
        os.chdir(self.working_dir)
            
        self.context = oxpy.Context()
        self.context.__enter__()

        input_file = oxpy.InputFile()
        input_file.init_from_filename(self.input_file)
        input_file["print_conf_interval"] = str(conf_interval)
        input_file["conf_file"] = "output/last.conf"
        
        self.manager = oxpy.OxpyManager(input_file)
        self.manager.load_options()
        self.manager.init()

    def run(self):
        self.initialise()
        
        while True:
            steps = self.queue.get()
            if steps is None:
                break
            self.manager.run(steps)
            self.queue.task_done()
            

class Estimator():
    def __init__(self,dX=0.1,sigma=0.2,A=0.01,dT=10,Niter=200,meta=True,wt=False,tau=int(1e5),
                    N_walkers=1,
                    p_dict = {},
                    dim = 1, ratio = 0,angle = 0,xmin=0,xmax = 20,conf_interval = int(1e3),
                    save_hills=1,continue_run = 0,T=295):

        print ("ANGLE:",angle)

        self.dX = dX
        self.sigma = sigma
        self.A = A
        self.dT = dT
        self.T = T
        self.Niter = Niter
        self.meta = meta
        self.wt = wt
        self.tau = tau
        self.N_walkers = N_walkers
        self.p_dict = p_dict
        self.dim = dim
        self.ratio = ratio
        self.angle = angle
        self.conf_interval = conf_interval
        self.save_hills = save_hills
        self.continue_run = continue_run

        self.dir_names = [f"run-meta_{index}" for index in range(N_walkers)]

        self.xmin = xmin 
        self.xmax = xmax 

        self.N_grid = int ((self.xmax - self.xmin)/self.dX) + 1
        self.dX = (self.xmax - self.xmin) / (self.N_grid - 1) 

        # we update the spacing for integer division 

        if not self.continue_run:
            if self.dim == 1:
                self.potential_grid = np.zeros(self.N_grid)

            if self.dim == 2:
                self.potential_grid = np.zeros((self.N_grid,self.N_grid))
        else: #load the most recent file
            center_fnames = glob.glob('./center_data/centers_*')
            self.max_index = max([int(i.split('_')[-1]) for i in center_fnames])
            max_fname = f"./center_data/centers_{self.max_index}" 
            print(f"restarting from centers file : {max_fname}")
            self.potential_grid = pkl.load(open(max_fname,'rb'))

        self.oxdivkT = 4.142e-20 /  (self.T*1.38e-23) 

        self.r_cut_int = 6*self.sigma / self.dX

        if self.ratio:
            if self.dim == 1:
                observable_string = f'''data_output_1 = {{
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
                observable_string = f'''data_output_1 = {{
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
                observable_string = f'''data_output_1 = {{
    name = pos.dat
    print_every = 100
    col_1 = {{
        type = distance
        particle_1 = {self.p_dict['p1a']}
        particle_2 = {self.p_dict['p2a']}
    }}
}}'''
            elif self.dim == 2:
                observable_string = f'''data_output_1 = {{
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

        # write from input to a meta input 
        with open('run-meta/input','r') as f:
            input_string = f.read()
        with open('run-meta/input-meta','w+') as f:
            f.write(input_string) 
            f.write(observable_string)

    def get_new_data(self,dir_name):
        os.system(f"tail -n 3 ./{dir_name}/pos.dat > ./{dir_name}/pos_min.dat")
        data = np.array(pd.read_csv(f'./{dir_name}/pos_min.dat',sep = '\s+',
                          header = None,
                          index_col = False))

        if self.ratio :
            if self.dim == 1:
                data[:,0] = data[:,1]/data[:,0] 
                data = np.arctan(data[:,0])
            elif self.dim == 2:
                data[:,0] = data[:,1]/data[:,0]
                data[:,1] = data[:,3]/data[:,2]
                data = np.arctan(data[:,:2])

        elif self.angle:
            #we're not even getting here?
            if self.dim == 1:
                A = data[:,0]
                B = data[:,1]
                C = data[:,2]
                val = (A**2+B**2-C**2)/(2*A*B)
                val[np.where(val >=1)] = 1 # just in case numerical precision takes us outside the domain 
                data_new = np.arccos(val)
                data = data_new
        else:
            if self.dim == 1:
                data = data.flatten()
            else:
                pass

        return data

    def save_old_trajectory(self,dir_name,index):
        os.system("cp -r ./%s/output/trajectory.dat %s/traj_files/trajectory_%s.dat"%(dir_name,dir_name,index,))
        #os.system("cp -r ./%s/ext %s/ext_files/ext_%s"%(dir_name,dir_name,index,))

    def update_grid_string(self):
        grid_string = ''
        oxDNA_potential_grid = self.potential_grid / self.oxdivkT
        #solve this later 
        if self.dim == 1:
            for i in oxDNA_potential_grid:
                grid_string += f"{i},"
        elif self.dim == 2: 
            for i in oxDNA_potential_grid:
                for j in i:
                    grid_string += f"{j},"
                grid_string += '|'
        self.grid_string = grid_string

    def write_new_external_forces(self,dir_name):
        # use our p-dictionary! 
        p_string = ""
        for i in self.p_dict:
            p_string += f"{i}={self.p_dict[i]}\n"
        with open("new_ext","w+") as f:
            if self.ratio == 1:
                if self.dim == 1:
                    for mode in range(1,5):
                        our_string = f'''{{\ntype=meta_atan_com_trap\n{p_string}\n
    xmin={self.xmin}\nxmax={self.xmax}\nN_grid={self.N_grid}\npotential_grid={self.grid_string}\nmode={mode}\nPBC=0\n}}\n'''
                        f.write(our_string)

            elif self.angle == 1:
                if self.dim == 1:
                    for mode in range(1,4):
                        our_string = f'''{{\ntype=meta_com_angle_trap\n{p_string}\n
    xmin={self.xmin}\nxmax={self.xmax}\nN_grid={self.N_grid}\npotential_grid={self.grid_string}\nmode={mode}\nPBC=0\n}}\n'''
                        f.write(our_string)

            else:
                if self.dim == 1:
                    for mode in range(1,3):
                        our_string = f'''{{\ntype=meta_com_trap\n{p_string}
    xmin={self.xmin}\nxmax={self.xmax}\nN_grid={self.N_grid}\npotential_grid={self.grid_string}\nPBC=0\n
    mode={mode}\n}}\n'''
                        f.write(our_string)

                elif self.dim == 2:
                    for mode in range(1,5):
                        our_string = f'''{{\ntype=meta_2D_com_trap\n{p_string}\n
    xmin={self.xmin}\nxmax={self.xmax}\nymin={self.xmin}\nymax={self.xmax}\nN_grid={self.N_grid}\npotential_grid={self.grid_string}\nmode={mode}\nPBC=0\n}}\n'''
                        f.write(our_string)

        os.system(f'cp new_ext ./{dir_name}/ext')

    def save_potential_grid(self,index):
        with open(f"center_data/centers_{index}",'wb+') as f:
            dump(self.potential_grid,f)

    def save_positions(self,new_data,index,walker_index):
        with open(f"run-meta_{walker_index}/positions/all-pos",'w+') as f:
            dump(new_data,f)

    def interpolatePotential1D(self,x,potential_grid):
        x_left = dX*np.floor(x/self.dX)
        x_right = x_left + self.dX
        ix_left = np.floor((x-self.xmin)/self.dX)
        ix_right = ix_left + 1
        f1 = potential_grid[int(ix_left)]
        f2 = potential_grid[int(ix_right)]
        fx = (x_right - x)/self.dX*f1 + (x - x_left)/self.dX*f2
        return fx

    def interpolatePotential2D(self,x,y,potential_grid):
        dX = self.dX
        xmin = self.xmin
        x_left  = dX*np.floor(x/dX)
        y_left  = dX*np.floor(y/dX)
        x_right = x_left + dX
        y_right = y_left + dX

        ix_left  = int((x-xmin)/dX)
        iy_left  = int((y-xmin)/dX)
        ix_right = ix_left + 1
        iy_right = iy_left + 1

        f11 = potential_grid[ix_left,iy_left]
        f12 = potential_grid[ix_left,iy_right]
        f21 = potential_grid[ix_right,iy_left]
        f22 = potential_grid[ix_right,iy_right]

        fxy1 =  (x_right - x) /dX * f11 + (x - x_left) /dX * f21 
        fxy2 =  (x_right - x) /dX * f12 + (x - x_left) /dX * f22
        local_potential = (y_right-y)/dX * fxy1 + (y-y_left)/dX * fxy2
        return local_potential

    def get_metad_potential(self,x,cx):
        U = self.A*np.exp(- ((x-cx)**2)/(2*self.sigma**2) )
        return U

    def get_metad_potential2D(self,x,cx,y,cy):
        U = self.A*np.exp(- ((x-cx)**2+(y-cy)**2)/(2*self.sigma**2) )
        return U

    def initialise(self):
        if not self.continue_run:
            os.system("rm -r run-meta_*")

            for dir_name in self.dir_names: 
                os.system(f"cp -r run-meta {dir_name}")
                os.system(f"mkdir {dir_name}/traj_files")
                os.system(f"mkdir {dir_name}/ext_files")
                os.system(f"mkdir {dir_name}/positions")

            os.system("rm -r center_data")
            os.system("mkdir center_data")
            
        self.queue = mp.JoinableQueue()
        self.processes = []
        for dir_name in self.dir_names:
            self.processes.append(oxDNARunner(dir_name, "input-meta", self.conf_interval, self.queue))
            self.processes[-1].start()

    def do_metadynamics_iteration(self,index):
        print("iteration %s"%(index,))
        
        # write forces
        self.update_grid_string()
        for dir_name in self.dir_names: 
            self.write_new_external_forces(dir_name)
            
        # run parallel computation
        print ("starting processes")
        for walker_index,dir_name in enumerate(self.dir_names): 
            self.queue.put(tau)
        self.queue.join()

        if index%self.save_hills == 0:
            self.save_potential_grid(index)
        '''
        for dir_name in self.dir_names:                     
            self.save_old_trajectory(dir_name,index)
        '''
       
        #load new data 
        new_collective_data = []
        for walker_index,dir_name in enumerate(self.dir_names):
            new_data = self.get_new_data(dir_name)
            #self.save_positions(new_data, index, walker_index)
            new_collective_data.append(copy(new_data))
           
        #update forces
        if self.meta:
            for data_selection in new_collective_data:
                x = data_selection[-1]
                if self.dim == 1:
                    local_potential = self.interpolatePotential1D(x,self.potential_grid)
                    print ('x',x,'U',local_potential);

                elif self.dim == 2:
                    local_potential = self.interpolatePotential2D(x[0],x[1],self.potential_grid)

                prefactor = np.exp(-local_potential / self.dT)

                if self.dim == 1:
                    ix_left = int((x-self.xmin)/self.dX);
                    ix_start = max(ix_left - self.r_cut_int,0)
                    ix_end = min(ix_left + self.r_cut_int+1,self.N_grid)
                    for ix in range(int(ix_start),int(ix_end)):
                        x_val = self.xmin + self.dX*ix
                        new_potential = self.get_metad_potential(x,x_val)
                        self.potential_grid[ix] += prefactor*new_potential

                elif self.dim == 2:
                    ix_left = int((x[0]-self.xmin)/self.dX);
                    ix_start = max(ix_left - self.r_cut_int,0)
                    ix_end = min(ix_left + self.r_cut_int+1,self.N_grid)

                    iy_left = int((x[1]-self.xmin)/self.dX);
                    iy_start = max(iy_left - self.r_cut_int,0)
                    iy_end = min(iy_left + self.r_cut_int+1,self.N_grid)

                    for ix in range(int(ix_start),int(ix_end)):
                        for iy in range(int(iy_start),int(iy_end)):
                            x_val = self.xmin + self.dX*ix
                            y_val = self.xmin + self.dX*iy
                            new_potential = self.get_metad_potential2D(x[0],x_val,x[1],y_val)
                            self.potential_grid[ix,iy] += prefactor*new_potential

    def do_run(self):
        if self.continue_run:
            for index in range(self.max_index+1,self.max_index+1+self.Niter):
                self.do_metadynamics_iteration(index)
        else:
            for index in range(self.Niter):
                self.do_metadynamics_iteration(index)
if __name__ == '__main__': 

	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument("--A",default = 0.1)
	parser.add_argument("--sigma",default = 0.05)
	parser.add_argument("--tau",default = 10000)
	parser.add_argument("--dX", default = 0.001)
	parser.add_argument("--N_walkers",default = 6)
	parser.add_argument("--dT",default = 3)
	parser.add_argument("--dim",default = 1)
	parser.add_argument("--p_fname",default = "locs.meta")
	parser.add_argument("--ratio",default = 0)
	parser.add_argument("--angle",default = 0)
	parser.add_argument("--Niter",default = 10000)
	parser.add_argument("--xmin",default = 0)
	parser.add_argument("--xmax",default = 30)
	parser.add_argument("--conf_interval",default = int(1e3))
	parser.add_argument("--save_hills",default = 1)
	parser.add_argument("--continue_run",default = 0)
	parser.add_argument("--T",default = 295)

	args = parser.parse_args()

	A = float(args.A)
	sigma = float(args.sigma)
	tau = int(args.tau )
	dX = float(args.dX)
	N_walkers = int(args.N_walkers)
	dT = float(args.dT)
	dim = int(args.dim)
	p_fname = str(args.p_fname)
	ratio = bool(int(args.ratio))
	Niter = int(args.Niter)
	angle = bool(int(args.angle)) # add defaults to these
	xmin = float(args.xmin) # add defaults to these
	xmax = float(args.xmax) # add defaults to these
	conf_interval = int(args.conf_interval)
	save_hills = int(args.save_hills) #how frequently we save
	continue_run = bool(int(args.continue_run)) #whether we continue
	temperature = float(args.T)

	#load in the pfile here.
	with open(p_fname) as f: p_string = f.readlines()
	p_dict = {}
	for i in p_string:
	    pair = i.replace('\n','').split(':')
	    p_dict[pair[0]] = pair[1]

	estimator = Estimator(    Niter = Niter, meta = True, wt = True,dT=dT,
				  sigma=sigma,dX=dX,A=A,tau = tau,
				  N_walkers=N_walkers,
				  p_dict = p_dict, dim=dim, ratio=ratio,angle=angle,xmin=xmin,xmax=xmax,
				  conf_interval = conf_interval,save_hills=save_hills,continue_run=continue_run,
				  T = temperature)

	estimator.initialise()
	estimator.do_run()
						      
