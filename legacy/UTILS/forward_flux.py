


import os
import random
import glob
import shutil
import readers
import base
import order_parameters
import subprocess
import time

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")


class ForwardFlux:
	def __init__(self,init,stop,Astate,seed,order_pars,init_template,launchcommand,prefix='TESTffs_'):
		self.init = init
		self.stop = stop
		self.Astate = Astate
		self.order_pars = order_pars
		self.launchcommand = launchcommand
		fin = open(init_template,'r')
		self.init_template = fin.readlines()
		steps = 0;
		dt = 0;
		for line in self.init_template:
  			if "topology" in line[0:len("topology")+2]:
    				self.topologyfile = line.split('=')[1].replace(' ','').replace('\n','')
			if "dt" in line[0:len("dt")+2]:
				dt = line.split('=')[1].replace(' ','').replace('\n','')
			if "steps" in line[0:len("steps")+2]:
				steps = line.split('=')[1].replace(' ','').replace('\n','')
					
		fin.close()
		self.iterations = 0
		self.completed = 0
		self.seed = seed
		self.prefix = prefix 

 		self.timeunit = float(dt) * float(steps)
		self.runningtime = 0
		print " Time between different order paramater checks is ",self.timeunit
 	

	def run_simulation(self,inputfile,last_conf,dbg=0):
		self.success = -1
		while self.success == -1 :
			#benchmark = time.time()
			#os.system(self.launchcommand + ' ' + inputfile + ' 2>/dev/null')
			launchargs = [self.launchcommand,inputfile]	
			myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
			myinput.wait()		
			self.runningtime += self.timeunit	
			#benchmark = time.time() - benchmark
			#print ' Running program took ', benchmark
 			#benchmark = time.time()
			reader = readers.LorenzoReader(last_conf,self.topologyfile)
			system = reader.get_system()
				
			system.map_nucleotides_to_strands()
			myinput = os.popen(PROCESSDIR + 'output_bonds ' + str(inputfile ) + " " + last_conf  + ' 2>/dev/null' )
			#print myinput.readlines()
			system.read_H_bonds(myinput.readlines())
			myinput.close()
			
			#self.order_pars.set_system(system)
			
			#self.order_pars.calculate_order_parameters()
			
			if(dbg != 0):
				print 'State ',dbg, ' is ',self.order_pars.get_state(dbg,system)
			#print self.stop, ' is ', self.order_pars.get_order_parameter(self.stop)
			#print self.Astate, ' is' , self.order_pars.get_order_parameter(self.Astate)
			#print self.init, ' is ', self.order_pars.get_order_parameter(self.init)
			
			#if(self.order_pars.get_order_parameter('MINDIST_0_3') > 10.0):
			#	print "These are the following parameters: "
			#	g = self.order_pars.get_all_parameters()
			#	print g
			#print self.stop, ' is ', self.order_pars.get_order_parameter(self.stop)

			if(self.order_pars.get_state(self.Astate,system)): #This means the system fell down to Astate before reaching stop order parameter	
					self.success = 0
				
			if(self.order_pars.get_state(self.stop,system)): #This means success 
					self.success = 1
			#print 'Success is ', self.success
			#benchmark = time.time() - benchmark
			#print 'Evalloop took ',benchmark
		return self.success

	def run_multiple_runs(self,logfileprefix,restarts=-1,dbg=0):
			self.iterations = 0
			self.completed = 0
			logfile = logfileprefix + '_' + self.init + '_to_' + self.stop + '.' + str(self.seed) + '.log'
			logout = open(logfile,'a')
			self.last_conf_name = self.prefix+'tempfile'+'_'+self.stop+'_'+str(self.seed)+'.dat'
			init_filename = self.prefix + 'it_' + str(self.iterations) + 'seed_'  + str(self.seed) 

			while(restarts == -1 or self.iterations < restarts):
				self.runningtime = 0
				
				possible_inits = glob.glob(self.prefix+'_lastconf_at_'+self.init+'_*.dat')	
				if(len(possible_inits) < 1):
					raise IOError('Error, cannot find any configuration file to start from '+ self.prefix+'_lastconf_at_'+self.init+'_*.dat' )

				 
				initfile = random.choice(possible_inits)
				print 'Starting simulation from ',initfile
				shutil.copy(initfile,self.last_conf_name)
				
				self.gen_init_file(self.seed,init_filename,self.last_conf_name,self.last_conf_name+'.trajectory')
				
				result = self.run_simulation(init_filename,self.last_conf_name,dbg  )
				
				if result == 1: #succes!
					self.completed += 1
					confname = self.prefix + '_lastconf_at_'+self.stop+'_'+str(self.seed)+'.'+str(self.iterations)+'.dat'
					k = 1
					while os.path.exists(confname):
						print 'Error, filename ' +confname+ ' exists, trying different one',
						confname = self.prefix + '_lastconf_at_'+self.stop+'_'+str(self.seed)+'.'+str(self.iterations)+'_sub'+str(k)+'.dat'
						k += 1
							
					shutil.copy(self.last_conf_name,confname)
					self.completed += 1
					print 'Sucess, reached order parameter ',self.stop ,' and configuration saved to ',confname	
				else:
					print 'Fell back to order parameter ',self.Astate
										
				self.seed += 1
				self.iterations +=1
				logout.write(self.init + ' ' + self.stop + ' '+ ' ' + str(result) + ' ' + str(self.runningtime) +'\n')	
				logout.flush()			
 			print 'Finished simulations, completed / runs = ', float(self.completed)/self.iterations

	def run_flux(self, order_parameter_to_reach, order_parameter_to_check_recross,start_configuration_filename ,number_of_crossings,dbg=0):
			input_filename = self.prefix + '_input_' + 'it_' + str('flux') + 'seed_'  + str(self.seed) 
			self.runningtime = 0
			self.last_conf_name = self.prefix+'fluxtempfile'+'_'+order_parameter_to_reach+'_'+str(self.seed)+'.dat' 
			shutil.copy(start_configuration_filename,self.last_conf_name)
			self.gen_init_file(self.seed,input_filename,self.last_conf_name,'/dev/null')
			
			logfile = self.prefix + '_' + self.init + '_fluxto_' + order_parameter_to_reach + '.log'
			logout = open(logfile,'a')
			crossings = 0
			time = 0
			
			while crossings < number_of_crossings:	
				self.stop = order_parameter_to_reach
				#self.last_conf_name = self.prefix+'fluxtempfile'+'_'+str(crossings)+'_'+self.stop+'_'+sefl.seed+'.dat'
				res = 0				
				while(res == 0):
					res = self.run_simulation(input_filename,self.last_conf_name,dbg)
				print 'Crossed ',order_parameter_to_reach
				logout.write(str(crossings)+':  crossed through '+ str(order_parameter_to_reach) + ' at '+str(self.runningtime) + '\n')
				logout.flush()
				crossings += 1
				confname = self.prefix + '_lastconf_at_'+self.stop+'_'+str(self.seed)+'.t'+str(self.runningtime)+'_N'+str(crossings)+'.dat'
				print 'Generated ',confname
				shutil.copy(self.last_conf_name,confname)
				if crossings >= number_of_crossings:
					return
				self.stop = order_parameter_to_check_recross
				res = 0
				while( res == 0):
					res = self.run_simulation(input_filename,self.last_conf_name,dbg)
				
				
				if(res == 1):
					print ' Recrossed to ',order_parameter_to_check_recross
			logout.close()

	def gen_init_file(self,seed,init_filename,last_conf_name,trajectory_name):
			#if(os.path.exists(init_filename)):
				#raise IOError('Error, trying to overwrite ' + init_filename)

			newfile = open(init_filename,'w')
			for line in self.init_template:
				newline = line.replace('$(conf_file)',last_conf_name).replace('$(seed)',str(seed)).replace('$(last_conf)',last_conf_name).replace('$(trajectory_name)',trajectory_name)
				newfile.write(newline)
								


