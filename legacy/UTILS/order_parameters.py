import base
try:
	import numpy as np
except:
	import mynumpy as np
import sys
import readers
import os
import subprocess
import math


def l2norm(vec):
	return math.sqrt(np.dot(vec,vec))


#PROCESSDIR = os.path.join(os.path.dirname(__file__),"process_data/")
PROCESSPROGRAM = os.path.join(os.path.dirname(__file__), "output_bonds.py")


class NEAR_HBONDS_Parameter:
	def __init__(self,param_command):
		self.command = param_command
		vals = param_command.split()
		if len(param_command) < 5:
			raise IOError('Order parameter ' + param_command + ' is wrong. Use id sid1 nucleotide1 sid2 nucleotide2')
		self.name = vals[0];
		self.sid1 = int(vals[1])
		self.nuc1 = int(vals[2])
		self.sid2 = int(vals[3])
		self.nuc2 = int(vals[4])
		self.state = 0 

	def get_state (self, system):
		self.nucleotide_id1 = system._strands[self.sid1]._nucleotides[self.nuc1].index
		self.nucleotide_id2 = system._strands[self.sid2]._nucleotides[self.nuc2].index
		if system.get_interaction (self.nucleotide_id1, self.nucleotide_id2, base.INT_HYDR_NEAR) == 'A':
			self.state = 1
		else:
			self.state = 0
		return self.state

	def get_type(self):
		return "NEAR_HBONDS"


class HBONDS_Parameter:
	def __init__ (self, param_command):
		self.command = param_command
		vals = param_command.split()
		if len(param_command) < 5:
			raise IOError('Order parameter ' + param_command + ' is wrong. Use id sid1 nucleotide1 sid2 nucleotide2')
		self.name = vals[0];
		self.sid1 = int(vals[1]);
		self.nuc1 = int(vals[2]);
		self.sid2 = int(vals[3]);
		self.nuc2 = int(vals[4]);
		self.state = 0

	def get_state(self, system):
		self.nucleotide_id1 = system._strands[self.sid1]._nucleotides[self.nuc1].index
		self.nucleotide_id2 = system._strands[self.sid2]._nucleotides[self.nuc2].index
		if system.check_H_interaction(self.nucleotide_id1,self.nucleotide_id2):
			self.state = 1
		else:
			self.state = 0
		return self.state

	def get_type(self):
		return "HBONDS"

class DISTANCE_Parameter:
	def __init__(self,param_command):
		self.command = param_command
		vals = param_command.split()
		if len(param_command) < 5:
			raise IOError('Order parameter ' + param_command + ' is wrong. Use id sid1 nucleotide1 sid2 nucleotide2 dist')
		self.name = vals[0];
		self.sid1 = int(vals[1]);
		self.nuc1 = int(vals[2]);
		self.sid2 = int(vals[3]);
		self.nuc2 = int(vals[4]);
		
		self.state = -1 

	def get_state(self,system):
		self.state = l2norm(system._strands[self.sid1]._nucleotides[self.nuc1].distance(system._strands[self.sid2]._nucleotides[self.nuc2],True,system._box))
		return self.state

	def get_type(self):
		return "DISTANCE"	 


class BASEDISTANCE_Parameter:
	def __init__(self,param_command):
		self.command = param_command
		vals = param_command.split()
		if len(param_command) < 5:
			raise IOError('Order parameter ' + param_command + ' is wrong. Use id sid1 nucleotide1 sid2 nucleotide2 dist')
		self.name = vals[0];
		self.sid1 = int(vals[1]);
		self.nuc1 = int(vals[2]);
		self.sid2 = int(vals[3]);
		self.nuc2 = int(vals[4]);
		
		self.state = -1 

	def get_state(self,system):
		
        	dr = system._strands[self.sid1]._nucleotides[self.nuc1].get_pos_base()  -  system._strands[self.sid2]._nucleotides[self.nuc2].get_pos_base()
       		dr -= system._box * np.rint (dr / system._box)
		self.state = l2norm(dr)
		return self.state

	def get_type(self):
		return "BASEDISTANCE"	 


class ALL_HBONDS_Parameter:
	def __init__(self,param_command):
		self.state = 0
		vals = param_command.split()
		self.name = vals[0];
		self.sid1 = int(vals[1])
		self.sid2 = int(vals[2])
	
	def get_state(self,system):
		inters =  system._strands[self.sid1].get_H_interactions()
		if self.sid2 in inters.keys():
			self.state = inters[self.sid2]
		else: 
			self.state = 0
		return self.state

	def get_type(self):
		return "ALL_BONDS"


class EVALUATED_Parameter:
	def __init__(self,param_command,par_dictionary):
		for par in par_dictionary.keys():
			param_command = param_command.replace('$('+par+')',par+'.get_state(_system)')
			#print par, par_dictionary[par], param_command
		self.command = param_command
		self.state  = 0
		self.par_dictionary = {}
		for key, val in par_dictionary.items():
			self.par_dictionary[key] = val 

	def get_type(self):
		return "EVALUATED"

	def get_state(self,system):
		self.par_dictionary['_system'] = system
		self.state = eval(self.command,{},self.par_dictionary)
		return self.state


class OrderParameters:
	def __init__(self, *args):
		self.parameters = {}
		for pars in args:
			if isinstance (pars, dict) or isinstance(pars, OrderParameters):
				for key, val in pars.items():
					if key in self.parameters.keys():
						raise IOError('Error, name conflict, '+key+', order paramter declared twice')
					else:
						self.parameters[key] = val
			elif isinstance(pars, str):
				self.read_order_parameters_from_file(pars)
			else:	
				raise IOError('Invalid format of input argument in OrderParameters init')

	def get_order_parameter(self, parameter_id):
		if(parameter_id in self.evaluated_parameters.keys()):
			return self.evaluated_parameters[parameter_id].get_state()
		else:
			return self.parameters[parameter_id].get_state()

	def get_all_states (self,system):
		self.all_evaluated  = {}
		for key,val in self.parameters.items ():
			self.all_evaluated[key] = val.get_state(system)
		return self.all_evaluated

	def get_state(self, key, system):
		return self.parameters[key].get_state(system)
	
	def __getitem__(self,key):
		return self.parameters[key]

	def items(self):
		return self.parameters.items()

	def keys(self): 
		return self.parameters.keys()

	#Those are parameter processing functions for loading from file	
	def add_order_parameter(self,ordertype,line):
		pamid = line.split()[0]
		if(pamid in self.parameters.keys()):
			raise IOError('Error ' + pamid + 'declared twice')

		if ordertype == "HBONDS":
			nparam = HBONDS_Parameter(line)
			self.parameters[pamid] = nparam
		elif ordertype == "DISTANCE":
			nparam = DISTANCE_Parameter(line)
			self.parameters[pamid] = nparam
		elif ordertype == "NEAR_HBONDS":
			nparam = NEAR_HBONDS_Parameter(line)
			self.parameters[pamid] = nparam
		elif ordertype == 'ALL_HBONDS':
			nparam = ALL_HBONDS_Parameter(line)
			self.parameters[pamid] = nparam
		elif ordertype == "EVALUATED":
			parameter = line.rstrip().replace(pamid,' ')
						
			for par in self.parameters.keys():
					parameter = parameter.replace('$('+par+')',par+'.get_state(_system)')
			nparam = EVALUATED_Parameter(parameter,pars)
			self.parameters[pamid] = nparam
		else:
			raise IOError('Wrong format of the order parameter file')
		#print 'Added ',line

	def read_order_parameters_from_file(self, filename):
		input = open(filename,'r')
		ordertype = 0
		for line in input.readlines():
			if( (not line.isspace()) and (line.strip()[0] != '#'    ) ):
					if( line.strip() == "HBONDS"):
						ordertype = "HBONDS"
					elif (line.strip() == "DISTANCE"):
						ordertype = "DISTANCE"
					elif (line.strip() == "EVALUATED" ):
						ordertype = "EVALUATED"
					elif (line.strip() =="NEAR_HBONDS"):
						ordertype = "NEAR_HBONDS"
					elif (linet.strip() == 'ALL_HBONDS'):
						ordertype = 'ALL_HBONDS'
					else:	
						self.add_order_parameter(ordertype,line)


class EVALUATED_Weight:
	def __init__(self,name,param_command,params):
		if not isinstance(params,OrderParameters) and not isinstance(params,dict):
			raise IOError('Incorrect weight')
		self.name = name
		self.params = params
		self.key_command = param_command	
		for par in params.keys():
			param_command = param_command.replace('$('+par+')',par+'.get_state(_system)')
			self.key_command = self.key_command.replace('$('+par+')',' ' + par +' ' )
		self.command = param_command
		self.par_dictionary = {}
		for key,val in params.items():
			self.par_dictionary[key] = val 

	def get_value(self,system,order_par_obj):
		self.par_dictionary['_system'] = system
		self.state = eval(self.command,{},self.par_dictionary)
		return self.state

	def get_value_from_key(self,keys):
		self.state = eval(self.key_command,{},keys)
		return self.state


class ARRAY_Weight:
	def __init__(self,name,order_par,vals):
		self.name = name
		self.order_par = order_par
		self.vals = vals
	def get_value(self,system,order_par_obj):
		return self.vals[order_par_obj.get_state(self.order_par,system)]
	def get_value_from_key(self,keys):
		return self.vals[keys[self.order_par]]

class DICT_Weight:
	def __init__ (self,name, dict_vals, index_names):
		self.name = name
		if isinstance (index_names, list) or isinstance (index_names, tuple):
			self.index_names = index_names
		elif isinstance (index_names, str):
			self.index_names = [index_names]
		else:
			raise ValueError ("Could not build DICT_Weight")
		self.dict_vals = dict_vals
	
	def get_value (self, system, order_par_obj):
		indextuple = []
		for par in self.index_names:
			indextuple.append((par, order_par_obj.get_state (par, system)))
		key = tuple (indextuple)
		#print "KEY:", key
		try:
			result = self.dict_vals [key]
		except KeyError:
			result = max (self.dict_vals.itervalues ()) # we return the maximum value
			#print order_par_obj.get_state (par, system)
			#print >> sys.stderr, "got unknown key", key
		if not (isinstance(result,int) or isinstance(result,float)):
			print 'Warning, weight ', self.name, '  inadmissible index detected ', indextuple, 'using weight 1'
			result = 1.0
		return result

	def get_value_from_key(self, keys):
		indextuple = []
		for par in self.index_names:
			#print par
			indextuple.append((par, keys[par]))
#		print '####', indextuple
		try:
			result = self.dict_vals[tuple(indextuple)]
		except KeyError:
			result = max (self.dict_vals.values())
		if not (isinstance(result,int) or isinstance(result,float)):
			print 'Warning, weight ', self.name, '  inadmissible index detected ', indextuple, 'using weight 1'
			result = 1.0
		return result


class Weights:
	def __init__(self, * weight_dicts):
		self.weights = {}
		for weidict in weight_dicts:
			if isinstance (weidict, dict):
				for key, val in weidict.items():
					self.weights[key] = val
			elif isinstance (weidict, list) or isinstance (weidict, tuple):
				for w in weidict:
					if not isinstance (w, dict):
						raise ValueError
					for key, val in w.items():
						self.weights[key] = val
			else:
				raise ValueError

	def get_value (self, weight, system, order_parameters):
		if weight is not None:
			return self.weights[weight].get_value(system,order_parameters)
		else:
			ret = 1.
			for key in self.weights.keys ():
				ret = ret * self.weights[key].get_value (system, order_parameters)
			return ret

	def get_value_from_key (self, weight, keys):
		if weight is not None:
			return self.weights[weight].get_value_from_key(keys)
		else:
			ret = 1.
			for key in self.weights.keys ():
#				print "###", key
				ret = ret * self.weights[key].get_value_from_key(keys)
			return ret

	def get_all_factors(self,system,order_parameters):
		values = 1.0
		for key,val in self.weights.items():
			values *= val.get_value(system,order_parameters)
		return values

	def get_all_values(self,system,order_param_obj):
		self.eval_weights = {}
		for key,val in self.weights.items():
			self.eval_weights[key] = val.get_value(system,order_param_obj)
		return self.eval_weights

	def keys(self):
		return self.weights.keys()
	def items(self):
		return self.weights.items()

def define_weight(weight_name, opindex_names, vals):
	if isinstance(vals,list):
		weight = ARRAY_Weight(weight_name, opindex_names, vals)
	elif isinstance(vals, dict):
		weight = DICT_Weight (weight_name, vals, opindex_names)
	elif isinstance(opindex_names,str) and (isinstance(vals,OrderParameters) or isinstance(vals,dict)):
		weight = EVALUATED_Weight(weight_name,opindex_names,vals)
	else:
		raise IOError ('Invalid format to specify weight')
	return {weight_name:weight}

def load_system (inputfile, conffile, topologyfile):
	reader = readers.LorenzoReader (conffile, topologyfile)
	system = reader.get_system ()
	if system == False:
		raise IOError('Crashed when trying to read configuration file '+str(conffile)+' and topology file '+str(topologyfile))
	system.map_nucleotides_to_strands()
	#print 'LAUNCHING',PROCESSDIR + "output_bonds.py"
	#os.system(PROCESSDIR + "output_bonds.py")
	if not os.path.isfile (PROCESSPROGRAM):
		raise IOError ("Cannot execute output_bonds program "+PROCESSPROGRAM)
	launchargs = [PROCESSPROGRAM,inputfile,conffile]
	myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if myinput.wait () != 0:
		for line in myinput.stderr.readlines():
			print line,
		raise IOError ('Critical Error')
	system.read_all_interactions(myinput.stdout.readlines())
	return system

def load_energies_to_system (inputfile, conffile,system,iteration_id=0):
		system.map_nucleotides_to_strands()
		if not os.path.isfile(PROCESSPROGRAM):
			raise IOError ("Cannot execute output_bonds program")
		launchargs = [PROCESSPROGRAM,inputfile,conffile,iteration_id]
		myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		myinput.wait ()
		system.read_all_interactions(myinput.stdout.readlines())
		for line in myinput.stderr.readlines():
			if "CRITICAL" in line:
				print "Critical error"
				raise IOError
		return system	

def complementary_bonds_between_strands(system,name,strand_id1,strand_id2,start=0,offset=0,include_all_parameters=False):
	#by default, this makes a bond between the first nucleotide of strand_id1 with the last bond of strand_id2
	subpars = {}
	#out.write('HBONDS \n')
	names = []	
	for i in range(len(system._strands[strand_id1]._nucleotides) - start - offset):
		index_a = i + start
		index_b = len(system._strands[strand_id2]._nucleotides) - offset - i - 1
		if(index_a >= len(system._strands[strand_id1]._nucleotides) or index_b < 0):
			print 'Reached impossible indexes', index_a, ' ' , index_b, ' not including them'
			continue
		bondname = name + 'cHB_S' + str(strand_id1) + 'N' + str(index_a) + 'S' + str(strand_id2) + 'N' + str(index_b)
		param = HBONDS_Parameter('%s %d %d %d %d' % (bondname,strand_id1,index_a,strand_id2,index_b) )
		subpars[bondname] = param
		names.append(bondname)
		if(system._strands[strand_id1]._nucleotides[index_a]._base + system._strands[strand_id2]._nucleotides[index_b]._base != 3):
			print 'Warning, bases ',index_a ,' and ', index_b, ' on the strands are not complementary'

	evalpam =  ' $(' + names[0]+')'	
	for bond in names[1:]:
		evalpam += ' + $('+ bond +')'
	if include_all_parameters:
		subpars[name] = EVALUATED_Parameter(evalpam,subpars)
		return subpars
	else:
		return {name:EVALUATED_Parameter(evalpam,subpars)}

def all_bonds_between_strands(system,name,strand_id1,strand_id2):
	if strand_id1 >= system._N or strand_id2 >= system._N:
		print 'Error while constructing order parameter, strand index is incorrect'
		raise IOError('Incorrect index')
	else:
		return {name:ALL_HBONDS_Parameter(name+' '+str(strand_id1) + ' ' +str(strand_id2))}


def array_of_complementary_bonds_between_strands(system,name,strand_id1,strand_id2,start=0,offset=0):
	#by default, this makes a bond between the first nucleotide of strand_id1 with the last bond of strand_id2
        subpars = {}
	#out.write('HBONDS \n')
	names = []
	params = []	
	for i in range(len(system._strands[strand_id1]._nucleotides) - start - offset):
		index_a = i + start
		index_b = len(system._strands[strand_id2]._nucleotides) - offset - i - 1
		if(index_a >= len(system._strands[strand_id1]._nucleotides) or index_b < 0):
			print 'Reached impossible indexes', index_a, ' ' , index_b, ' not including them'
			continue
		bondname = name + 'cHB_S' + str(strand_id1) + 'N' + str(index_a) + 'S' + str(strand_id2) + 'N' + str(index_b)
		param = HBONDS_Parameter('%s %d %d %d %d' % (bondname,strand_id1,index_a,strand_id2,index_b) )
		subpars[bondname] = param
		params.append(param)
		names.append(bondname)
		if(system._strands[strand_id1]._nucleotides[index_a]._base + system._strands[strand_id2]._nucleotides[index_b]._base != 3):
			print 'Warning, bases ',index_a ,' and ', index_b, ' on the strands are not complementary'


	evalpam =  ' $(' + names[0]+')'	
	for bond in names[1:]:
		evalpam += ' + $('+ bond +')' 
	
	subpars[name] = EVALUATED_Parameter(evalpam,subpars)
	return subpars,params		

	

def minimum_distance_between_num_bases(system,name,strand_id1,strand_id2,start_A=0,offset_B=0,number_of_bases=1,include_all_parameters=False):
	#by default, this makes a distance between the first nucleotide of strand_id1 with the last bond of strand_id2
	names = []
	subpars = {}	
	for i in range(number_of_bases):
		index_a = i + start_A
		index_b = len(system._strands[strand_id2]._nucleotides) - offset_B - i - 1
		if(index_a >= len(system._strands[strand_id1]._nucleotides) or index_b < 0):
			print 'Reached impossible indexes', index_a, ' ' , index_b, ' not including them'
			continue
		bondname = name + '_D_S' + str(strand_id1) + 'N' + str(index_a) + 'S' + str(strand_id2) + 'N' + str(index_b)
		subpars[bondname] = BASEDISTANCE_Parameter('%s %d %d %d %d ' % (bondname,strand_id1,index_a,strand_id2,index_b) ) 
		names.append(bondname)
		#print strand_id1, ' ' , strand_id2, ' ',index_a, index_b
		if(system._strands[strand_id1]._nucleotides[index_a]._base + system._strands[strand_id2]._nucleotides[index_b]._base != 3):
			print 'Warning, bases ',index_a ,' and ', index_b, ' on the strands are not complementary'

	evalpam =  ' min( $(' + names[0]+')'	
	for bond in names[1:]:
		evalpam += ' , $('+ bond +')'
	evalpam += ')'
 
	if include_all_parameters:
		subpars[name] = EVALUATED_Parameter(evalpam,subpars)
		return subpars
	else:
		return {name:EVALUATED_Parameter(evalpam,subpars)}


def minimum_distance_between_complementary_bases(system,name,strand_id1,strand_id2,start=0,offset=0,include_all_parameters=False):
	#by default, this makes a distance between the first nucleotide of strand_id1 with the last bond of strand_id2
	names = []
	subpars = {}	
	for i in range(len(system._strands[strand_id1]._nucleotides) - start - offset):
		index_a = i + start
		index_b = len(system._strands[strand_id2]._nucleotides) - offset - i - 1
		if(index_a >= len(system._strands[strand_id1]._nucleotides) or index_b < 0):
			print 'Reached impossible indexes', index_a, ' ' , index_b, ' not including them'
			continue
		bondname = name + '_D_S' + str(strand_id1) + 'N' + str(index_a) + 'S' + str(strand_id2) + 'N' + str(index_b)
		subpars[bondname] = DISTANCE_Parameter('%s %d %d %d %d ' % (bondname,strand_id1,index_a,strand_id2,index_b) ) 
		names.append(bondname)
		#print strand_id1, ' ' , strand_id2, ' ',index_a, index_b
		if(system._strands[strand_id1]._nucleotides[index_a]._base + system._strands[strand_id2]._nucleotides[index_b]._base != 3):
			print 'Warning, bases ',index_a ,' and ', index_b, ' on the strands are not complementary'

	evalpam =  ' min( $(' + names[0]+')'	
	for bond in names[1:]:
		evalpam += ' , $('+ bond +')'
	evalpam += ')'
 
	if include_all_parameters:
		subpars[name] = EVALUATED_Parameter(evalpam,subpars)
		return subpars
	else:
		return {name:EVALUATED_Parameter(evalpam,subpars)}


def minimum_distance_between_bases(system,name,indices,include_all_parameters=False):
	#by default, this makes a distance between the first nucleotide of strand_id1 with the last bond of strand_id2
	names = []
	subpars = {}	
	for i in indices:
		strand_id1 = i[0]
		index_a = i[1]
		strand_id2 = i[2]
		index_b = i[3]
		if(index_a >= len(system._strands[strand_id1]._nucleotides) or index_b < 0):
			print 'Reached impossible indexes', index_a, ' ' , index_b, ' not including them'
			continue
		bondname = name + '_D_S' + str(strand_id1) + 'N' + str(index_a) + 'S' + str(strand_id2) + 'N' + str(index_b)
		subpars[bondname] = BASEDISTANCE_Parameter('%s %d %d %d %d ' % (bondname,strand_id1,index_a,strand_id2,index_b) ) 
		names.append(bondname)
		#print strand_id1, ' ' , strand_id2, ' ',index_a, index_b
		if(system._strands[strand_id1]._nucleotides[index_a]._base + system._strands[strand_id2]._nucleotides[index_b]._base != 3):
			print 'Warning, bases ',index_a ,' and ', index_b, ' on the strands are not complementary'

	evalpam =  ' min( $(' + names[0]+')'	
	for bond in names[1:]:
		evalpam += ' , $('+ bond +')'
	evalpam += ')'
 
	if include_all_parameters:
		subpars[name] = EVALUATED_Parameter(evalpam,subpars)
		return subpars
	else:
		return {name:EVALUATED_Parameter(evalpam,subpars)}



###############################################################################
# Histogram class                                                             #
###############################################################################
class Histo:
    def __init__ (self):
        self.data = {}
        self.unbiased = {}
        self.dsum = 0.
        self.usum = 0.

    def flatten (self):
        flat = []
        for k in self.data.keys():
            row = []
            for s1, s2 in k:
                if isinstance (s2, tuple):
                    for x in s2:
                        row.append (x)
                else:
                    row.append (s2)
            row.append (self.data[k])
            row.append (self.unbiased[k])
            flat.append (row)
        
        # a bit of sorting
        flat.sort()
        return flat
     
    def flatten_old (self):
        # we need an array with the names of the order parameters
        unames = []
        for key in self.data.keys():
            for subkey in key:
                if unames.count (subkey[0]) == 0:
                    unames.append (subkey[0])
        unames.sort () # we always sort them
        
        # let's build a dictionary with the maximum and minimum
        # values assumed by the order parameters
        maxs, mins = {}, {}
        for u in unames:
            maxs[u] = None
        #maxs = [None for u in unames] # anything is > None
        for key, val in self.data.items():
            for subkey in key:
                print '@@', subkey
                if subkey[1] > maxs[unames.index(subkey[0])]:
                    maxs[unames.index(subkey[0])] = subkey[1]
        for u in unames:
            mins[u] = maxs[u]
        #mins = [i for i in maxs]
        for key, val in self.data.items():
            for subkey in key:
                if subkey[1] < mins[unames.index(subkey[0])]:
                    mins[unames.index(subkey[0])] = subkey[1]
        nrows = 1
        for i in unames:
            nrows *= (maxs[u] - mins[u] + 1)
        
        flat = []
        
        for i in xrange (nrows):
            toappend = []
            larger = 1
            for j in xrange (len (unames)):
                interval = maxs[j] - mins[j] + 1
                toappend.append (mins[j] + (i/larger) % interval)
                larger *= interval
            
            mykey = tuple([(unames[k], toappend[k]) for k in xrange (len (unames))])
            try:
                toappend.append(self.data[mykey])
                toappend.append(self.unbiased[mykey])
            except:
                toappend.append(0.)
                toappend.append(0.)
            
            flat.append(toappend)
        
        return flat
    
    def reset (self):
        # we do not reset the unbiased data here
        # and we do not want to erase the keys (just the values)
        for key in self.data:
            self.data[key] = 0.
        self.dsum = 0
    
    def deep_reset (self):
        # we do reset the unbiased data here
        self.data = {}
        self.unbiased = {}
        self.dsum = 0
        self.usum = 0
    
    def update (self, mydict, we=1.):
        if not isinstance (mydict, dict):
            raise ValueError
        hkeys = mydict.keys()
        hkeys.sort()
        hkey = tuple([(key, mydict[key]) for key in hkeys])
        if hkey not in self.data:
            self.data[hkey] = 0.
            self.unbiased[hkey] = 0.
        self.data[hkey] += 1.
        self.unbiased[hkey] += 1./we
        self.dsum += 1.
        self.usum += 1./we
    
    def print_n (self, outfile=sys.stdout, mode='w'):
        if type (outfile) not in [str, file]:
            print >> sys.stderr, "Cowardly refusing to print histogram to random stuff. Soprassiedo..."
        if isinstance (outfile, str):
            out = open (outfile, mode)
        else:
            out = outfile
        
        flat = self.flatten ()
        
        # we weant to normalise the histogram
        for row in flat:
            for e in row[0:-2]:
                print >> out, e,
            print >> out, long(row[-2]), row[-2] / self.dsum, row[-1] / self.usum
        
        return
    
    def project (self, keylist):
        # here we project the histogram integrating along
        # all dimensions that are not listed in the keylist
        if type (keylist) not in [str, list]:
            raise ValueError
        
        keylist = list (keylist)
        
        unames = []
        for key in self.data.keys():
            for subkey in key:
                if unames.count (subkey[0]) == 0:
                    unames.append (subkey[0])
        unames.sort () # we always sort them
        maxs = [None for u in unames] # anything is > None
        for key, val in self.data.items():
            for subkey in key:
                if subkey[1] > maxs[unames.index(subkey[0])]:
                    maxs[unames.index(subkey[0])] = subkey[1]
        mins = [i for i in maxs]
        for key, val in self.data.items():
            for subkey in key:
                if subkey[1] < mins[unames.index(subkey[0])]:
                    mins[unames.index(subkey[0])] = subkey[1]
        
        prdata = {}
        for key, val in self.data.items ():
            prkey = tuple([(key[i][0], key[i][1]) for i in tokeep])
            prdata[prkey] += val
        
        return prdata

    def __getitem__ (self, key):
        return self.data[key]

    def items (self):
        return self.data.items()

    def keys (self):
        return self.data.keys()

## builds histogram from trajectory
def buildHisto (inputfile, conffile, topologyfile, ops, skip=1):
    import random
    output_bonds_exe = PROCESSPROGRAM
    # check for presence of files
    if not os.path.isfile (output_bonds_exe):
        raise IOError ("Cannot find output_bonds program")
    for myfile in [inputfile, conffile, topologyfile]:
        if not os.path.isfile (myfile):
            raise ValueError ('Histo.build(): could not open ' + myfile)
    if not isinstance (ops, OrderParameters):
        raise TypeError ('Histo.buidl(): last argument must be OrderParameters instance')
    if not isinstance (skip, int):
        raise TypeError ('Histo.buidl(): skip optional argument must be integer >= 1')
    if not skip >= 1:
        raise ValueError ('Histo.buidl(): skip optional argument must be integer >=1')
    # we scan through the trajectory, rebuild each configuration,
    # feed it to load_system and than process it with order_parameters
    rstr = ''.join(random.choice("abcdefg") for x in range(10))
    tmpcnf = '__' + rstr + '__.cnf'
    tmptop = '__' + rstr + '__.top'
    while os.path.isfile (tmpcnf) or os.path.isfile (tmptop):
        rstr = ''.join(random.choice("abcdefg") for x in range(10))
        tmpcnf = '__' + rstr + '__.cnf'
        tmptop = '__' + rstr + '__.top'

    ret = Histo()
    l = readers.LorenzoReader(conffile, topologyfile)
    s = l.get_system()
    nconf = 1
    while s:
        if nconf % skip == 0:
            s.print_lorenzo_output (tmpcnf, tmptop)
            s.map_nucleotides_to_strands ()
            launchargs = [output_bonds_exe, inputfile, tmpcnf]
            myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if myinput.wait() != 0:
                raise IOError ('Critical Error')
            s.read_all_interactions (myinput.stdout.readlines ())
            
            state = ops.get_all_states (s)
            ret.update (state, 1.)

            #print >> sys.stderr, '#', nconf, 'done'
        
        nconf += 1
        # get next system
        s = l.get_system ()
    os.remove (tmpcnf)
    os.remove (tmptop)
    return ret


####################################################
# here we adapt our weights; given one histogram, 
# we compute the new weights and implement them as 
# a dictionary
####################################################
def adaptive_weights (histo, oldweights=None, carried=[]):
    # we check that the arguments are sensible
    if oldweights:
        try:
            oldweights.get_value
        except AttributeError:
            raise TypeError ('Second argument to adaptive_weights is not a Weights instance')
    if not isinstance (carried, list):
        raise TypeError ('Third argument to adaptive_weights must be a list of Weights instances')
    if carried:
        for c in carried:
            try:
                c.get_value
                c.get_value_from_key
            except AttributeError:
                raise TypeError ('Third argument to adaptive_weights must be a list of Weights instances')
    
    unames = []
    newdict = {}
    for key, val in histo.unbiased.items():
        # val is never zero, otherwise key does not exist
        newdict[key] = (1./ val) * histo.usum
        for subkey in key:
            if unames.count (subkey[0]) == 0:
                unames.append (subkey[0])

        # we don't want to change the weights too much
        if oldweights is not None:
            mystate = {}
            for subkey in key: #keys are sorted, no need to sort them again
                mystate[subkey[0]] = subkey[1]
            owe = oldweights.get_value_from_key(None, mystate)
            if newdict[key] > 2. * owe:
                newdict[key] = 2. * owe
            if newdict[key] < 0.5 * owe:
                newdict[key] = 0.5 * owe
        
    unames.sort () # we always sort the names

    # we check that if the carried weight have something to say about a given
    # state, we respect it
    for key in newdict:
        mystate = {}
        for c in carried:
            mystate = {}
            for subkey in key: #keys are sorted, no need to sort them again
                mystate[subkey[0]] = subkey[1]
            owe = c.get_value_from_key(None, mystate)
            newdict[key] *= owe
    
    w1 = define_weight('W1', unames, newdict)

    return Weights (w1, [c.weights for c in carried])

