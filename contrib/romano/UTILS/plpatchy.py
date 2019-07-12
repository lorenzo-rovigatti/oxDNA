#!/usr/bin/env python

#This file loads patchy particle file from topology and Configuration
import sys
import numpy as np
import copy
import random
from __builtin__ import str
import math




myepsilon = 0.00001
def l2norm(v):
    return np.sqrt(np.dot(v,v))

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2)
    b, c, d = -axis*math.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


class Patch:
  
    def __init__(self,type = None,color=None,relposition=None,a1=None,a2=None):
        self._position = relposition
        self._a1 = a1
        self._a2 = a2
        self._type = type
        self._color = color
       
    def get_abs_position(self,r):
        return r + self._position
    
    def save_to_string(self):
        #print self._type,self._type,self._color,1.0,self._position,self._a1,self._a2
        
        outs = 'patch_%d = {\n  id = %d \n color = %d \n strength = %f \n position = %f,%f,%f  \n a1 = %f,%f,%f \n  a2 = %f,%f,%f \n } \n' % (self._type,self._type,self._color,1.0,self._position[0],self._position[1],self._position[2],self._a1[0],self._a1[1],self._a1[2],self._a2[0],self._a2[1],self._a2[2])
        return outs
       
    def init_from_string(self,lines):
        for line in lines:
            line = line.strip()
            if len(line) > 1 and line[0] != '#':
               if "id" in line:
                   vals = int(line.split('=')[1])
                   self._type = vals
               if "color" in line:
                   vals = int(line.split('=')[1])
                   self._color = vals    
               elif "a1" in line:
                    vals = line.split('=')[1]
                    x,y,z = [float(g) for g in vals.split(',')]
                    self._a1 = np.array([x,y,z])
               elif "a2" in line:
                    vals = line.split('=')[1]
                    x,y,z = [float(g) for g in vals.split(',')]
                    self._a2 = np.array([x,y,z])
               elif "position" in line:
                    vals = line.split('=')[1]
                    x,y,z = [float(g) for g in vals.split(',')]
                    self._position = np.array([x,y,z])
                         

class PLPatchyParticle:
   all_letters = ['C','H','O','N','P','S','F','K','I','Y']
   def  __init__(self,type=0,index_=0,position=np.array([0.,0.,0.]),radius=0.5):
    self.cm_pos = position
    self.index = index_
    self._radius = radius
    self._patches = None
    self._type = type
    self._patch_ids = None
    self.v = np.array([0.,0.,0.])
    self.L = np.array([0.,0.,0.])
    self.a1 = None
    self.a3 = None
   
   def set_radius(self,radius):
       self._radius = radius
 
   def add_patches(self,patches):
       self._patches = patches
   
   def translate(self,dir):
       self.cm_pos += dir
       
   def rotate(self,rot_matrix):
       self.a1 = np.dot(rot_matrix,self.a1)
       self.a3 = np.dot(rot_matrix,self.a3)
   
   def set_random_orientation(self):
       self.a1 = np.array(np.random.random(3))
       self.a1 = self.a1 / np.sqrt (np.dot(self.a1,self.a1))
       x = np.random.random(3)
       self.a3 = x - np.dot(self.a1,x) * self.a1
       self.a3 = self.a3 /  np.sqrt (np.dot(self.a3,self.a3))
       if (abs(np.dot(self.a1,self.a3)) > 1e-10):
           raise IOError("Could not generate random orientation?")
       
       
   def get_patch_position(self,patchid):
       p = self._patches[patchid]
       return p._position[0] * self.a1 + p._position[1] * self.a2 + p._position[2] * self.a3
   
   def get_patch_orientation_a1(self,patchid):
       p = self._patches[patchid]
       v = np.array( self.a1 * p._a1[0] + self.a2 * p._a1[1] + self.a3 * p._a1[2]  )
       return v
   
   def get_patch_orientation_a2(self,patchid):
       p = self._patches[patchid]
       if np.dot(p._a2,p._a2) < 0.9999:
           print "PRIMARY MAX CRITICAL ERRROR MAX CRITICAL ERROR", p._a2,  np.dot(p._a2,p._a2)
       v = np.array( self.a1 * p._a2[0] + self.a2 * p._a2[1] + self.a3 * p._a2[2]  )
       if np.dot(v,v) < 0.9999:
           print "MAX CRITICAL ERROR MAX CRITICAL ERROR", v, np.dot(v,v)
       return v
   
   def set_patch_a2_orientation(self,patchid,new_a2):
       coor = np.array([np.dot(new_a2,self.a1), np.dot(new_a2,self.a2), np.dot(new_a2,self.a3)])
       self._patches[patchid]._a2 = coor / np.sqrt(np.dot(coor,coor))
       
   def aligned_with(self,p):
      #checks if paticles are aligned
      print 'Verification of alignment ', self.cm_pos, p.cm_pos
      if len(self._patches) != len(p._patches):
           return None
      correspondence = {} 
      for i,patchA in enumerate(self._patches):
          positionA = self.get_patch_position(i)
          for j,patchB in enumerate(p._patches):
                 positionB = p.get_patch_position(j)
                
                 val = np.dot(positionA/l2norm(positionA),positionB/l2norm(positionB))
                 if val > 1.0 - myepsilon:
                     if j in correspondence.values():
                         #print 'Error two patches would correspond to the same patch'
                         return None
                     else:
                         correspondence[i] = j
                         print 'CHECKING patch positions, we have MATCH ' ,i,j, positionA, positionB, np.dot(positionA/l2norm(positionA),positionB/l2norm(positionB))
                         break
          if i not in correspondence.keys():
              print 'Could not match patch ',i  
              return None
       
      print 'Found perfect correspondence',correspondence    
      return correspondence
   
   
   def better_align_with(self,part2):
      all_pos = [x._position for x  in self._patches]
      all2_pos =  [x._position for x  in part2._patches]
      print 'Trying to align FROM:', all_pos , '\n TO: ', all2_pos   
      for pi,p in enumerate(self._patches):
           for qi,q in enumerate(self._patches):
               if qi != pi:
                   for li,l in enumerate(part2._patches):
                       for fi,f in enumerate(part2._patches):
                           if li != fi:
                               print 'Aligning patch %d with %d; and %d with %d' % (pi,li,qi,fi)
                               v1 = p._position / l2norm(p._position)
                               v2 = l._position / l2norm(l._position)
                               b1 = q._position / l2norm(q._position)
                               b2 = f._position / l2norm(f._position)
                               v1 = np.matrix(v1).transpose()
                               b1 = np.matrix(b1).transpose()
                               B = v1  * v2 + b1 * b2
                               U, s, V = np.linalg.svd(B, full_matrices=True)
                               M = np.diag([1, 1, np.linalg.det(U) * np.linalg.det(V)])
                               R = U * M * V
                               rot = np.asarray(R)
                               test = copy.deepcopy(part2)
                               test.rotate(rot)
                               c = self.aligned_with(test)
                               if c != None:
                                   print 'Success! '
                                   xxx = [test.get_patch_position(i) for i in xrange(len(test._patches))]
                                   print 'Using rotation \n', rot
                                   print 'After rotatoin patches change to ',xxx
                                   return c,rot
      print 'MAXIMUM ERROR'
      raise IOError('Cannot align patches')
      return None,[]
                    
                               
   def align_with(self,part2,rot_matrix=None):
       #this method tries to align particle with particle p, so that their patches overlap
       if len(self._patches) != len(part2._patches):
           return False
       
       all_pos = [x._position for x  in self._patches]
       all2_pos =  [x._position for x  in part2._patches]
       print 'Trying to align FROM:', all_pos , '\n TO: ', all2_pos   
       rot = copy.deepcopy(rot_matrix)
       if rot_matrix != None:
           np.transpose(rot)
           test = copy.deepcopy(part2)
           test.rotate(rot)
           xxx = [test.get_patch_position(i) for i in xrange(len(test._patches))]
           print 'Using rotation \n', rot_matrix
           print 'After rotatoin patches change to ',xxx
           c = self.aligned_with(test)
           if c != None:
                return c,rot
           else:
                print 'More detailed test of alignment failed'
                raise IOError ('Could not align')
       
       #TOHL JE BLBE!!    
       all_pos = [x._position for x  in self._patches]
       all2_pos =  [x._position for x  in part2._patches]
       print 'Trying to align ', all_pos , ' to ', all2_pos   
       print 'Of particles whose positions are', self.cm_pos, part2.cm_pos
       for pi,p in enumerate(self._patches):
           for qi,q in enumerate(self._patches):
               if qi != pi:
                   for li,l in enumerate(part2._patches):
                       for fi,f in enumerate(part2._patches):
                           if li != fi:
                               print 'Aligning patch %d with %d; and %d with %d' % (pi,li,qi,fi)
                               v1 = p._position / l2norm(p._position)
                               v2 = l._position / l2norm(l._position)
                               b1 = q._position / l2norm(q._position)
                               b2 = f._position / l2norm(f._position)
                               #print 'Positions are', v1,v2,b1,b2
                               #print 'NP.dot is ',np.dot(v1,v2)
                               theta = np.arccos(np.dot(v1,v2))
                               print 'Theta is' , theta
                               if abs(theta) < myepsilon:
                                   r1 = np.eye(3)
                               else:
                                   u1 = np.cross(v1,v2)
                                   #print 'U1 is',u1
                                   if l2norm(u1) < myepsilon: #this means u1 = -u2, we pick u1 as perpendicular to v1
                                          u1 = np.array([v1[1],-v1[0],0])
                                          if(l2norm(u1) == 0):
                                               u1 = np.array([0,-v1[1],v1[2]])
                                   u1 = u1/l2norm(u1)
                                   r1 = rotation_matrix(u1, theta)
                               v1 = np.dot(r1,v1)
                               b1 = np.dot(r1,b1)
                               #print r1
                               #print v1
                               print 'Po natoceni',np.dot(v1,v2)
                               b1proj = b1 - np.dot(b1,v1)
                               b2proj = b2 - np.dot(b2,v1)
                               b1proj = b1proj / l2norm(b1proj)
                               b2proj = b2proj / l2norm(b2proj)
                               print 'Dot and Theta2 is ',np.dot(b1proj,b2proj)
                               if np.dot(b1proj,b2proj) < 1 - myepsilon :
                                  theta2 = np.arccos(np.dot(b1proj,b2proj))
                                  u2 = v1
                                  u2 = u2/l2norm(u2)
                                  r2 = rotation_matrix(u2,theta2)
                                  r2PI = rotation_matrix(u2,2.*math.pi-theta2)
                               else:
                                  r2 = np.eye(3)
                                  r2PI = r2
                               v1 = np.dot(r2,v1)
                               b1old = b1
                               b1 = np.dot(r2,b1)
                               print 'After final alignment', np.dot(v1,v2), np.dot(b1,b2)
                               
                               xxx = copy.deepcopy(part2)
                               rot = np.dot(r1,r2)
                               print "Trying rotation matrix ", rot
                               np.transpose(rot)
                               xxx.rotate(rot)
                               print xxx.export_to_mgl(['blue','green','red','black','yellow','cyan','magenta','orange','violet'],'blue')
                               
                               if(  np.dot(b1,b2) < 1-myepsilon ):
                                   b1 = b1old
                                   r2 = r2PI
                                   b1 = np.dot(r2,b1)
                                   if np.dot(b1,b2) < 1-myepsilon:
                                       print 'Alignment double failed',np.dot(b1,b2)
                                       continue  
                                   else:
                                       print 'Second PI alignment was successful'
                               test = copy.deepcopy(part2)
                               rot = np.dot(r1,r2)
                               print "Trying rotation matrix ", rot
                               np.transpose(rot)
                               test.rotate(rot)
                               c = self.aligned_with(test)
                               if c != None:
                                   return c,rot
                               else:
                                   print 'More detailed test of alignment failed'
                               
        #if we got all the way here         
       return None,None            
                               
                               
           
   def distance_from(self,particle,box_size):
       d = self.cm_pos - particle.cm_pos
       d[0] -= np.rint(d[0] / float(box_size[0])) * box_size[0]
       d[1] -= np.rint(d[1] / float(box_size[1])) * box_size[1]
       d[2] -= np.rint(d[2] / float(box_size[2])) * box_size[2]
       return np.sqrt(np.dot(d,d))
                            
   def add_patch(self,patch):
       if(self._patches == None):
           self._patches = []
       self._patches.append(patch)
       

   def fill_patches(self,patch_array):
       if self._patch_ids != None:
           self._patches = []
           for i,id in enumerate(self._patch_ids):
               self._patches.append( copy.deepcopy(patch_array[id]))
    
   def fill_configuration(self,ls):
        self.cm_pos = np.array([float(x) for x in ls[0:3]])
        self.a1 = np.array([float(x) for x in ls[3:6]])
        self.a3 = np.array([float(x) for x in ls[6:9]])
        self.v = np.array([float(x) for x in ls[9:12]])
        self.L = np.array([float(x) for x in ls[12:15]])
   
   def save_type_to_string(self):
       outs = 'particle_%d = { \n type = %d \n ' % (self._type,self._type)
       outs = outs + 'patches = '
       for i,p in enumerate(self._patches):
           outs = outs + str(p._type) 
           if i < len(self._patches)-1:
               outs = outs +','
       outs = outs + ' \n } \n'
       return outs
   
   def save_conf_to_string(self):
       str = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' % (self.cm_pos[0],self.cm_pos[1],self.cm_pos[2],self.a1[0],self.a1[1],self.a1[2],self.a3[0],self.a3[1],self.a3[2],self.v[0],self.v[1],self.v[2],self.L[0],self.L[1],self.L[2] )       
       return str +'\n'
   
   def init_from_string(self,lines):
       for line in lines:
            line = line.strip()
            if len(line) > 1 and line[0] != '#':
               if "type" in line:
                   vals = int(line.split('=')[1])
                   self._type = vals
               if 'patches' in line:
                   vals = line.split('=')[1]
                   self._patch_ids = [int(g) for g in vals.split(',')]
                    
   def get_a2 (self):
        return np.cross (self.a3, self.a1)

   a2 = property (get_a2)       
  
    
   def export_to_mgl(self,patch_colors,particle_color,patch_width=0.1,patch_extension=0.2): 
       #sout = '%f %f %f @ %f C[%s] ' % (self.cm_pos[0],self.cm_pos[1],self.cm_pos[2],self._radius,particle_color)

       sout = '%f %f %f @ %f C[%s] ' % (self.cm_pos[0],self.cm_pos[1],self.cm_pos[2],0.4,particle_color)
       if len(self._patches) > 0:
           sout = sout + 'M '
       for i,p in enumerate(self._patches):
           pos = p._position[0] * self.a1 + p._position[1] * self.a2 + p._position[2] * self.a3
           pos *= (1.0+patch_extension) 
           #print 'Critical: ',p._type,patch_colors
           g = '%f %f %f %f C[%s] ' %(pos[0],pos[1],pos[2],patch_width,patch_colors[i])  
           sout = sout + g
       return sout

   def export_to_xyz(self,patch_width=0.4,patch_extension=0.2): 
       letter = PLPatchyParticle.all_letters[self._type]
       sout = '%s %f %f %f ' % (letter,self.cm_pos[0],self.cm_pos[1],self.cm_pos[2])
       
       return sout
   
   
   
class PLPSimulation:
    def __init__(self,seed=88):
        self._topology_file = ''
        self._config_file = ''
        self._input_file = ''
        self._box_size = np.array([0.,0.,0.])
        self.N = 0
        self._N_particle_types = 0
        self._N_patch_types = 0
        self._patch_types = []
        self._particle_types = []
        self._E_tot = 0.
        self._E_pot = 0.
        self._E_kin = 0.
        
        self._particle_colors = ['blue','green','red','black','yellow','cyan','magenta','orange','violet']
        self._patch_colors = ['green','violet','pink','brown','orange','red','black']
        self._complementary_colors = {20: 'blue', -20: 'cyan', 21: 'red', -21: 'yellow',22: 'black', -22: 'brown' }
        
        
        self.particles = []
        random.seed(seed)
        
        self._colorbank = []


    def set_radius(self,radius):
        for p in self.particles:
               p.set_radius(radius)
    
    def generate_random_color(self,id=-1):
        #random.seed(seed)
        if(True or len(self._colorbank) < 1 or id == -1):
          r = [random.random() for i in xrange(3)]
          color = '%f,%f,%f' % (r[0],r[1],r[2])
          return color
        else:
          return self._colorbank[id]
        
    def translate(self,dir):
        for p in self.particles:
            p.translate(dir)
            
    def set_box_size(self,box):
        self._box_size = np.array(box)
        
    def load_input_file(self,input_file):
        patchy_file = None
        particle_file = None
        
        handle = open(input_file,'r')
        for line in handle.readlines():
            line = line.strip()
            if(len(line) > 1 and line[0] != '#'):
                if 'patchy_file' in line:
                  val = line.split('=')[1].strip()
                  patchy_file = val
                elif 'particle_file' in line:
                  val = line.split('=')[1].strip()
                  particle_file = val
                elif 'particle_types_N' in line:
                  val = line.split('=')[1].strip()
                  self._N_particle_types = int(val)
                elif 'patch_types_N' in line:
                  val = line.split('=')[1].strip()
                  self._N_patch_types = int(val)  
        
        handle.close()
        self._colorbank = [self.generate_random_color() for x in range(self._N_patch_types)]
        #print >> sys.stderr, "Loaded patchy_file: %s, particle_file: %s, N_part_types: %d, n_patch_types: %d " % (patchy_file,particle_file,self._N_particle_types, self._N_patch_types)
        #now process the patch types file:
        patches = [None for x in range(self._N_patch_types)]
        handlep = open(patchy_file,'r')
        j = 0
        Np = 0
        lines = handlep.readlines()
        for line in lines:
            line = line.strip()
            if len(line) > 1 and line[0] != '#':
                if 'patch_' and '{' in line:
                    strargs = []
                    k = j+1
                    while '}' not in lines[k]:
                        strargs.append(lines[k].strip())
                        k = k + 1
                    patch = Patch()   
                    #print 'Loaded patch',strargs
                    patch.init_from_string(strargs)
                    index = patch._type
                    patches[index]= patch    
                    Np += 1
            j = j + 1
            
            
        if Np != self._N_patch_types:
             raise IOError('Loaded %d patches, as opposed to the desired %d types ' % (Np, self._N_patch_types))
        
        self._patch_types = patches
            
        particles = [None for x in range(self._N_particle_types)]    
        handlep.close()
        handlep = open(particle_file,'r')
        lines = handlep.readlines()
        j = 0
        for line in lines:
            line = line.strip()
            if len(line) > 1 and line[0] != '#':
                if 'particle_' and '{' in line:
                    strargs = []
                    k = j+1
                    while '}' not in lines[k]:
                        strargs.append(lines[k].strip())
                        k = k + 1
                    particle = PLPatchyParticle()    
                    #print 'Loaded particle ',strargs
                    particle.init_from_string(strargs)
                    particle.fill_patches(self._patch_types)
                    index = particle._type
                    particles[index] = copy.deepcopy(particle)    
                    Np += 1
            j = j + 1
        self._particle_types = particles
        #print 'Critical', self._particle_types
        #print "LOADED", len(self._particle_types)
    
    def load_topology(self,topology_file):
        handle = open(topology_file,'r')
        line = handle.readlines()
        vals = line[0].split()
        N = int(vals[0])
        Ntypes = int(vals[1])
        types = [int(x) for x in line[1].split()]
        #print 'Critical', line[1].split()
        #print types
        #print self._particle_types
        #print 'THERE ARE', len(self._particle_types), ' particle types '
        self.particles = []
        for index,type in enumerate(types):
            p = copy.deepcopy(self._particle_types[type])
            p.index = index
            self.particles.append(p)
    
        handle.close()
        self.N = N
        if N != len(self.particles):
            raise IOError('Particle number mismatch while reading topology file')
        
    def save_topology(self,topology_file):
        handle = open(topology_file,'w')
        handle.write('%d %d\n' %( len(self.particles), self._N_particle_types))
        outstr = ''
        for p in self.particles:
            outstr = outstr +  str(p._type) + ' '
        handle.write(outstr + '\n')
        handle.close()
        
    def save_patchy_types_file(self,ptypes_file): 
        handle = open(ptypes_file,'w')
        for p in self._patch_types:
            outs = p.save_to_string() 
            handle.write(outs)   
        handle.close()
    
    def save_particle_types_file(self,ptypes_file):
        handle = open(ptypes_file,'w')
        for p in self._particle_types:
            outs = p.save_type_to_string() 
            handle.write(outs)   
        handle.close()
            
    
    def check_for_particle_overlap(self,particle,dist_cutoff=1.0):  
        #print 'Adding ', particle.cm_pos      
        for p in self.particles:
            dist = p.distance_from(particle,self._box_size)
            #print ' Looking at distance from ', p.cm_pos,dist
            if dist <= dist_cutoff:
                #print 'Unfortunately overlaps with ',p.cm_pos,dist
                return True    
        #print 'Check is fine!'    
        return False
                      
    def save_configuration(self,conf_name,t=0.):
        handle = open(conf_name,'w')
        handle.write('t = %f\nb = %f %f %f\nE = %f %f %f\n' % (t,self._box_size[0],self._box_size[1],self._box_size[2],self._E_pot,self._E_kin,self._E_tot))
        for p in self.particles:
            outs = p.save_conf_to_string()
            handle.write(outs)
        handle.close()    
            
    def add_particles(self,particles,strict_check=True):
        #adds particles to the field, also initializes paricle types and patchy types based on these data
        #it overwrites any previosuly stored particle!!
        self.particles = copy.deepcopy(particles)
        self.N = len(particles)
        #now treat types:
        saved_types = {}
        for p in self.particles:
            if not p._type in saved_types.keys():
                saved_types[p._type] = copy.deepcopy(p)
        self._particle_types = []
        for i,key_id in enumerate(sorted(saved_types.keys())):
            if key_id != i and strict_check:
                raise IOError("Error while adding particles to the PLPSimulation class, indices of types are not correctly ordered")
            self._particle_types.append(copy.deepcopy(saved_types[key_id]))
        self._N_particle_types = len(self._particle_types)    
        
        #now treat patches    
        saved_patch_types = {}
        for p in self._particle_types:
            for patch in p._patches:
                if not patch._type in saved_patch_types.keys():
                    saved_patch_types[patch._type] = copy.deepcopy(patch)
        self._patch_types = []
        for i,key_id in enumerate(sorted(saved_patch_types.keys())):
            if key_id != i and strict_check:
                raise IOError("Error while adding patches to the PLPSimulation class, indices of types are not correctly ordered")
            self._patch_types.append(copy.deepcopy(saved_patch_types[key_id]))            
        self._N_patch_types = len(self._patch_types)
                    
    def insert_particle(self,particle,check_overlap=False):
        if check_overlap: 
         if  check_for_particle_overlap(particle) == True:
             return False
        self.particles.append(particle)
        self.N += 1
        if particle._type not in [x._type for x in self._particle_types]:
             self._particle_types.append(copy.deepcopy(particle))
             self._N_particle_types += 1
        return True
        
    def load_configuration(self,configuration_file,conf_to_skip=0,close_file=True):
        _conf = open(configuration_file,'r')
        if conf_to_skip > 0:
            conf_lines = 3 + self.N
            for j in xrange(conf_lines*conf_to_skip):
                _conf.readline()
                
        self.read_next_configuration(_conf)
        
        if close_file:
          _conf.close()   
        
        return _conf
       
    def load_from_files(self,input_file,topology_file,config_file,conf_to_skip=0):
        self.load_input_file(input_file)
        self.load_topology(topology_file)
        self.load_configuration(config_file,conf_to_skip)
        
    def read_next_configuration(self,file_handle):
        _conf = file_handle
        
        timeline = _conf.readline()
        time = 0.
        if  len(timeline) == 0:
            return False
        else:
            time = float(timeline.split()[2])

        box = np.array([float(x) for x in _conf.readline().split()[2:]])
        [E_tot, E_pot, E_kin] = [float(x) for x in _conf.readline().split()[2:5]]
        
        self._box_size = box
        self._E_tot = E_tot
        self._E_pot = E_pot
        self._E_kin = E_kin
        
        
        for i in xrange(self.N):
            ls = _conf.readline().split()
            self.particles[i].fill_configuration(ls)   
            
        return _conf     
        
        
    def bring_in_box(self,all_positive=False):
        for p in self.particles:
            nx = np.rint(p.cm_pos[0] / float(self._box_size[0])) * self._box_size[0]
            ny = np.rint(p.cm_pos[1] / float(self._box_size[1])) * self._box_size[1]
            nz = np.rint(p.cm_pos[2] / float(self._box_size[2])) * self._box_size[2]
            #print np.array([nx,ny,nz])
            p.cm_pos -= np.array([nx,ny,nz])
            if all_positive:
                for i in xrange(3):
                    if p.cm_pos[i] < 0:
                        p.cm_pos[i] += self._box_size[i]          
    
    def _get_color(self,index,ispatch=False):
        if abs(index) >= 20:
            index = abs(index)-20
        #print index
        #print self._patch_colors
        return self._patch_colors[index]


        if not ispatch:
           if index <0 or index  >= len(self._particle_colors):
               if index in self._complementary_colors.keys():
                   return self._complementary_colors[index]
               else:
                   return self.generate_random_color()
           else:
               return self._particle_colors[index]
        else:
           if index < 0 or index >= len(self._patch_colors):
               if index in self._complementary_colors.keys():
                   return self._complementary_colors[index]
               else:
                   return self.generate_random_color(index)
           else:
               return self._patch_colors[index]
   
        
    def export_to_mgl(self,filename,regime='w',icosahedron=True):
        out = open(filename,regime)
        sout = ".Box:%f,%f,%f\n" % (self._box_size[0],self._box_size[1],self._box_size[2])
        for p in self.particles:
            patch_colors = [self._get_color(pat._color,True) for pat in p._patches]
            particle_color = self._get_color(p._type)
            sout = sout + p.export_to_mgl(patch_colors,particle_color) + '\n'
            if icosahedron:
                 line = "%f %f %f @ 0.5 C[%s] I %f %f %f %f %f %f \n" % (p.cm_pos[0],p.cm_pos[1],p.cm_pos[2],particle_color , p.a1[0],p.a1[1],p.a1[2], p.a2[0],p.a2[1],p.a2[2] )
                 sout = sout + line
        
        out.write(sout)
        out.close()    
                                                                                   
    def export_to_xyz(self,filename,regime='w'):
        out = open(filename,regime)
          
        sout = str(len(self.particles))+'\n'
        sout += "Box:%f,%f,%f\n" % (self._box_size[0],self._box_size[1],self._box_size[2])
        for p in self.particles:
            sout = sout + p.export_to_xyz() + '\n'
        
        out.write(sout)
        out.close()    
                                                                                   
                     

    


