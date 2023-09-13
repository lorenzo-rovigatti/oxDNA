import os
import sys
import numpy as np

if len(sys.argv) not in [5]:
	print("\033[0;31mUsage is %s omegas e3 step a \033[0m" % sys.argv[0])
	sys.exit()


Omega3Binary = sys.argv[1]
f_1 = os.path.splitext(sys.argv[1])
e3Binary = sys.argv[2]
step = int(sys.argv[3])
a = np.float(sys.argv[4])

#loading binary files
omega3 = np.load(Omega3Binary)
e3 = np.load(e3Binary)
shape = np.shape(omega3)
#print shape
s = shape[0]
N_bp   = shape[1]

#l_t/2
for m in range(1,step+1):
        om3 = np.zeros((s,N_bp-3-m-20),dtype=np.float64)
        for k in range (s):    
            for i in range (10,N_bp-3-m-10):
                om3[k][i-10] = np.sum(omega3[k][i:i+m])
        cosomega = np.cos(a*om3)                            
        l_t_half = (-0.5*m*a)/(np.log(np.mean(cosomega)))   # use exact value of 'a' here, calulated in triad3.py stored in file 'a'
        with open (f_1[0] +'_lt' , "a") as f:
            print >> f ,  '{:^8}{:^8}'.format(m,round(l_t_half,1))  
 
#lb
for m in range(1,step+1):
        costheta = np.zeros((s,N_bp-2-m-20),dtype=np.float64)
        for k in range (s):         
            for i in range (10,N_bp-2-m-10):                  # 10 bps from start and end of strand are not taken into account for analysis as they show splaying
                costheta[k][i-10] = np.dot(e3[k][i],e3[k][i+m])        
        l_b = (-m*a)/np.log(np.mean(costheta))             # use exact value of a here, calulated in triad3.py stored in file 'a'
        with open (f_1[0] + '_lb' , "a") as f:
            print >> f ,  '{:^8}{:^8}'.format(m,round(l_b,1))
        
#e3
for m in range(1,step+1):
        costheta = np.zeros((s,N_bp-2-m-20),dtype=np.float64)
        for k in range (s):
            for i in range (10,N_bp-2-m-10):
                costheta[k][i-10] = np.dot(e3[k][i],e3[k][i+m])
        cos_theta = (np.mean(costheta))
        with open (f_1[0] + '_cos_theta' , "a") as f:
            print >> f ,  '{:^8}{:^8}'.format(m,(cos_theta))
