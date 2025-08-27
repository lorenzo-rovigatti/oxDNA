#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:15:12 2024

@author: yqb22156
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 12:43:04 2024

@author: yqb22156

generates N sequences with bps base pairs in such a way that the melting temperature is uniformly reprsented


"""
import sys
import numpy as np
import random
import matplotlib.pyplot as plt
import SantaLucia as SL

"""
if len(sys.argv) != 6 :
    print("Unknown argument format.")
    print("Usage: python3 gen_seq.py bps Ct Cs N gen_type")
    print("bps = number of base pairs")
    print("Ct = total single strand concentration in M")
    print("Cs = salt concentration in M")
    print("For 1 duplex in a box of size l ox units, Ct = 2/l^3*2.6868 M")
    print("N = number of sequebces to generate")
    print("gen_type = 0 for accounting for conf space population vs mT")
    print("gen_type = anything else for uniformly spaced mTs")
    
    sys.exit()
    
bps = int(sys.argv[1])
Ct = float(sys.argv[2])
Cs = float(sys.argv[3])
N = int(sys.argv[4])
uni_T = True
if sys.argv[4] == 0:
    uni_T = False
    
"""
Ngen = 1000000
N = 40
bps = 15
Ct = 0.000672
Cs = 0.15
uni_T = True

class Sequence :
    
    def __init__(self,seq_int) :
        
        SEQ = ""
        
        for i in range(len(seq_int)) :
            SEQ += SL.id_to_base(seq_int[i])
            
        MT = SL.melting_temperature(SEQ,Ct,Cs)
        
        self.seq = SEQ
        self.mT = MT
        

def loop_seqs(seq,n,SEQS) :
    if n>0 :
        for i in range(4) :
           seq[n-1] = i
           loop_seqs(seq,n-1,SEQS)
    else:        
        SEQS.append(Sequence(seq))
        
        
def search_complementary(s,seqs) :
    
    s_rev= ""                                                   
    s_rev = "".join(reversed(s.seq))

    for i in sampled :
        if i.seq == s_rev :
            return True
        
    return False

SEQS = []

seq_min = []
seq_max = []

#find min max

for i in range(bps) :
    if i%2 == 0:
        seq_min.append(0)
        seq_max.append(1)
    else :
        seq_min.append(3)
        seq_max.append(2)
        
SMin = Sequence(seq_min)
SMax = Sequence(seq_max)



#generate Ngen random sequences

for i in range(Ngen) :
    seq = []
    for j in range(bps):
        seq.append(random.randint(0,3))
        
    SEQS.append(Sequence(seq))



#sort sequences in ascending melting temperature
SEQS.sort(key=lambda x: x.mT)

print("max mT: "+ SMax.seq + " " + str(SMax.mT))
print("min mT: "+ SMin.seq + " "+ str(SMin.mT))


#crate histogram

# Generate some data
data = []

for sq in SEQS :
    data.append(sq.mT)

# Create histogram with 30 bins and a black edge color
plt.hist(data, bins=30, color='skyblue', edgecolor='black')

# Add labels and title
plt.xlabel('melting T [C˚]')
plt.ylabel('Frequency')
plt.title('8mers melting T')

# Display the plot
plt.show()




delta = int(len(SEQS)/N)

sampled = []

if uni_T :  
    
    #uniform mT sampling  
    
    delta = (SEQS[len(SEQS)-1].mT - SEQS[0].mT)/N
    
    k0 = 0
    
    for n in range(N) :
        
        #print(n)
        
        tmT = SEQS[0].mT+(n+1)*delta
        
        #print(tmT)
        
        for k in range(k0,len(SEQS)-1) :
                
            if SEQS[k].mT <= tmT and (SEQS[k+1].mT >= tmT or k == len(SEQS)-2):
                
                r = random.randint(k0, k)
                
                #print(r)
                
                sampled.append(SEQS[r])
                
                k0 = k+1

else:
    
    #sampling accounting for population
    
    for n in range(N) :
        
        while 1 == 1:
            r = n*delta + random.randint(0, delta)
            
            if SEQS[r] in sampled or search_complementary(SEQS[r],sampled) :
                continue
            else :
                sampled.append(SEQS[r])
                break
            
            
ofile = open("gen_seqs_n"+str(bps)+".txt", 'w')

for seq in sampled :
    print(seq.seq + " " + str(seq.mT), file = ofile)

ofile.close()
