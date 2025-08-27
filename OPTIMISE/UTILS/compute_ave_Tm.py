#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:55:44 2024

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

bps = 15
Ct = 0.000672
Cs = 0.5
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
        
"""
SEQS = []

seq = np.zeros(bps,dtype=int)

#generate all sequences
loop_seqs(seq,bps,SEQS)


#sort sequences in ascending melting temperature
SEQS.sort(key=lambda x: x.mT)

print("max mT: "+ SEQS[len(SEQS)-1].seq + " " +str(SEQS[len(SEQS)-1].mT))
print("min mT: "+ SEQS[0].seq + " "+ str(SEQS[0].mT))
"""

SEQ = ""

for i in range(bps) :
    SEQ += "C"


ave_mT =SL.melting_temperature_ave(SEQ,Ct,Cs)
    
print("ave mT: "+str(ave_mT))


