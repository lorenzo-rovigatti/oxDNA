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

N = 100
bps = 5
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
        


SEQS = []
def loop_seqs(seq,n) :
    global SEQS
    if n>0 :
        for i in range(4) :
           seq[n-1] = i
           loop_seqs(seq,n-1)
    else:        
        SEQS.append(Sequence(seq))


def loop_seqs_cutoff(seq,n) :
    global SEQS
    if n>0 :
        for i in range(4) :
           seq[n-1] = i
           loop_seqs_cutoff(seq,n-1)
    else:
        tmp = Sequence(seq)
        if tmp.mT > 12:
            #print(tmp.mT)
            SEQS.append(tmp)


def search_reverse(s,seqs) :

    s_rev= ""
    s_rev = "".join(reversed(s.seq))

    for i in sampled :
        if i.seq == s_rev :
            return True

    return False

comple = {
 'A':'T',
 'C':'G',
 'G':'C',
 'T':'A'
}

def search_complementary(s) :
    s_rev= ""
    for i in range(len(s.seq)-1,-1,-1):
        s_rev+=comple[s.seq[i]]

    for i in sampled :
        if i.seq == s_rev :
            return True

    return False

seq = np.zeros(bps,dtype=int)

#generate all sequences
loop_seqs_cutoff(seq,bps)

print("using cutoff 12 C.")

#sort sequences in ascending melting temperature
SEQS.sort(key=lambda x: x.mT)

print("max mT: "+ SEQS[len(SEQS)-1].seq + " " +str(SEQS[len(SEQS)-1].mT))
print("min mT: "+ SEQS[0].seq + " "+ str(SEQS[0].mT))


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

coarse = 50

if uni_T :  
    
    #uniform mT sampling  
    
    delta = coarse*(SEQS[len(SEQS)-1].mT - SEQS[0].mT)/N
    
    k0 = 0
    
    for n in range(int(N/coarse)) :
        
        print(n)
        
        tmT = SEQS[0].mT+(n+1)*delta
        
        #print(tmT)
        
        for k in range(k0,len(SEQS)-1) :
            print(k0,k)
            if SEQS[k].mT <= tmT and (SEQS[k+1].mT >= tmT or k == len(SEQS)-2):
                
                sampled_ids = []
                while len(sampled_ids) < coarse:
                    r = random.randint(k0, k)
                    if r in sampled_ids or search_complementary(SEQS[r]):
                        continue
                    else:
                        sampled_ids.append(r)
                        sampled.append(SEQS[r])
                #print(r)
                
                #sampled.append(SEQS[r])
                
                k0 = k
                break

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
            
sampled.sort(key=lambda x: x.mT)
            
ofile = open("gen_seqs_n"+str(bps)+".txt", 'w')

for seq in sampled :
    print(seq.seq + " " + str(seq.mT), file = ofile)

ofile.close()

ofile = open("for_melting_script.txt", 'w')
line1 = "seqs=( "
line2 = "OP_files=( "
line3 = "W_files=( "
line4 = "nbps=( "
line5 = "mTs=( "
for seq in sampled:
    line1+=seq.seq+" "
    line2+="op_n5.txt "
    line3+="wfile_n5.txt "
    line4+="5 "
    line5+=str(round(seq.mT))+" "
line1+=")"
line2+=")"
line3+=")"
line4+=")"
line5+=")"
print(line1,file=ofile)
print(line2,file=ofile)
print(line3,file=ofile)
print(line4,file=ofile)
print(line5,file=ofile)

ofile.close()
