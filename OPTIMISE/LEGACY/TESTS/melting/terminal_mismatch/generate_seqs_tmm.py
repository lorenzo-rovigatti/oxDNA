import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import sys
import random

if len(sys.argv) != 3:
    print("Syntax")
    print(sys.argv[0] + " seq0 seq1")
    print("seq1-2: 2 sequences (same length) starting with AT, CG respectively")
    exit(1)

def get_mismatches(base):
    mms = []
    if base == 'A':
        mms = ['A', 'G', 'C']
    elif base == 'G':
        mms = ['A', 'G', 'T']
    elif base == 'C':
        mms = ['A', 'C', 'T']
    else :
        mms = ['G', 'C', 'T']
    return mms

Ct=0.000336*1000000000 #in nM
Cs=500 #0.5M

bases = ['A', 'G', 'C', 'T']
c_bases = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}

def get_random_base():
    random_num = random.randint(0, 3)
    return bases[random_num]

def get_mT(seq, seq_mm, comp) :
    mT = mt.Tm_NN(seq, Na=500, dnac1=Ct, dnac2=Ct, nn_table=mt.DNA_NN4) #SantaLucia duplex
    mT_mm = mt.Tm_NN(seq_mm, c_seq=comp, Na=500, dnac1=Ct, dnac2=Ct, nn_table=mt.DNA_NN4, tmm_table=mt.DNA_TMM1) #mismatches SantaLucia
    
    return mT, mT_mm

print("Ct (single strand) = 0.000336 M")
print("Salt concentration (Na) = 0.5M")

seq0 = []

seq0.append(sys.argv[1])
seq0.append(sys.argv[2])

Nbps = len(seq0[0])

if len(seq0[0]) != len(seq0[1]) :
    print("sequences are not of the same length")
    exit(1)

seqs = []
seqs_mm = []
c_seqs_mm = []
mTs = []
mTs_mm = []
mm_ids = []

for i in range(4):
    for j in range(2) :
        b1 = bases[i]
        seq = seq0[j]
        seq_mm = b1+seq
        
        mms = get_mismatches(b1)
        #print(mms)
        
        for z in range(len(mms)):
            seqs.append(seq)
            seqs_mm.append(seq_mm)
            mm_ids.append(z)
            #reverse because oxdna reads 3'->5', whereas biopython 5'->3'
            seq_r = seq[::-1]
            seq_mm_r = seq_mm[::-1]
            c_seq_r=""
            for k in range(0,len(seq_mm_r)-1): c_seq_r+=c_bases[seq_mm_r[k]]
            c_seq_r+=mms[z]
            c_seq=c_seq_r[::-1]
            c_seqs_mm.append(c_seq)
            #print(seq_r, seq_mm_r, c_seq_r)
            mT, mT_mm = get_mT(seq_r, seq_mm_r, c_seq_r)
            mTs.append(mT)
            mTs_mm.append(mT_mm)
        

#generate op file
op_file_name = "op_n"+str(Nbps)+".txt"
ofile=open(op_file_name,'w')
print("{",file=ofile)
print("order_parameter = bond",file=ofile)
print("name = all_native_bonds",file=ofile)

for i in range(1,Nbps+1) :
    line = "pair"+str(i)+" = "+str(i)+" "+str(2*Nbps+1-i)
    print(line,file=ofile)
print("}",file=ofile)
ofile.close()

#print strings for prepare_sims

seqs_str = "seqs=( "
op_str = "OP_files=( "
mTs_str = "mTs=( "
nbps_str = "nbps=( "
mmid_str = "mm_id=( "
        
for i in range(len(seqs_mm)) :
    op_str+=op_file_name+" "
    nbps_str+=str(Nbps)+" "
    mTs_str+=str(mTs_mm[i])[:5]+" "
    seqs_str+=seqs_mm[i]+" "
    mmid_str+=str(mm_ids[i])+" "
    
op_str+=")"
nbps_str+=")"
mTs_str+=")"
seqs_str+=")"
mmid_str+=")"

ofile = open("prepare_term_mm_headers.txt", 'w')
print(seqs_str,file=ofile)
print(op_str,file=ofile)
print(nbps_str,file=ofile)
print(mTs_str,file=ofile)
print(mmid_str,file=ofile)

ofile.close()

ofile = open("out_gen_term_mm.txt", 'w')

print("#seq_mm c_seq_mm mT_mm mT", file=ofile)
delta_T=""
for i in range(len(seqs_mm)) :
    line=seqs_mm[i] + " " + c_seqs_mm[i] + " " + str(mTs_mm[i])[:5] + " " + str(mTs[i])[:5]
    print(line, file=ofile)
    delta_T+=str(mTs_mm[i]-mTs[i])[:5]+" "
print("DT: ",delta_T)
ofile.close()

