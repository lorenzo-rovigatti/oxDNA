import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import sys
import random

if len(sys.argv) != 5:
    print("Syntax")
    print(sys.argv[0] + " seq0 seq1 seq2 seq3")
    print("seq1-4: 4 sequences (sme length) starting with A,C,G,T respectively")
    exit(1)

Ct=0.000336*1000000000 #in nM
Cs=500 #0.5M

bases = ['A', 'G', 'C', 'T']
c_bases = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}

def get_random_base():
    random_num = random.randint(0, 3)
    return bases[random_num]

def get_mT(seq, seq_de, comp) :
    mT = mt.Tm_NN(seq, Na=500, dnac1=Ct, dnac2=Ct, nn_table=mt.DNA_NN4) #SantaLucia duplex
    mT_de = mt.Tm_NN(seq_de, c_seq=comp, Na=500, dnac1=Ct, dnac2=Ct, nn_table=mt.DNA_NN4, de_table=mt.DNA_DE1) #dangling end NN4=SantaLucia
    
    return mT, mT_de

print("Ct (single strand) = 0.000336 M")
print("Salt concentration (Na) = 0.5M")

seq0 = []

seq0.append(sys.argv[1])
seq0.append(sys.argv[2])
seq0.append(sys.argv[3])
seq0.append(sys.argv[4])

Nbps = len(seq0[0])

if len(seq0[0]) != len(seq0[1]) or len(seq0[0]) != len(seq0[2]) or len(seq0[0]) != len(seq0[3]) :
    print("sequences are not of the same length")
    exit(1)

seqs = []
seqs_de = []
mTs = []
mTs_de = []

for i in range(4) :
    for j in range(4) :
        b1 = bases[i]
        seq = seq0[j]
        seq_de = b1+seq
     
        seqs.append(seq)
        seqs_de.append(seq_de)
        #reverse because oxdna reads 3'->5', whereas biopython 5'->3'
        seq_r = seq[::-1]
        seq_de_r = seq_de[::-1]
        c_seq_r=""
        for z in range(0,len(seq_de_r)-1): c_seq_r+=c_bases[seq_de_r[z]]
        #print(seq_r, seq_de_r, c_seq_r)
        mT, mT_de = get_mT(seq_r, seq_de_r, c_seq_r)
        mTs.append(mT)
        mTs_de.append(mT_de)
        

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
        
for i in range(len(seqs_de)) :
    op_str+=op_file_name+" "
    nbps_str+=str(Nbps)+" "
    mTs_str+=str(mTs_de[i])[:5]+" "
    seqs_str+=seqs_de[i]+" "
    
op_str+=")"
nbps_str+=")"
mTs_str+=")"
seqs_str+=")"

ofile = open("prepare_dangling_headers.txt", 'w')
print(seqs_str,file=ofile)
print(op_str,file=ofile)
print(nbps_str,file=ofile)
print(mTs_str,file=ofile)

ofile.close()

ofile = open("out_gen_dangling.txt", 'w')

print("#seq_de seq mT_de mT", file=ofile)
delta_T=""
for i in range(len(seqs_de)) :
    line=seqs_de[i] + " " + seqs[i] + " " + str(mTs_de[i])[:5] + " " + str(mTs[i])[:5]
    print(line, file=ofile)
    delta_T+=str(mTs_de[i]-mTs[i])[:5]+" "
print(delta_T)
ofile.close()

