import numpy as np
import sys

bases = ['A','C','G','T']

bases_id = {'A': 0, 'G': 1, 'C': 2, 'T': 3}

complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def id_to_tetra(TY) :
    ty0 = TY%4
    ty1 = (TY//4)%4
    ty2 = (TY//4//4)%4
    ty3 = (TY//4//4//4)%4

    return bases[ty0]+bases[ty1]+bases[ty2]+bases[ty3]

def id_to_bp(TY) :
    ty0 = TY%4
    ty1 = (TY//4)%4

    return bases[ty0]+bases[ty1]


if len(sys.argv)!=4:
    print("usage:", sys.argv[0], "file_n5_seqs file_n8_seqs file_n15_seqs")
    exit(1)

base_pairs = np.zeros(16,dtype=int)
trimers = np.zeros(64,dtype=int)
tetramers = np.zeros(256,dtype=int)

def bp_to_id(bp):
    return bases_id[bp[0]]+4*bases_id[bp[1]]

def tri_to_id(tri):
    return bases_id[tri[0]]+4*bases_id[tri[1]]+16*bases_id[tri[2]]

def tetra_to_id(tetra):
    return bases_id[tetra[0]]+4*bases_id[tetra[1]]+16*bases_id[tetra[2]]+64*bases_id[tetra[3]]

def reverse(seq):
    return ''.join(complement[base] for base in reversed(seq))


file_n5=open(sys.argv[1], 'r')

counts_bp = 0
counts_tetra = 0

for line in file_n5.readlines():
    seq=line.strip().split(" ")[0]
    seq_rev=reverse(seq)

    for i in range(len(seq)-1):
        bp=seq[i]+seq[i+1]
        base_pairs[bp_to_id(bp)]+=1
        bp=seq_rev[i]+seq_rev[i+1]
        base_pairs[bp_to_id(bp)]+=1
        counts_bp+=2

    for i in range(len(seq)-3):
        tetra=seq[i]+seq[i+1]+seq[i+2]+seq[i+3]
        tetramers[tetra_to_id(tetra)]+=1
        tetra=seq_rev[i]+seq_rev[i+1]+seq_rev[i+2]+seq_rev[i+3]
        tetramers[tetra_to_id(tetra)]+=1
        counts_tetra+=2

print("Counts bps n5:",counts_bp)
print("Counts tetra n5:",counts_tetra)
file_n5.close()

file_n8=open(sys.argv[2], 'r')

counts_bp = 0
counts_tetra = 0

for line in file_n8.readlines():
    seq=line.strip().split(" ")[0]
    seq_rev=reverse(seq)

    for i in range(len(seq)-1):
        bp=seq[i]+seq[i+1]
        base_pairs[bp_to_id(bp)]+=1
        bp=seq_rev[i]+seq_rev[i+1]
        base_pairs[bp_to_id(bp)]+=1
        counts_bp+=2

    for i in range(len(seq)-3):
        tetra=seq[i]+seq[i+1]+seq[i+2]+seq[i+3]
        tetramers[tetra_to_id(tetra)]+=1
        tetra=seq_rev[i]+seq_rev[i+1]+seq_rev[i+2]+seq_rev[i+3]
        tetramers[tetra_to_id(tetra)]+=1
        counts_tetra+=2

file_n8.close()

print("Counts bps n8:",counts_bp)
print("Counts tetra n8:",counts_tetra)


counts_bp = 0
counts_tetra = 0

file_n15=open(sys.argv[3], 'r')

for line in file_n15.readlines():
    seq=line.strip().split(" ")[0]
    seq_rev=reverse(seq)

    for i in range(len(seq)-1):
        bp=seq[i]+seq[i+1]
        base_pairs[bp_to_id(bp)]+=1
        bp=seq_rev[i]+seq_rev[i+1]
        base_pairs[bp_to_id(bp)]+=1
        counts_bp+=2

    for i in range(len(seq)-3):
        tetra=seq[i]+seq[i+1]+seq[i+2]+seq[i+3]
        tetramers[tetra_to_id(tetra)]+=1
        tetra=seq_rev[i]+seq_rev[i+1]+seq_rev[i+2]+seq_rev[i+3]
        tetramers[tetra_to_id(tetra)]+=1
        counts_tetra+=2

print("Counts bps n15:",counts_bp)
print("Counts tetra n15:",counts_tetra)

file_n15.close()

ofile=open("base_pairs_counts.txt",'w')
for i in range(len(base_pairs)):
    print(id_to_bp(i),base_pairs[i],file=ofile)
ofile.close()

ofile=open("tetramers_counts.txt",'w')
for i in range(len(tetramers)):
    print(id_to_tetra(i),tetramers[i],file=ofile)
ofile.close()
