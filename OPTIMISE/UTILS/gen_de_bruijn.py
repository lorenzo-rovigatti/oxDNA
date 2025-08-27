import random
import copy
import SantaLucia as SL

bases = ['A','C','G','T']

bases_id = {'A': 0, 'G': 1, 'C': 2, 'T': 3}

complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

N_n5 = 130
N_n8 = 16
N_n15 = 16

Ct = 0.000672
Cs = 0.15

class Sequence :
    
    def __init__(self,SEQ) :
            
        MT = SL.melting_temperature(SEQ,Ct,Cs)
        
        self.seq = SEQ
        self.mT = MT

def id_to_tetra(TY) :
    ty0 = TY%4
    ty1 = (TY//4)%4
    ty2 = (TY//4//4)%4
    ty3 = (TY//4//4//4)%4

    return bases[ty0]+bases[ty1]+bases[ty2]+bases[ty3]

def tetra_to_id(tetra):
    return bases_id[tetra[0]]+4*bases_id[tetra[1]]+16*bases_id[tetra[2]]+64*bases_id[tetra[3]]

def reverse(seq):
    return ''.join(complement[base] for base in reversed(seq))

def random_de_bruijn(n: int) -> str:
    alphabet = copy.deepcopy(bases)
    k = len(alphabet)
    random.shuffle(alphabet)  # Randomize the alphabet order
    
    a = [0] * k * n
    sequence = []

    def db(t, p):
        if t > n:
            if n % p == 0:
                sequence.extend(a[1 : p + 1])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)

    db(1, 1)
    return "".join(alphabet[i] for i in sequence)


#n5

#check how many bps we need
nbps_n5 = N_n5*2
#generate de bruijn sequence
seq_n5 = random_de_bruijn(4)
#append random bases to cover required length
nbps_DB = len(seq_n5)

print(len(seq_n5))
print(seq_n5)

for i in range(nbps_n5-nbps_DB+3):
    seq_n5 += bases[random.randint(0, 3)]

print(len(seq_n5))

#patch into sequences

SEQS_n5 = []

start = 0
delta = 5-3;
cutoff_T=10

while (start+5)<len(seq_n5):
    seq=""
    for i in range(5):
        seq+=seq_n5[start+i]
    print(start, seq)
    tmp = Sequence(seq)
    if tmp.mT>=cutoff_T:
        SEQS_n5.append(tmp)
    start+=delta



SEQS_n5.sort(key=lambda x: x.mT)

#n8+n15

padding_ls = []
for i in range(N_n8):
    padding_ls.append(8)
for i in range(N_n15):
    padding_ls.append(15)
    
random.shuffle(padding_ls)


#check how many bps we need
nbps_n8_n15 = N_n8*5+N_n15*12+3
#generate de bruijn sequence
seq_n8_n15 = random_de_bruijn(4)
#append random bases to cover required length
nbps_DB = len(seq_n8_n15)

print(len(seq_n8_n15))
print(seq_n8_n15)

for i in range(nbps_n8_n15-nbps_DB):
    seq_n8_n15 += bases[random.randint(0, 3)]

print(len(seq_n8_n15))

#patch into sequences

SEQS_n8 = []
SEQS_n15 = []

start = 0
cutoff_T=10

start = 0

for n in range(len(padding_ls)):
    seq=""
    for i in range(padding_ls[n]):
        seq+=seq_n8_n15[start+i]
    print(start, seq)
    tmp = Sequence(seq)
    if tmp.mT>=cutoff_T:
        if padding_ls[n] == 8:
            SEQS_n8.append(tmp)
        if padding_ls[n] == 15:
            SEQS_n15.append(tmp)
    start+=padding_ls[n]-3



SEQS_n8.sort(key=lambda x: x.mT)
SEQS_n15.sort(key=lambda x: x.mT)

print("Number of n5 sequences: ",len(SEQS_n5))
for i in range(len(SEQS_n5)):
    print(SEQS_n5[i].seq, SEQS_n5[i].mT)

print("Number of n8 sequences: ",len(SEQS_n8))
for i in range(len(SEQS_n8)):
    print(SEQS_n8[i].seq, SEQS_n8[i].mT)
    
print("Number of n15 sequences: ",len(SEQS_n15))
for i in range(len(SEQS_n15)):
    print(SEQS_n15[i].seq, SEQS_n15[i].mT)








