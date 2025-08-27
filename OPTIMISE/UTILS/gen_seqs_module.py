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

max_gen_N=200000
Ct = 0.000672
Cs = 0.5

class Sequence :

    def __init__(self,seq_int) :

        SEQ = ""

        for i in range(len(seq_int)) :
            SEQ += SL.id_to_base(seq_int[i])

        MT = SL.melting_temperature(SEQ,Ct,Cs)

        self.seq = SEQ
        self.mT = MT


def set_concentartions(ct,cs) :
    global Ct
    global Cs
    Ct = ct
    Cs = cs
    return


SEQS = []

def gen_seqs(bps) :
    for n in range(max_gen_N):
        seq = np.zeros(bps)
        if n%100000 == 0:
            print("Generated",n,"sequences")
        for i in range(bps):
            seq[i]=np.random.randint(0,4)
        SEQS.append(Sequence(seq))

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

def search_complementary(sampled,s) :
    s_rev= ""
    for i in range(len(s.seq)-1,-1,-1):
        s_rev+=comple[s.seq[i]]

    for i in sampled :
        if i.seq == s_rev :
            return True

    return False


def generate_sequences(bps) :

    seq = np.zeros(bps,dtype=int)

    tot_nseqs = pow(4,bps)

    print("Total number of sequences with bps: ",bps," = ",tot_nseqs)

    if tot_nseqs > max_gen_N :
        gen_seqs(bps)
    else:
        #generate all sequences
        loop_seqs_cutoff(seq,bps)

        print("using cutoff 12 C.")

    #sort sequences in ascending melting temperature
    SEQS.sort(key=lambda x: x.mT)

    print(len(SEQS))

    print("max mT: "+ SEQS[len(SEQS)-1].seq + " " +str(SEQS[len(SEQS)-1].mT))
    print("min mT: "+ SEQS[0].seq + " "+ str(SEQS[0].mT))

    return

def sample(seqs,N,coarse,uni_T=True) :

    delta = int(len(seqs)/N)

    sampled = []

    if uni_T :

        #uniform mT sampling

        delta = coarse*(seqs[len(seqs)-1].mT - seqs[0].mT)/N

        k0 = 0

        for n in range(int(N/coarse)) :

            tmT = SEQS[0].mT+(n+1)*delta

            for k in range(k0,len(seqs)-1) :
                if seqs[k].mT <= tmT and (seqs[k+1].mT >= tmT or k == len(seqs)-2):

                    sampled_ids = []
                    while len(sampled_ids) < coarse:
                        r = random.randint(k0, k)
                        if r in sampled_ids or search_complementary(sampled,seqs[r]):
                            continue
                        else:
                            sampled_ids.append(r)
                            sampled.append(seqs[r])

                    k0 = k
                    break

    else:

        #sampling accounting for population

        for n in range(N) :

            while 1 == 1:
                r = n*delta + random.randint(0, delta)

                if SEQS[r] in sampled or search_complementary(seqs[r],sampled) :
                    continue
                else :
                    sampled.append(seqs[r])
                    break

    sampled.sort(key=lambda x: x.mT)

    return sampled


def print_sequences(sampled,bps) :

    ofile = open("gen_seqs_n"+str(bps)+".txt", 'w')

    for seq in sampled :
        print(seq.seq + " " + str(seq.mT), file = ofile)

    ofile.close()

    ofile = open("for_melting_script_n"+str(bps)+".txt", 'w')
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
    return
