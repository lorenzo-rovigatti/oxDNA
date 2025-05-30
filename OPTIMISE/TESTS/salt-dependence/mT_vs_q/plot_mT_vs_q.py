import os
import numpy as np
import math
import matplotlib.pyplot as plt
import "../../colours.py"
import SantaLucia as SL

Cs=[0.1, 0.2, 0.5]

Ct = 0.0006717 #Important: we are assuming the box size is 20 oxdna units.

#style

plt.style.use("../../style.sty")

seqs_n5 = []
seqs_n8 = []
seqs_n15 = []

mTs_n5 = []
mTs_n8 = []
mTs_n15 = []

qs = []

ave_DmT_n5 = None
ave_DmT_n8 = None
ave_DmT_n15 = None

first = True

#data

for i in len(Cs) :
    for filename in os.listdir("Results/"+str(Cs[i])):
	with open(filename) as f:
        headers = f.readline().strip().split()
        seq = headers[3]
        mT = SL.melting_temperature(seq,Ct,Cs[i])
        if len(seq)==5:
            seqs_n5.append(seq)
            mTs_n5.append(mT)
            data = np.readtxt(filename)
            if first:
                qs = data[:,0]
                ave_DmT_n5 = np.zeros(len(qs))
                ave_DmT_n8 = np.zeros(len(qs))
                ave_DmT_n15 = np.zeros(len(qs))
                first = False
                
            ave_DmT_n5 += data[:,1]
        if len(seq)==8:
            seqs_n8.append(seq)
            mTs_n8.append(mT)
            data = np.readtxt(filename)
            if first:
                qs = data[:,0]
                ave_DmT_n5 = np.zeros(len(qs))
                ave_DmT_n8 = np.zeros(len(qs))
                ave_DmT_n15 = np.zeros(len(qs))
                first = False
                
            ave_DmT_n8 += data[:,1]
        if len(seq)==15:
            seqs_n15.append(seq)
            mTs_n15.append(mT)
            data = np.readtxt(filename)
            if first:
                qs = data[:,0]
                ave_DmT_n5 = np.zeros(len(qs))
                ave_DmT_n8 = np.zeros(len(qs))
                ave_DmT_n15 = np.zeros(len(qs))
                first = False
                
            ave_DmT_n15 += data[:,1]
            
        if len(seq) != 5 and len(seq) != 8 and len(seq) != 15:
            print("Warning: found a sequence of length ", len(seq))
            print("Only available length are: 5, 8, 15.")
            print("If you want to add another length, change the code.")
            


ave_DmT_n5 /= len(seqs_n5)
ave_DmT_n8 /= len(seqs_n8)
ave_DmT_n15 /= len(seqs_n15)


#find optimal qeff
tot_DmT = ave_DmT_n5+ave_DmT_n8+ave_DmT_n15

qopt = min(tot_DmT)

print("Optimal qeff: ", qopt)


#plot
fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$q_eff$",fontsize=20)
ax.set_ylabel(r"Av. melting temp. deviation $\Delta T_m$ [C]",fontsize=20)
#ax.set_ylim(30,60)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(qs,ave_DmT_n5,color=CCS_G[1],label="nbs = 5")
ax.plot(qs,ave_DmT_n8,color=CCS_G[7],label="nbs = 8")
ax.plot(qs,ave_DmT_n15,color=CCS_G[12],label="nbs = 15")
ax.axvline(x=qopt, color=CCS_G[0], linestyle='--', label="optimal qeff")
ax.legend(fontsize = 20)
plt.savefig("DmT_vs_q.pdf",bbox_inches='tight',pad_inches=0.05)

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$q_eff$",fontsize=20)
ax.set_ylabel(r"Av. melting temp. deviation $\Delta T_m$ [C]",fontsize=20)
#ax.set_ylim(30,60)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(qs,tot_DmT,color=CCS_G[1],label="nbs = 5,8,15")
ax.axvline(x=qopt, color=CCS_G[0], linestyle='--', label="optimal qeff")
ax.legend(fontsize = 20)
plt.savefig("tot_DmT_vs_q.pdf",bbox_inches='tight',pad_inches=0.05)

