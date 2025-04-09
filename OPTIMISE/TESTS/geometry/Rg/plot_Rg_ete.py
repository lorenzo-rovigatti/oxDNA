import numpy as np
import math
import matplotlib.pyplot as plt
import "../../colours.py"

Ns = [20, 40, 60, 80, 100]

#style

plt.style.use("../../style.sty")

long_range = 100

def get_long_range(Rgs) :
    Rg = 0.
    Rg2 = 0.
    ete = 0.
    ete2 = 0.

    if long_range >= len(Rgs): 
        print("Long range too long.")
        exit(1)
    entries = len(Rgs)-long_range
    for i in range(long_range, len(Rgs)) :
        Rg += Rgs[i,1]/entries
        Rg2 += Rgs[i,1]*Rgs[i,1]/entries
        ete += Rgs[i,2]/entries
        ete2 += Rgs[i,2]*Rg[i,2]/entries
    return Rg, math.sqrt(Rg2-Rg*Rg), ete, math.sqrt(ete2-ete*ete) 

#data

Rgs_20 = np.loadtxt("0/Rg_ete.txt")
Rgs_40 = np.loadtxt("1/Rg_ete.txt")
Rgs_60 = np.loadtxt("2/Rg_ete.txt")
Rgs_80 = np.loadtxt("3/Rg_ete.txt")
Rgs_100 = np.loadtxt("4/Rg_ete.txt")


#compute long range Rgs and ete

lr_Rgs = []
lr_Rgs_err = []
lr_ete = []
lr_ete_err = []

rg, rge, ete, etee = get_long_range(Rgs_20)
lr_Rgs.append(rg)
lr_Rgs_err.append(rge)
lr_ete.append(ete)
lr_ete_err.append(etee)

rg, rge, ete, etee = get_long_range(Rgs_40)
lr_Rgs.append(rg)
lr_Rgs_err.append(rge)
lr_ete.append(ete)
lr_ete_err.append(etee)

rg, rge, ete, etee = get_long_range(Rgs_60)
lr_Rgs.append(rg)
lr_Rgs_err.append(rge)
lr_ete.append(ete)
lr_ete_err.append(etee)

rg, rge, ete, etee = get_long_range(Rgs_80)
lr_Rgs.append(rg)
lr_Rgs_err.append(rge)
lr_ete.append(ete)
lr_ete_err.append(etee)

rg, rge, ete, etee = get_long_range(Rgs_100)
lr_Rgs.append(rg)
lr_Rgs_err.append(rge)
lr_ete.append(ete)
lr_ete_err.append(etee)


#versus time

#Rg

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$m$",fontsize=20)
ax.set_ylabel(r"radius of gyration, $R_g$ [nm]",fontsize=20)
#ax.set_ylim(30,60)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(Rgs_20[:,0],Rgs_20[:,1],color=CCS_G[1],label="20mers")
ax.plot(Rgs_40[:,0],Rgs_40[:,1],color=CCS_G[12],label="40mers")
ax.plot(Rgs_60[:,0],Rgs_60[:,1],color=CCS_G[10],label="60mers")
ax.plot(Rgs_80[:,0],Rgs_80[:,1],color=CCS_G[7],label="80mers")
ax.plot(Rgs_100[:,0],Rgs_100[:,1],color=CCS_G[11],label="100mers")
ax.legend(fontsize = 20)
plt.savefig("Rg_v_time.pdf",bbox_inches='tight',pad_inches=0.05)

#end to end

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$m$",fontsize=20)
ax.set_ylabel(r"end to end distance [nm]",fontsize=20)
#ax.set_ylim(80,200)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(Rgs_20[:,0],Rgs_20[:,2],color=CCS_G[1],label="20mers")
ax.plot(Rgs_40[:,0],Rgs_40[:,2],color=CCS_G[12],label="40mers")
ax.plot(Rgs_60[:,0],Rgs_60[:,2],color=CCS_G[10],label="60mers")
ax.plot(Rgs_80[:,0],Rgs_80[:,2],color=CCS_G[7],label="80mers")
ax.plot(Rgs_100[:,0],Rgs_100[:,2],color=CCS_G[11],label="100mers")
ax.legend(fontsize = 20)
plt.savefig("Ete_vs_time.pdf",bbox_inches='tight',pad_inches=0.05)

#long times

#Rg

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$m$",fontsize=20)
ax.set_ylabel(r"radius of gyration, $R_g$ [nm]",fontsize=20)
#ax.set_ylim(30,60)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.errorbars(Ns,lr_Rgs,yerr=lr_Rgs_err,color=CCS_G[1],label=None)
plt.savefig("Rg_long_times.pdf",bbox_inches='tight',pad_inches=0.05)

#ete

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$m$",fontsize=20)
ax.set_ylabel(r"end to end distance [nm]",fontsize=20)
#ax.set_ylim(30,60)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.errorbars(Ns,lr_ete,yerr=lr_ete_err,color=CCS_G[1],label=None)
plt.savefig("Ete_long_times.pdf",bbox_inches='tight',pad_inches=0.05)

