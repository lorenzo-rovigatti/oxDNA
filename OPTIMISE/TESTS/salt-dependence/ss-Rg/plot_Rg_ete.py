import numpy as np
import math
import matplotlib.pyplot as plt
import "../../colours.py"

Ns = [20, 40, 60, 80, 100]
Cs=[0.1, 0.15, 0.2, 0.3, 0.5, 1.1]

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

Rgs = [[[] for _ in range(Cs)] for _ in range(Ns)]  #[salt][N]

lr_Rgs = [[[] for _ in range(Cs)] for _ in range(Ns)]  #[salt][N] = [rg, rge, ete, etee], long range

for i in range(len(Cs)):
    for j in range(len(Ns)):
        Rgs[i][j] = np.loadtxt(str(i)+"/"+str(j)+"/Rg_ete.txt")

        #compute long range Rgs and ete
        rg, rge, ete, etee = get_long_range(Rgs[i][j])
        lr_Rgs[i][j].append(rg)
        lr_Rgs[i][j].append(rge)
        lr_Rgs[i][j].append(ete)
        lr_Rgs[i][j].append(etee)



#versus time

#Rg

#N=20

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
ax.plot(Rgs[0][0][:,0],Rgs[0][0][:,1],color=CCS_G[1],label="0.1 M")
ax.plot(Rgs[1][0][:,0],Rgs[1][0][:,1],color=CCS_G[12],label="0.15 M")
ax.plot(Rgs[2][0][:,0],Rgs[2][0][:,1],color=CCS_G[10],label="0.2 M")
ax.plot(Rgs[3][0][:,0],Rgs[3][0][:,1],color=CCS_G[7],label="0.3 M")
ax.plot(Rgs[4][0][:,0],Rgs[4][0][:,1],color=CCS_G[11],label="0.5 M")
ax.plot(Rgs[5][0][:,0],Rgs[5][0][:,1],color=CCS_G[8],label="1.1 M")
ax.legend(fontsize = 20)
plt.savefig("Rg_v_time_N20.pdf",bbox_inches='tight',pad_inches=0.05)

#N=40

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
ax.plot(Rgs[0][1][:,0],Rgs[0][1][:,1],color=CCS_G[1],label="0.1 M")
ax.plot(Rgs[1][1][:,0],Rgs[1][1][:,1],color=CCS_G[12],label="0.15 M")
ax.plot(Rgs[2][1][:,0],Rgs[2][1][:,1],color=CCS_G[10],label="0.2 M")
ax.plot(Rgs[3][1][:,0],Rgs[3][1][:,1],color=CCS_G[7],label="0.3 M")
ax.plot(Rgs[4][1][:,0],Rgs[4][1][:,1],color=CCS_G[11],label="0.5 M")
ax.plot(Rgs[5][1][:,0],Rgs[5][1][:,1],color=CCS_G[8],label="1.1 M")
ax.legend(fontsize = 20)
plt.savefig("Rg_v_time_N40.pdf",bbox_inches='tight',pad_inches=0.05)

#N=60

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
ax.plot(Rgs[0][2][:,0],Rgs[0][2][:,1],color=CCS_G[1],label="0.1 M")
ax.plot(Rgs[1][2][:,0],Rgs[1][2][:,1],color=CCS_G[12],label="0.15 M")
ax.plot(Rgs[2][2][:,0],Rgs[2][2][:,1],color=CCS_G[10],label="0.2 M")
ax.plot(Rgs[3][2][:,0],Rgs[3][2][:,1],color=CCS_G[7],label="0.3 M")
ax.plot(Rgs[4][2][:,0],Rgs[4][2][:,1],color=CCS_G[11],label="0.5 M")
ax.plot(Rgs[5][2][:,0],Rgs[5][2][:,1],color=CCS_G[8],label="1.1 M")
ax.legend(fontsize = 20)
plt.savefig("Rg_v_time_N60.pdf",bbox_inches='tight',pad_inches=0.05)

#end to end

#N=20

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
ax.plot(Rgs[0][0][:,0],Rgs[0][0][:,2],color=CCS_G[1],label="0.1 M")
ax.plot(Rgs[1][0][:,0],Rgs[1][0][:,2],color=CCS_G[12],label="0.15 M")
ax.plot(Rgs[2][0][:,0],Rgs[2][0][:,2],color=CCS_G[10],label="0.2 M")
ax.plot(Rgs[3][0][:,0],Rgs[3][0][:,2],color=CCS_G[7],label="0.3 M")
ax.plot(Rgs[4][0][:,0],Rgs[4][0][:,2],color=CCS_G[11],label="0.5 M")
ax.plot(Rgs[5][0][:,0],Rgs[5][0][:,2],color=CCS_G[8],label="1.1 M")
ax.legend(fontsize = 20)
plt.savefig("ete_v_time_cs_N20.pdf",bbox_inches='tight',pad_inches=0.05)

#N=40

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
ax.plot(Rgs[0][1][:,0],Rgs[0][1][:,2],color=CCS_G[1],label="0.1 M")
ax.plot(Rgs[1][1][:,0],Rgs[1][1][:,2],color=CCS_G[12],label="0.15 M")
ax.plot(Rgs[2][1][:,0],Rgs[2][1][:,2],color=CCS_G[10],label="0.2 M")
ax.plot(Rgs[3][1][:,0],Rgs[3][1][:,2],color=CCS_G[7],label="0.3 M")
ax.plot(Rgs[4][1][:,0],Rgs[4][1][:,2],color=CCS_G[11],label="0.5 M")
ax.plot(Rgs[5][1][:,0],Rgs[5][1][:,2],color=CCS_G[8],label="1.1 M")
ax.legend(fontsize = 20)
plt.savefig("ete_v_time_cs_N40.pdf",bbox_inches='tight',pad_inches=0.05)

#N=60

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
ax.plot(Rgs[0][2][:,0],Rgs[0][2][:,2],color=CCS_G[1],label="0.1 M")
ax.plot(Rgs[1][2][:,0],Rgs[1][2][:,2],color=CCS_G[12],label="0.15 M")
ax.plot(Rgs[2][2][:,0],Rgs[2][2][:,2],color=CCS_G[10],label="0.2 M")
ax.plot(Rgs[3][2][:,0],Rgs[3][2][:,2],color=CCS_G[7],label="0.3 M")
ax.plot(Rgs[4][2][:,0],Rgs[4][2][:,2],color=CCS_G[11],label="0.5 M")
ax.plot(Rgs[5][2][:,0],Rgs[5][2][:,2],color=CCS_G[8],label="1.1 M")
ax.legend(fontsize = 20)
plt.savefig("ete_v_time_cs_N60.pdf",bbox_inches='tight',pad_inches=0.05)


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
ax.errorbars(Ns,lr_Rgs[0][:][0],yerr=lr_Rgs[0][:][1],color=CCS_G[1],label="0.1 M")
ax.errorbars(Ns,lr_Rgs[1][:][0],yerr=lr_Rgs[1][:][1],color=CCS_G[12],label="0.15 M")
ax.errorbars(Ns,lr_Rgs[2][:][0],yerr=lr_Rgs[2][:][1],color=CCS_G[1]0,label="0.2 M")
ax.errorbars(Ns,lr_Rgs[3][:][0],yerr=lr_Rgs[3][:][1],color=CCS_G[7],label="0.3 M")
ax.errorbars(Ns,lr_Rgs[4][:][0],yerr=lr_Rgs[4][:][1],color=CCS_G[11],label="0.5 M")
ax.errorbars(Ns,lr_Rgs[5][:][0],yerr=lr_Rgs[5][:][1],color=CCS_G[8],label="1.1 M")
ax.legend(fontsize = 20)
plt.savefig("Rg_long_times_cs.pdf",bbox_inches='tight',pad_inches=0.05)

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
ax.errorbars(Ns,lr_Rgs[0][:][2],yerr=lr_Rgs[0][:][3],color=CCS_G[1],label="0.1 M")
ax.errorbars(Ns,lr_Rgs[1][:][2],yerr=lr_Rgs[1][:][3],color=CCS_G[12],label="0.15 M")
ax.errorbars(Ns,lr_Rgs[2][:][2],yerr=lr_Rgs[2][:][3],color=CCS_G[1]0,label="0.2 M")
ax.errorbars(Ns,lr_Rgs[3][:][2],yerr=lr_Rgs[3][:][3],color=CCS_G[7],label="0.3 M")
ax.errorbars(Ns,lr_Rgs[4][:][2],yerr=lr_Rgs[4][:][3],color=CCS_G[11],label="0.5 M")
ax.errorbars(Ns,lr_Rgs[5][:][2],yerr=lr_Rgs[5][:][3],color=CCS_G[8],label="1.1 M")
ax.legend(fontsize = 20)
plt.savefig("ete_long_times_cs.pdf",bbox_inches='tight',pad_inches=0.05)

