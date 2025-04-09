import numpy as np
import math
import matplotlib.pyplot as plt
import "../../colours.py"

Cs=[0.1, 0.15, 0.2, 0.3, 0.5, 1.1]

#style

plt.style.use("../../style.sty")


#data

lps_01 = np.loadtxt("0/lps_vs_m.txt")
lps_015 = np.loadtxt("1/lps_vs_m.txt")
lps_02 = np.loadtxt("2/lps_vs_m.txt")
lps_03 = np.loadtxt("3/lps_vs_m.txt")
lps_05 = np.loadtxt("4/lps_vs_m.txt")
lps_11 = np.loadtxt("5/lps_vs_m.txt")

#salt_concentration
#bending

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$m$",fontsize=20)
ax.set_ylabel(r"bending persistence length, $l_b$ [nm]",fontsize=20)
ax.set_ylim(30,60)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(lps_01[:,0],lps_01[,1],color=CCS_G[1],label="0.1")
ax.plot(lps_015[:,0],lps_015[:,1],color=CCS_G[12],label="0.15")
ax.plot(lps_02[:,0],lps_02[:,1],color=CCS_G[10],label="0.2")
ax.plot(lps_03[:,0],lps_03[:,1],color=CCS_G[7],label="0.3")
ax.plot(lps_05[:,0],lps_05[:,1],color=CCS_G[11],label="0.5")
ax.plot(lps_11[:,0],lps_11[:,1],color=CCS_G[11],label="1.1")
ax.legend(fontsize = 20)
plt.savefig("lbs_Cs.pdf",bbox_inches='tight',pad_inches=0.05)

#tortional

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$m$",fontsize=20)
ax.set_ylabel(r"half tortional per. length, $l_t/2$ [nm]",fontsize=20)
ax.set_ylim(80,200)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(lps_01[:,0],lps_01[,2]/2,color=CCS_G[1],label="0.1")
ax.plot(lps_015[:,0],lps_015[:,2]/2,color=CCS_G[12],label="0.15")
ax.plot(lps_02[:,0],lps_02[:,2]/2,color=CCS_G[10],label="0.2")
ax.plot(lps_03[:,0],lps_03[:,2]/2,color=CCS_G[7],label="0.3")
ax.plot(lps_05[:,0],lps_05[:,2]/2,color=CCS_G[11],label="0.5")
ax.plot(lps_11[:,0],lps_11[:,2]/2,color=CCS_G[11],label="1.1")
ax.legend(fontsize = 20)
plt.savefig("lts_Cs.pdf",bbox_inches='tight',pad_inches=0.05)


#limit values

limits = np.loadtxt("limit_lps.txt")

#bending

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"salt concentration, $\left[ Na^{+} \right]/M",fontsize=20)
ax.set_ylabel(r"bending persistence length, $l_b$ [nm]",fontsize=20)
ax.set_ylim(35,55)
ax.set_xlim(0.05,1.2)
plt.xscale('log')
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.errorbar(limits[:,0],limits[,1],yerr=limits[,2],color=CCS_G[1],label=None)
ax.legend(fontsize = 20)
plt.savefig("limit_lbs_Cs.pdf",bbox_inches='tight',pad_inches=0.05)

#tortional

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"salt concentration, $\left[ Na^{+} \right]/M",fontsize=20)
ax.set_ylabel(r"half tortional per. length, $l_t/2$ [nm]",fontsize=20)
ax.set_ylim(110,140)
ax.set_xlim(0.05,1.2)
plt.xscale('log')
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.errorbar(limits[:,0],limits[,2],yerr=limits[,3],color=CCS_G[1],label=None)
ax.legend(fontsize = 20)
plt.savefig("limit_lts_Cs.pdf",bbox_inches='tight',pad_inches=0.05)

