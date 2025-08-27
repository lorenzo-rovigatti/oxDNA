import numpy as np
import math
import matplotlib.pyplot as plt
import "../../colours.py"

poly_steps = ['AA', 'AG', 'AC', 'AT', 'CG', 'CC']
richness = [10, 30, 50, 70, 90]

#style

plt.style.use("../../style.sty")


#data

lps_AA = np.loadtxt("0/lps_vs_m.txt")
lps_AG = np.loadtxt("1/lps_vs_m.txt")
lps_AC = np.loadtxt("2/lps_vs_m.txt")
lps_AT = np.loadtxt("3/lps_vs_m.txt")
lps_CG = np.loadtxt("4/lps_vs_m.txt")
lps_CC = np.loadtxt("5/lps_vs_m.txt")

lps_10 = np.loadtxt("6/lps_vs_m.txt")
lps_30 = np.loadtxt("7/lps_vs_m.txt")
lps_50 = np.loadtxt("8/lps_vs_m.txt")
lps_70 = np.loadtxt("9/lps_vs_m.txt")
lps_90 = np.loadtxt("10/lps_vs_m.txt")

#repeats

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
ax.plot(lps_AA[:,0],lps_AA[:,1],color=CCS_G[1],label="AA")
ax.plot(lps_AG[:,0],lps_AG[:,1],color=CCS_G[12],label="AG")
ax.plot(lps_AC[:,0],lps_AC[:,1],color=CCS_G[10],label="AC")
ax.plot(lps_AT[:,0],lps_AT[:,1],color=CCS_G[7],label="AT")
ax.plot(lps_CC[:,0],lps_CC[:,1],color=CCS_G[11],label="CC")
ax.plot(lps_CG[:,0],lps_CG[:,1],color=CCS_G[8],label="CG")
ax.legend(fontsize = 20)
plt.savefig("lbs_repeats.pdf",bbox_inches='tight',pad_inches=0.05)

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
ax.plot(lps_AA[:,0],lps_AA[:,2]/2,color=CCS_G[1],label="AA")
ax.plot(lps_AG[:,0],lps_AG[:,2]/2,color=CCS_G[12],label="AG")
ax.plot(lps_AC[:,0],lps_AC[:,2]/2,color=CCS_G[10],label="AC")
ax.plot(lps_AT[:,0],lps_AT[:,2]/2,color=CCS_G[7],label="AT")
ax.plot(lps_CC[:,0],lps_CC[:,2]/2,color=CCS_G[11],label="CC")
ax.plot(lps_CG[:,0],lps_CG[:,2]/2,color=CCS_G[8],label="CG")
ax.legend(fontsize = 20)
plt.savefig("lts_repeats.pdf",bbox_inches='tight',pad_inches=0.05)


#richness
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
ax.plot(lps_10[:,0],lps_10[:,1],color=CCS_G[1],label="10% AT")
ax.plot(lps_30[:,0],lps_30[:,1],color=CCS_G[12],label="30% AT")
ax.plot(lps_50[:,0],lps_50[:,1],color=CCS_G[10],label="50% AT")
ax.plot(lps_70[:,0],lps_70[:,1],color=CCS_G[7],label="70% AT")
ax.plot(lps_90[:,0],lps_90[:,1],color=CCS_G[11],label="90% AT")
ax.legend(fontsize = 20)
plt.savefig("lbs_AT_rich.pdf",bbox_inches='tight',pad_inches=0.05)

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
ax.plot(lps_10[:,0],lps_10[:,2]/2,color=CCS_G[1],label="10% AT")
ax.plot(lps_30[:,0],lps_30[:,2]/2,color=CCS_G[12],label="30% AT")
ax.plot(lps_50[:,0],lps_50[:,2]/2,color=CCS_G[10],label="50% AT")
ax.plot(lps_70[:,0],lps_70[:,2]/2,color=CCS_G[7],label="70% AT")
ax.plot(lps_90[:,0],lps_90[:,2]/2,color=CCS_G[11],label="90% AT")
ax.legend(fontsize = 20)
plt.savefig("lts_AT_rich.pdf",bbox_inches='tight',pad_inches=0.05)


