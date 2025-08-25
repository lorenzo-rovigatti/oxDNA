import numpy as np
import matplotlib.pyplot as plt
plt.style.use("../../style.sty")
import "../../colours.py"

grooves = None
pitch_lengths = None

Cs=[0.1, 0.15, 0.2, 0.3, 0.5, 1.1]

Nseqs=len(Cs)

grooves_ave = []
grooves_ave_err = []
pitch_ave = []
pitch_ave_err = []

for i in range(Nseqs):
    grooves = np.loadtxt(str(i)+"/grooves.txt")-5.8
    pitch_lengths = np.loadtxt(str(i)+"/pitch.dat",dtype=str)[:,2:-1].astype(np.float64)

    grooves_ave.append(np.mean(grooves, axis=0))
    grooves_ave_err.append(np.std(grooves, axis=0))
    pitch_ave.append(np.mean(pitch_lengths, axis=0))
    pitch_ave_err.append(np.std(pitch_lengths, axis=0))


#grooves

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_ylabel(r"groove width [A]",fontsize=20)
ax.set_xlabel(r"salt concentration, $\left[ Na^{+} \right]/M$",fontsize=20)
#ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

plt.errorbars(Cs, grooves_ave[:,0], yerr=grooves_ave_err[:,0], color=CCS_G[1],label="minor", density=True)
plt.errorbars(Cs, grooves_ave[:,1], yerr=grooves_ave_err[:,1], color=CCS_G[7],label="major", density=True)

ax.legend(fontsize = 20)
plt.savefig("grooves_vs_cs.pdf",bbox_inches='tight',pad_inches=0.05)


#pitch

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_ylabel(r"pitch length [A]",fontsize=20)
ax.set_xlabel(r"salt concentration, $\left[ Na^{+} \right]/M$",fontsize=20)
#ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

plt.errorbars(Cs, pitch_ave[:,0], yerr=pitch_ave_err[:,0], color=CCS_G[1],label=None, density=True)

ax.legend(fontsize = 20)
plt.savefig("pitch_vs_cs.pdf",bbox_inches='tight',pad_inches=0.05)

#rise

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_ylabel(r"rise [A]",fontsize=20)
ax.set_xlabel(r"salt concentration, $\left[ Na^{+} \right]/M$",fontsize=20)
#ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

plt.errorbars(Cs, pitch_ave[:,1], yerr=pitch_ave_err[:,1], color=CCS_G[1], label=None, density=True)

ax.legend(fontsize = 20)
plt.savefig("rise_ave_vs_cs.pdf",bbox_inches='tight',pad_inches=0.05)

#twist

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_ylabel(r"twist [deg]",fontsize=20)
ax.set_xlabel(r"salt concentration, $\left[ Na^{+} \right]/M$",fontsize=20)
#ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

plt.errorbars(Cs, pitch_ave[:,2], yerr=pitch_ave_err[:,2], color=CCS_G[1], label=None, density=True)

ax.legend(fontsize = 20)
plt.savefig("twist_ave_vs_cs.pdf",bbox_inches='tight',pad_inches=0.05)

