import numpy as np
import matplotlib.pyplot as plt

grooves = None
pitch_lengths = None


for i in range(5):
    for j in range(5):
        if i == 0 and j == 0: continue
        #get grooves
        if i ==0 and j == 1: 
        	grooves = np.loadtxt(str(i)+str(j)+"/grooves.txt")
        	pitch_lengths = np.loadtxt(str(i)+str(j)+"/pitch.dat",dtype=str)[:,2:4].astype(np.float64)
        else :
        	grooves=np.vstack((grooves,np.loadtxt(str(i)+str(j)+"/grooves.txt")))
        	pitch_lengths = np.vstack((pitch_lengths,np.loadtxt(str(i)+str(j)+"/pitch.dat",dtype=str)[:,2:4].astype(np.float64)))
        	
        	
grooves = grooves - 5.8


print(pitch_lengths)


#colours

CCS_G = []
CCS_G.append('#000000')
CCS_G.append('#006884')
CCS_G.append('#00909E')
CCS_G.append('#5B5B5B')
CCS_G.append('#6E006C')
CCS_G.append('#89DBEC')
CCS_G.append('#91278F')
CCS_G.append('#B00051')
CCS_G.append('#CF97D7')
CCS_G.append('#D4D4D4')
CCS_G.append('#ED0026')
CCS_G.append('#F68370')
CCS_G.append('#FA9D00')
CCS_G.append('#FEABB9')
CCS_G.append('#FFD08D')

#style

plt.style.use("/work/Desktop/Topo_Roux/Analysis/style.sty")


fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"groove width [A]",fontsize=20)
ax.set_ylabel(r"distribution",fontsize=20)
#ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

plt.hist(grooves[:,0], bins=20, color=CCS_G[1],label="minor", density=True)
plt.hist(grooves[:,1], bins=20, color=CCS_G[7],label="major", density=True)

ax.legend(fontsize = 20)
plt.savefig("grooves_hist.pdf",bbox_inches='tight',pad_inches=0.05)


fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"pitch length [A]",fontsize=20)
ax.set_ylabel(r"distribution",fontsize=20)
#ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

plt.hist(pitch_lengths[:,1], bins=20, color=CCS_G[1], density=True)

ax.legend(fontsize = 20)
plt.savefig("pitch_hist.pdf",bbox_inches='tight',pad_inches=0.05)

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"pitch length [bps]",fontsize=20)
ax.set_ylabel(r"counts",fontsize=20)
#ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

plt.hist(pitch_lengths[:,0], bins=20, color=CCS_G[1], density=True)

ax.legend(fontsize = 20)
plt.savefig("pitch_hist_bps.pdf",bbox_inches='tight',pad_inches=0.05)


