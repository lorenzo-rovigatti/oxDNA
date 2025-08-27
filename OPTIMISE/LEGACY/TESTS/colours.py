import numpy as np
import math
import matplotlib.pyplot as plt

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

if __name__ == "__main__" :
    Ncols=len(CCS_G)

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    plt.title(r"CCS_G palette")
    ax.title.set_fontsize(20)
    ax.set_ylabel(r"colour id",fontsize=20)
    ax.set_ylim(-1,Ncols)
    ax.set_xlim(0,10)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=8)

    x = []
    for j in range(11): x.append(j)

    for i in range(Ncols) :

        y = []
        for j in range(11): y.append(i)
    
        ax.plot(x,y,color=CCS_G[i], label=None)
    
    plt.savefig("palette.pdf",bbox_inches='tight',pad_inches=0.05)






