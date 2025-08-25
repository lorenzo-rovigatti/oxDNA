import numpy as np
import sys
from pathlib import Path


if len(sys.argv) != 5 :
    print("Invalid syntax")
    print("Usage: "+sys.argv[0]+" T q dq nbs")
    print("T must be in C!")
    exit(1)
    

backwards = False

T = sys.argv[1] #in C
Tox = (T+273.15)/3000. #oxdna units
q = float(sys.argv[2])
dq = float(sys.argv[3])
nbs = int(sys.argv[4])

if dq < 0: backwards = True

En_bn = np.loadtxt("split_energy_bn.dat")
En_un = np.loadtxt("split_energy_un.dat")

dh_ave_bn = np.mean(En_bn[:,8])*nbs
pot_ave_bn = np.mean(En_bn[:,9])*nbs

dh_ave_un = np.mean(En_un[:,8])*nbs
pot_ave_un = np.mean(En_un[:,9])*nbs

ofile = open(integrated_averages.dat, 'a')
print(q, T, dh_ave_bn, pot_ave_bn, dh_ave_un, pot_ave_un, file=ofile)
ofile.close()


dTox = Tox*dq*(2/q)*(dh_ave_un-dh_ave_bn)/(pot_ave_un-pot_ave_bn) #integrate

dT = dTox*3000-273.15 #in C

if backwards == False : 
    ofile = open(integrated_mT_vs_q.dat, 'a') #append
    print(q+dq, T+dT, file=ofile)
    ofile.close()
else : 
    ofile = Path(integrated_mT_vs_q.dat)
    line = str(q+dq) + " " + str(T+dT)
    ofile.write_text(line+ofile.read_text()) #prepend

