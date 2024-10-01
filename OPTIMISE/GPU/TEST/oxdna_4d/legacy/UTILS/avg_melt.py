#!/usr/bin/env python


import sys
import math
saltcon = 0.50
if len(sys.argv) < 3:
	print 'Usage: ./program LENGTH box [salt, default 0.5]'
	sys.exit(1)

box = float(sys.argv[2])
length = int(sys.argv[1]) - 1
if len(sys.argv) >= 4:
  saltcon = float(sys.argv[3])



molcon = (2/(8.5179*10.0**(-9)*box)**3)/(6.0221415*10**23)
# See end of while for explanation of why we use length - 1
saltcorr =  0.368*(length -1 ) * math.log(saltcon)

DH = (-8.2375*length + 1.2*2)
DS = (-22.0188*length + 1.2) + saltcorr

RealTemp = (DH )*1000 / (DS + 1.9859* math.log(molcon/4) )

print 'Salt concentration = ', str(saltcon)
print 'T = ',str(RealTemp)+'K'
print RealTemp - 273.15

# ----- why do we use length - 1?
# In the SantaLucia model (described in the paper "THE THERMODYNAMICS OF
# DNA STRUCTURAL MOTIFS"), the salt concentration modifies the entropic term
# in eq. 5. of the paper. The term is proportional to N/2, where N is the
# number of phosphates. If we had one phosphate per nucleotide we will
# therefore have to use 0.368*length* math.log(saltcon), but since terminal
# nucleotides only have half a phosphate in oxDNA2 in total we lack 4*0.5 = 2
# phosphates. THerefore we write 0.368*(length - 1)*math.log(saltcon)
