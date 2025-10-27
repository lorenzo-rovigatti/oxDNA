import numpy
import sys
import re

box = 20
Sc = 0.5
deltaT = 7


print("Box size = ",box)
print("Salt concentration = ",Sc)

seqs=[]
op_files=[]
w_files=[]
mTs=[]

if len(sys.argv) != 2:
    print("Usage: "+sys.argv[0]+ " input_file")
    exit(1)

def prepare_script() :

   ifile=open(sys.argv[1], 'r')
   for line in ifile.readlines():
       entries=re.split(r'\s+',line.strip())
       if entries[0] == "seqs=(":
           for n in range(1,len(entries)):
               if entries[n] != ")":
                   seqs.append(entries[n])
       if entries[0] == "OP_files=(":
           for n in range(1,len(entries)):
               if entries[n] != ")":
                   op_files.append(entries[n])
       if entries[0] == "W_files=(":
           for n in range(1,len(entries)):
               if entries[n] != ")":
                   w_files.append(entries[n])
       if entries[0] == "mTs=(":
           for n in range(1,len(entries)):
               if entries[n] != ")":
                   mTs.append(int(entries[n]))

   ifile.close()
   ofile=open("melting_script_SEQs.txt", 'w')
   for n in range(len(seqs)) :
       line = "SEQ " + seqs[n] + " " + str(box) + " " + str(Sc) + " 5 "
       T0 = mTs[n]-deltaT*2
       for j in range(5):
           T = T0+j*deltaT
           line += str(T) + " "
       line += "2.0 " + op_files[n] + " "
       for j in range(5):
           line += w_files[n]+str(j)+" "

       print(line,file=ofile)

   ofile.close()

   return


prepare_script()
