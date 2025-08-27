import sys

if len(sys.argv) != 2:
    print("Invalid syntax.")
    print("Usage:",sys.argv[0],"file_with_sequences")
    exit(1)

class Sequence:
   def __init__(self,SEQ,MT):
       self.seq = SEQ
       self.mT = MT

seqs_n5 = []
seqs_n8 = []
seqs_n15 = []

ifile = open(sys.argv[1], 'r')

for line in ifile.readlines():
    vals = line.strip().split(" ")
    nbp = len(vals[0])
    if nbp == 5:
        seqs_n5.append( Sequence(vals[0],round(float(vals[1]))) )
    elif nbp == 8:
        seqs_n8.append( Sequence(vals[0],round(float(vals[1]))) )
    elif nbp == 15:
        seqs_n15.append( Sequence(vals[0],round(float(vals[1]))) )

ifile.close()

ofile = open("seqs_for_sampling_melting_script.txt",'w')

line_s = "seqs=( "
line_t = "mTs=( "
line_o = "OP_files=( "
line_w = "W_files=( "
line_b = "nbps=( "
for s in seqs_n5:
    line_s+=s.seq+" "
    line_t+=str(s.mT)+" "
    line_o+="op_n5.txt "
    line_w+="wfile_n5.txt "
    line_b+="5 "
line_s+=")"
line_t+=")"
line_o+=")"
line_w+=")"
line_b+=")"

print(line_s,file=ofile)
print(line_t,file=ofile)
print(line_o,file=ofile)
print(line_w,file=ofile)
print(line_b,file=ofile)

line_s = "seqs=( "
line_t = "mTs=( "
line_o = "OP_files=( "
line_w = "W_files=( "
line_b = "nbps=( "
for s in seqs_n8:
    line_s+=s.seq+" "
    line_t+=str(s.mT)+" "
    line_o+="op_n8.txt "
    line_w+="wfile_n8.txt "
    line_b+="8 "
line_s+=")"
line_t+=")"
line_o+=")"
line_w+=")"
line_b+=")"

print(line_s,file=ofile)
print(line_t,file=ofile)
print(line_o,file=ofile)
print(line_w,file=ofile)
print(line_b,file=ofile)

line_s = "seqs=( "
line_t = "mTs=( "
line_o = "OP_files=( "
line_w = "W_files=( "
line_b = "nbps=( "
for s in seqs_n15:
    line_s+=s.seq+" "
    line_t+=str(s.mT)+" "
    line_o+="op_n15.txt "
    line_w+="wfile_n15.txt "
    line_b+="15 "
line_s+=")"
line_t+=")"
line_o+=")"
line_w+=")"
line_b+=")"

print(line_s,file=ofile)
print(line_t,file=ofile)
print(line_o,file=ofile)
print(line_w,file=ofile)
print(line_b,file=ofile)

ofile.close()
