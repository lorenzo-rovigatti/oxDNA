
import sys

if len(sys.argv) != 5 :
    print("Unknown argument format.")
    print("Usage: python3 par_name par_value dim symm_type")
    print("dim = 2 or 4")
    print("sym_type = 0 for complementary symm (e.g. stacking)")
    print("sym_type = 1 for reverse symmetry (e.g. cross stack)")
    print("sym_type = 2 no symmetry")
    sys.exit()

parname = sys.argv[1]
parvalue = sys.argv[2]
dim = int(sys.argv[3])
sy_ty = int(sys.argv[4])

ofile = open("in_pars_value.txt", 'w')
ofile_opt = open("in_opt_par.txt", 'w')

bases = ['A','C','G','T']


jx_ind = []

if dim == 2 and sy_ty == 0:
    for i in range(4):
        for j in range(4):
            string = parname+"_"+bases[i]+"_"+bases[j]
            print(string+" = " + parvalue,file=ofile)
            if [3-j,3-i] in jx_ind : continue
            else:
                print("cname "+string, file=ofile_opt)
            jx_ind.append([i,j])
elif dim == 2 and sy_ty == 1:
    for i in range(4):
        for j in range(4):
            string = parname+"_"+bases[i]+"_"+bases[j]
            print(string+" = " + parvalue,file=ofile)
            if [j,i] in jx_ind : continue
            else:
                print("cname "+string, file=ofile_opt)
            jx_ind.append([i,j])
elif dim == 2 and sy_ty == 2:
    for i in range(4):
        for j in range(4):
            string = parname+"_"+bases[i]+"_"+bases[j]
            print(string+" = " + parvalue,file=ofile)
            print("cname "+string, file=ofile_opt)

elif dim == 4 and sy_ty == 0:
    for i in range(4):
        for j in range(4):
            for l in range(4):
                for m in range(4):
                    string = parname+"_"+bases[i]+"_"+bases[j]+"_"+bases[l]+"_"+bases[m]
                    print(string+" = " + parvalue,file=ofile)
                    if [3-m,3-l,3-j,3-i] in jx_ind : continue
                    else:
                        print("cname "+string, file=ofile_opt)
                    jx_ind.append([i,j,l,m])
elif dim == 4 and sy_ty == 1:
    for i in range(4):
        for j in range(4):
            for l in range(4):
                for m in range(4):
                    string = parname+"_"+bases[i]+"_"+bases[j]+"_"+bases[l]+"_"+bases[m]
                    print(string+" = " + parvalue,file=ofile)
                    if [m,l,j,i] in jx_ind : continue
                    else:
                        print("cname "+string, file=ofile_opt)
                    jx_ind.append([i,j,l,m])
elif dim == 4 and sy_ty == 2:
    for i in range(4):
        for j in range(4):
            for l in range(4):   
                for m in range(4):
                    string = parname+"_"+bases[i]+"_"+bases[j]+"_"+bases[l]+"_"+bases[m]
                    print(string+" = " + parvalue,file=ofile)
                    print("cname "+string, file=ofile_opt)

else:
	print("Invalid parameter rank.")

ofile.close()
ofile_opt.close()

