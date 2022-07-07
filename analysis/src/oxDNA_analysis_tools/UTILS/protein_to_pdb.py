import Bio
import Bio.PDB
import numpy as np
import numpy.linalg as la
import copy


conv={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

                
def getalphacarbons_pdb(pdbid,pdbfile):
        structure = Bio.PDB.PDBParser().get_structure(pdbid, pdbfile)
        model = structure[0]
        # #chain index, residue index, residue identitiy, CA coordinates
        # cindx, rindx, rids, coord = [], [], [], []
        # if len(model) ==1:
        #     self.single_chain=True
        # else:
        #     self.single_chain=False
        calphacoord, ocoord, c_rids, c_atomids = [], [], [], []
        for cid, chain in enumerate(model):
            rids = []
            atomids = []
            res_indx =0
            for residue in chain.get_residues():
                rids.append(res_indx)
                tags= residue.get_full_id()
                if tags[3][0] == " ":
                    atoms=residue.get_atoms()
                    # Check for alpha carbon in residue
                    ac = [x for x in atoms if x.get_id() == 'CA']
                    if ac:
                        #Shouldn't have to recall this but won't work unless i do
                        atoms=residue.get_atoms()
                        calpha = ac[0].get_coord()
                        calphacoord.append(calpha)
                        for atom in atoms:
                            atomids.append(res_indx)
                            if atom.get_id() == 'CA':
                                ocoord.append(calpha)
                            else:
                                other = atom.get_coord()
                                ocoord.append(other) 
                        res_indx += 1   
            c_rids.append(rids)
            c_atomids.append(atomids) 
        npca = np.asarray(calphacoord)
        npoa = np.asarray(ocoord)
        #print(npca.shape, npoa.shape, len(rids))
        return npca, npoa, c_rids, c_atomids

def get_positions_oxfiles(datfile, conv=True):
    m = open(datfile, 'r')
    datdata = m.readlines()[3:]
    N = len(datdata)
    m.close()
    pos = []
    for line in datdata:
        if conv:
            ptmp = [8.518 * float(x) for x in line.split()[0:3] ]
        else:
            ptmp = [ float(x) for x in line.split()[0:2] ]
        pos.append(ptmp)
    nppos = np.asarray(pos)
    return nppos, N

# All units are Angstroms
# oca, N = get_positions_oxfiles('rotated.dat')
# pca, pfull, residue_ids = getalphacarbons_pdb('sgfp', 'sgfp.pdb')

def get_centroid(npcoords):
    centroid = np.mean(npcoords, axis=0)
    return centroid

# cs = get_centroid(oca)
# ct = get_centroid(pca)

def offset_coord(npcoords, centroid):
    off = np.subtract(npcoords, centroid)
    return off

#Adjust Coordinates to center at (0,0,0)
# oadj = offset_coord(oca, cs)
# padj = offset_coord(pca, ct)

# fulladj = offset_coord(pfull, ct)

# # Cross Covariance Matrix
# M = np.dot(padj.T,oadj)
# # Singular Value Decomposition
# u, s, v = la.svd(M, compute_uv=True)
# # Get Rotation Matrix
# R = np.dot(v.T, u.T)

def update_fullpdb_positions(R, fulladj, pcent, ocent):
    N = fulladj.shape[0]
    fullfit = np.full((N, 3), 0.0)
    for i in range(N):
        fullfit[i] = np.add(ocent, np.dot(R, fulladj[i]))
    return fullfit

def update_core_positions(R, padj, pcent, ocent):
    N = padj.shape[0]
    corefit = np.full((N, 3), 0.0)
    for i in range(N):
        corefit[i] = np.dot(R, padj[i])
    return corefit

# fullp = update_fullpdb_positions(R, fulladj, ct, cs)
# corep = update_core_positions(R, padj, ct, cs)

def slide_pdb_to_oxcoordinates(corep, oadj, fullp, rids):
    # rids is the residue each atom belongs to
    #corepositions are the unshifted calpha coordinates with centroid at (0,0,0)
    fulladj = np.full(fullp.shape, 0.0)
    # print(fulladj.shape)
    # print(fullp.shape)
    # print(oadj.shape, corep.shape)
    adj = np.subtract(oadj, corep)  # shift of c alpha position
    # print(adj)
    # print(adj.shape)
    
    for idx, i in enumerate(rids):
        # if i < 10:
        # print(i, idx)
        fulladj[idx] = fullp[idx] + adj[i]
    return fulladj

# print(residue_ids)
# fcoord = slide_pdb_to_oxcoordinates(corep, oadj, fullp, residue_ids)

def prep_pdb_sys(fullcoordinates, original_pdbfile, original_pdbid, old_reading_position, new_reading_position, offset=0, wipe=True):
    structure = Bio.PDB.PDBParser().get_structure(original_pdbid, original_pdbfile)
    model = structure[0]
    chain_total = len(model)
    # print('ctotal', chain_total)
    atm_idx = 0
    removal, chainremoval = [], []
    
    
    chainkeep  = []
    if new_reading_position == -1 and old_reading_position == 0:
        chainkeep += [x for x in range(chain_total)]  # entire file
    elif new_reading_position - old_reading_position == 1 or new_reading_position == -1: #just a specific chain
        chainkeep.append(old_reading_position)
    
    # print('chainkeep', chainkeep)
    for cid, chain in enumerate(model):
        # print(chain.get_id())
        if cid not in chainkeep:
            chainremoval.append(chain.get_id())
        else:
            for residue in chain.get_residues():
                tags = residue.get_full_id()
                if tags[3][0] == " ":
                    atoms = residue.get_atoms()
                    ac = [x for x in atoms if x.get_id() == 'CA']
                    if ac:
                        atoms=residue.get_atoms()
                        for atom in atoms:
                            atom.set_coord(list(fullcoordinates[atm_idx]))
                            atm_idx += 1
                    # elif wipe:
                    #     removal.append(residue.get_id())
                elif wipe:
                    removal.append(residue.get_id())
            if wipe:
                for x in removal:
                    chain.detach_child(x)
    if wipe:
        for x in chainremoval:
            model.detach_child(x)

    return model

def write_pdb(model, filehandle='', pdb_out=''):
    #Save File
    if pdb_out:
        io = Bio.PDB.PDBIO()
        io.set_structure(model)
        io.save(pdb_out+'.pdb')
    elif filehandle:
    # Returns PDB in Text for use in oxdna_to_PDB
        io = Bio.PDB.PDBIO()
        io.set_structure(model)
        io.save(filehandle, write_end=False)


# write_pdb(fcoord, 'sgfp.pdb', 'sgfp', 'sgfp_ox', wipe=True)
# ca, trash, trash2 =  getalphacarbons_pdb('sgfp_ox','sgfp_ox.pdb')
# print(ca[0], oca[0])



#For use in oxDNA_PDB.py Script
#will write when given the filehandle

def oxdna_to_pdb(filehandle, ox_protein_positions, pdbfile, nuccom, reading_position, offset=0):
    if "/" in pdbfile:
        pdbid=pdbfile.rsplit('/',1)[1].split('.')[0]
    else:
        pdbid=pdbfile.split('.')[0]
    # Get Coordinates
    oca = np.multiply(np.asarray(ox_protein_positions), 8.518)  # simulation to Angstroms
    pca, pfull, chain_separated_residue_ids, chain_separated_atom_ids = getalphacarbons_pdb(pdbid, pdbfile)
    
    #print('reading', reading_position)
    tmp = chain_separated_residue_ids[int(reading_position)]
    tmp_atom = chain_separated_atom_ids[int(reading_position)]
    #print('chains', len(chain_separated_residue_ids))
    
    res_ids = []
    new_reading_position = copy.deepcopy(reading_position)

    if reading_position >= len(chain_separated_residue_ids):
        print('trying to read data that doesnt exist')
        return
    elif reading_position == len(chain_separated_residue_ids)-1:
        abs_pos_res = sum([len(x) for x in chain_separated_residue_ids[:reading_position]])
        abs_pos_atom = sum([len(x) for x in chain_separated_atom_ids[:reading_position]])
        res_ids += chain_separated_atom_ids[reading_position]
        # print("ab res", abs_pos_res)
        # print("ab_atom", abs_pos_atom)
        pca = pca[abs_pos_res:abs_pos_res+len(tmp)]
        pfull = pfull[abs_pos_atom:abs_pos_atom+len(tmp_atom)]
        new_reading_position = -1 # done
    elif reading_position == 0:
        # print(len(tmp), len(oca))
        if sum([len(x) for x in chain_separated_residue_ids]) == len(oca):
             pca = pca
             pfull = pfull
             res_ids += [item for chain in chain_separated_atom_ids for item in chain]
             new_reading_position = -1
        elif len(tmp) == len(oca):
             # print("here")
             pca = pca[:len(tmp)]
             pfull = pfull[:len(tmp_atom)]
             res_ids += chain_separated_atom_ids[reading_position]
             new_reading_position += 1
        else:
             print('trying to read data that doesnt exist')
             return
    else:
        abs_pos_res = sum([len(x) for x in chain_separated_residue_ids[:reading_position]])
        abs_pos_atom = sum([len(x) for x in chain_separated_atom_ids[:reading_position]])
        # print("ab res", abs_pos_res)
        # print("ab_atom", abs_pos_atom)
        res_ids += chain_separated_atom_ids[reading_position]
        pca = pca[abs_pos_res:abs_pos_res+len(tmp)]
        pfull = pfull[abs_pos_atom:abs_pos_atom+len(tmp_atom)]
        new_reading_position += 1


    # pfull -> all atoms belonging to chain/model of PDB
    # pca -> CA atoms belonging to chain/model of PDB

    # Get Centroids
    cs = get_centroid(oca)   # COM of CA oxdna coordinates
    ct = get_centroid(pca)   # COM of CA DB coordinates
    
    # Center Centroid (0, 0, 0) in Coordinates for Rotation Matrix Calculation
    oadj = offset_coord(oca, cs)   # centered CA positions of oxDNA coordinates
    padj = offset_coord(pca, ct)   # centered CA positions of PDB coordinates
    
    #All atoms
    #print(pfull)
    fulladj = offset_coord(pfull, ct)
    # print('fulladj', fulladj)

    # Cross Covariance Matrix
    M = np.dot(padj.T,oadj)
    # Singular Value Decomposition
    u, s, v = la.svd(M, compute_uv=True)
    # Get Rotation Matrix
    R = np.dot(v.T, u.T)

    # Rotates and Centers all atoms of PDB to oxDNA Coordinates
    fullp = update_fullpdb_positions(R, fulladj, ct, cs)
    # Rotates and Centers ALphacarbons of PDB to oxDNA Coordinates
    corep = update_core_positions(R, padj, ct, cs)

    #Adjust position from PDB coordinate to oxCoordinate
    fcoord = slide_pdb_to_oxcoordinates(corep, oadj, fullp, res_ids)
    
    #Adjust position to line up with Nucleotides
    fcoord2 = np.subtract(fcoord, 8.518*nuccom)
    
    #Removes Everything from pdbfile but atoms of residues
    sys = prep_pdb_sys(fcoord2, pdbfile, pdbid, reading_position, new_reading_position, offset=offset)

    #write_pdb
    write_pdb(sys, filehandle)
    
    return new_reading_position


