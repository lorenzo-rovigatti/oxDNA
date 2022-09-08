#!/usr/bin/env python

# WORKS ONLY WITH AT MOST 2 REGIONS PER SCAFFOLD (details: competing staples as defined below will not be correctly counted)

# Measure assembly order parameter. This parameter is defined as follows:
# It measures the number of correctly bonded regions
# - a 'correctly bonded region' is defined as a continuous region between the staple and scaffold that has at least 8 correctly "hydrogen-bonded" nucleotides (condition is that there is at least 1 staple with 8 hb in that region
# - competing staples that have both formed 'correctly bonded regions' are only counted once
# - NB as the program is currently, a region with 8 correct hb from 2 different strands will be counted as being occupied by one of them in an unpredictable way (this is very unlikely to happen)

# the program requires a map of the correct bonds

# if "duplicate" option is enabled - assumes staples with identical sequences are copies and are therefore interchangeable when considering what counts as a correctly bonded nucleotide pair BUT IT PROBABLY DOESN'T WORK CORRECTLY
# otherwise it just counts the nucleotides it is given in the reference file as correctly bonded pairs
# NB the current procedure is to make a reference file with every possible correctly bonded pair (i.e. the duplication of staple strands is accounted for with this file). The file is generated using extract_ref.py 29/10/12

import readers
import base
import numpy as np
import sys

region_width = 16 # number of base pairs in a region
hb_threshold = 8 # this number of hbonds or more between a staple region and a scaffold region to count as attached

if len(sys.argv) < 5:
    base.Logger.log("Usage is %s configuration topology reference_h_bonds_file h_bonds_file [duplicates]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

base.Logger.log("Using region width %d, hydrogen bond count threshold %d" %(region_width, hb_threshold))
r = readers.LorenzoReader(sys.argv[1], sys.argv[2])
ss = r.get_system()
rawref = sys.argv[3]
rawdat = sys.argv[4]

check_duplicates = False
if len(sys.argv) > 5:
    if sys.argv[5] == "duplicates":
        base.Logger.log("checking for duplicates", base.Logger.INFO)
        check_duplicates = True

outfile = "assembly_op.dat"

lens = []
for strand in ss._strands:
    lens.append(strand.get_length())
scafid = lens.index(max(lens))
base.Logger.log("scaffold strand index %d" % scafid, base.Logger.INFO)
if (ss._strands[scafid].get_length() % region_width) != 0:
    base.Logger.log("scaffold length %d not commensurate with specified region width %d, dying" % (ss._strands[scafid].get_length(), region_width))
    sys.exit()

# get list of duplicate staple strands from topology file sequences
if check_duplicates:
    sequences = []
    duplicates = []
    first_nuc_id = [-1 for x in range(ss._N_strands)]
    for strand_id in range(len(ss._strands)):
        first_nuc_id[strand_id] = ss._strands[strand_id]._nucleotides[0].index
        if strand_id != scafid:
            strand = ss._strands[strand_id]
            sqs = strand.sequence
            # either start a new element in duplicates list
            if sqs not in sequences:
                sequences.append(sqs)
                duplicates.append([strand_id])
            # or add an index to an existing element
            else:
                this_id = sequences.index(sqs)
                duplicates[this_id].append(strand_id)

fref = open(rawref, "r")
fdat = open(rawdat, "r")
fout = open(outfile, "w")
    
base.Logger.log("writing to %s" % (outfile), base.Logger.INFO)

# get correct bonding data from reference file
partner = [-1 for x in range(ss._N)]
ss.map_nucleotides_to_strands()
scaf_nuc0id = ss._strands[scafid]._nucleotides[0].index
nuc2region = [-1 for x in range(ss._N)]

line = fref.readline()
while line != "":
    while line[0] == "#":
        line = fref.readline()

    while len(line) > 0 and line[0] != "#":
        cols = line.split(" ")
        cols = [x.strip() for x in cols]
        nuc1id = int(cols[1])
        nuc2id = int(cols[4])
        # partner[] records which nucleotide is bound to which
        partner[nuc1id] = nuc2id
        partner[nuc2id] = nuc1id

        # deal with staple duplicates
        if check_duplicates:
            base.Logger.log("guessing duplicates; partially developed feature that probably doesn't work", base.Logger.WARNING)
            strand_id = ss._nucleotide_to_strand[nuc1id]
            if strand_id != scafid:
                hits = 0
                for ii in duplicates:
                    if strand_id in ii:
                        dup_id = duplicates.index(ii)
                        hits += 1
                if hits != 1:
                    base.Logger.log("duplicate hits should be 1, but found %d hits, dying" % hits, base.Logger.CRITICAL)
                for strand_dup_id in duplicates[ii]:
                    nuc1id_dup = nuc1id + first_nuc_id[strand_dup_id] - first_nuc_id[strand_id]
                    partner[nuc1id_dup] = nuc2id 
                    partner[nuc2id] = nuc1id_dup # this bit is not sound; we are overwriting partners for the scaffold strand

            strand_id = ss._nucleotide_to_strand[nuc2id]
            if strand_id != scafid:
                hits = 0
                for ii in duplicates:
                    if strand_id in ii:
                        dup_id = duplicates.index(ii)
                        hits += 1
                if hits != 1:
                    base.Logger.log("duplicate hits should be 1, but found %d hits, dying" % hits, base.Logger.CRITICAL)
                for strand_dup_id in duplicates[ii]:
                    nuc2id_dup = nuc2id + first_nuc_id[strand_dup_id] - first_nuc_id[strand_id]
                    partner[nuc2id_dup] = nuc1id
                    partner[nuc1id] = nuc2id_dup

        # nuc2region maps nucleotides to regions
        if ss._nucleotide_to_strand[nuc1id] == scafid:
            scaf_nucid = nuc1id - scaf_nuc0id
        elif ss._nucleotide_to_strand[nuc2id] == scafid:
            scaf_nucid = nuc2id - scaf_nuc0id
        else:
            base.Logger.log("while reading reference bonds file: neither of partnered nucleotides %d, %d belong to scaffold, dying" % (nuc1id, nuc2id), base.Logger.CRITICAL)
            sys.exit(1)
        nuc2region[nuc1id] = int(np.floor(scaf_nucid/region_width))
        nuc2region[nuc2id] = int(np.floor(scaf_nucid/region_width))

        line = fref.readline()

# get region association map - verified to give correct answer for 7000 nuc excess assembly
region_assoc = []
for strand_id in range(ss.get_N_strands()):
    if strand_id != scafid:
        if ss._strands[strand_id].get_length() % region_width != 0:
            base.Logger.log("strand %d is not an integer number of region_width %d in length" % (strand_id, region_width), base.Logger.WARNING)
        for ii in range(int(ss._strands[strand_id].get_length()/region_width)):
            thisnuc = ss._strands[strand_id]._nucleotides[ii*region_width].index
            thisregion = nuc2region[thisnuc]
            hit = False
            for rr in region_assoc:
                if thisregion in rr:
                    hit = True
            if not hit:
                if ii == 0:
                    region_assoc.append([thisregion])
                else:
                    region_assoc[-1].append(thisregion)


line = fdat.readline()
conf_count = 0
while line != "":
    conf_count += 1
    base.Logger.log("read conf %d" % conf_count, base.Logger.INFO)
    regions = [0 for x in range(int(ss._strands[scafid].get_length()/region_width))] # check s._N/region_width integer !!!!!!!!!!!!
    correctly_bound_regions = 0
    double_count_check = []
    region_strand_counter = [{} for x in range(len(regions))]
    region_strand_id = [-1 for x in range(len(regions))]

    while line[0] == "#":
        line = fdat.readline()
        
    # we have reached a new configuration
    while len(line) > 0 and line[0] != "#":
        cols = line.split(" ")
        cols = [x.strip() for x in cols]
        # record the hydrogen bond
        nuc1id = int(cols[1])
        nuc2id = int(cols[4])
        #print nuc1id, nuc2id
        #print partner[nuc2id]
        if nuc1id == partner[nuc2id]:
            if (nuc1id not in double_count_check) and (nuc2id not in double_count_check):
                double_count_check.extend([nuc1id, nuc2id])
                this_region = nuc2region[nuc1id]
                regions[this_region] += 1
                # get the strand id for the staple strand
                strandid = ss._nucleotide_to_strand[ss._nucleotides[nuc1id].index]
                if strandid == scafid:
                    strandid = ss._nucleotide_to_strand[ss._nucleotides[nuc2id].index]
                if strandid not in region_strand_counter[this_region]:
                    #if this_region in [20, 21]:
                        #print nuc1id, this_region, strandid
                    region_strand_counter[this_region][strandid] = 1
                else:
                    region_strand_counter[this_region][strandid] += 1

        line = fdat.readline()

    # if a region has at least 1 staple strand with at least hb_threshold hydrogen bonds with the scaffold, that region counts as correctly bound
    for ii in range(len(region_strand_counter)):
        region = region_strand_counter[ii]
        for strand_id in iter(region):
            if region[strand_id] >= hb_threshold:
                region_strand_id[ii] = strand_id
                print ii, strand_id
                break

    for ra in region_assoc:
        strands_in_assoc_regions = []
        for ii in ra:
            if region_strand_id[ii] >= 0: # check whether there is actually a strand here
                strands_in_assoc_regions.append(region_strand_id[ii])
        # this part assumes at most 2 regions per staple
        if len(list(set(strands_in_assoc_regions))) > 1: # count number of unique staple ids in the associated regions
            #print list(set(strands_in_assoc_regions))
            base.Logger.log("found associated regions with competing staples", base.Logger.INFO)
            correctly_bound_regions += 1
        else:
            correctly_bound_regions += len(strands_in_assoc_regions)

    # write coarse grained info
    fout.write("%d %d\n" % (conf_count, correctly_bound_regions))
fref.close()
fdat.close()
fout.close()

base.Logger.log("read %d configurations" % conf_count, base.Logger.INFO)
