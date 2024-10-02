#!/usr/bin/env python
import math
import sys, os, re
import argparse
'''
Implements the sequence-dependent SantaLucia model of DNA duplex melting temperatures,
as discussed in the paper:

SantaLucia, J. & Hicks, D. The thermodynamics of DNA structural motifs. Ann. Rev. Biophys. Biomol. Struct. 33, 415-440 (2004)

The values of Delta H and Delta S for internal mismatches are taken from references 1-4, 64 of the paper above. They are:
	# Biochemistry 1997, 36, 10581-10594 "Thermodynamics and NMR of Internal G,T Mismatches in DNA", by Hallawy and SantaLucia
	# Biochemistry 1998, 37, 9435-9444 "Nearest-Neighbor Thermodynamics of Internal A,C Mismatches in DNA: Sequence Dependence and pH Effects", by Allawy and SantaLucia
	# Biochemistry 1998, 37, 2170-2179 "Nearest Neighbor Thermodynamic Parameters for Internal G,A Mismatches in DNA", by Allawy and SantaLucia
	# Nucleic Acids Research, 1998, Vol. 26, No. 11 "Thermodynamics of internal C,T mismatches in DNA", by Allawy and SantaLucia
	# Biochemistry, Vol. 38, No. 12, 1999 "Nearest-Neighbor Thermodynamics and NMR of DNA Sequences with Internal A,A, C,C, G,G, and T,T Mismatches", by Peyret, Seneviratne, Allawy and SantaLucia

Output/input units are:
delta H			cal mol^-1
delta S			cal mol^-1 K^-1
T						K
salt_conc		molarity(mol L^-1)
box					see get_TM
'''
complementary_base = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
complementary_base.update({k.lower():v.lower() for k,v in complementary_base.iteritems()})
R_cal_mol_K = 1.9858775 

def pretty_print_duplex(seq, comp_seq, seq_is_53 = True):
	'''
	Returns a string containing a pretty-print a DNA duplex.
	
	seq: strand 1 sequence
	comp_seq: strand 2 sequence, in parallel order to strand 1, so that seq[0] is complementary to comp_seq[0], and so on.
	seq_is_53: if True, the order of the strand 1 is assumed to be 5'-3'. Otherwise it's 3'-5'. Defaults to True.
	'''
	mismatch_string = ''
	for i in range(len(seq)):
		if seq[i] != complementary_base[comp_seq[i]]: 
			mismatch_string += "-"
		else:
			mismatch_string += " "
	mismatch_mark = ''.join(mismatch_string)
	output_string = ''
	if args.seq_53 == None:
		output_string += "    "+mismatch_mark[::-1]+'\n'
		output_string += "3' -"+seq[::-1]+"- 5'"+'\n'
		output_string += "5' -"+comp_seq[::-1]+"- 3'"+'\n'
		output_string += "    "+mismatch_mark[::-1]
	else:
		output_string += "    "+mismatch_mark+'\n'
		output_string += "5' -"+seq+"- 3'"+'\n'
		output_string += "3' -"+comp_seq+"- 5'"+'\n'
		output_string += "    "+mismatch_mark
	return output_string


def complementary(seq, keep_direction = False):
	'''
	Returns a sequence complementary to the given sequence seq, and keeps the orientation if keep_direction is True (defaults to False).
	
	This means that if keep_direction is True, the complementary of AC will be GT.
	If keep_direction is false, the complementary of AC will be TG.
	
	
	'''
	comp = [complementary_base[i] for i in seq]
	if keep_direction: comp = comp[::-1]
	if type(seq) is str:
		comp = ''.join(comp)
	return comp
	
		

def is_self_complementary(seq):
	'''
	Returns true if seq is a self-complementary DNA sequence, and returns false otherwise.
	'''
	length = len(seq)
	# assumes that a sequence is self-complementary, and checks wheter that's actually the case.
	self_complementary = True
	for i in range(length):
		i_c = length - 1 - i
		if seq[i] != complementary_base[seq[i_c]]:
			self_complementary = False
			break

	return self_complementary

oxDNA_box_unit = 0.8518e-9 #box unit in meters, from the oxDNA wiki
N_Avogadro = 6.022140857e23 #Avogadro costant
def get_box_from_duplex_concentration(duplex_concentration, box_unit = oxDNA_box_unit):
	'''
	Given the duplex concentration, return the length of the side of the cubic box in which to simulate one duplex.
	'''
	conc_factor = (0.1 / box_unit)**3 / N_Avogadro # the factor of 0.1 comes from the fact that molarity is mol/L
	box_side = (conc_factor / duplex_concentration)**(1./3.)
	return box_side

def get_duplex_concentration_from_box(box_size, box_unit = oxDNA_box_unit):
	'''
	Given the side of the cubic box in oxDNA units (or any other, if box_unit is changed) and return the duplex concentration.
	'''
	conc_factor = (0.1 / box_unit)**3 / N_Avogadro # the factor of 0.1 comes from the fact that molarity is mol/L
	duplex_concentration = conc_factor / box_size ** 3
	# alternate, old way to compute duplex molar concentration
	# (yields a slightly different value)
	# I'm keeping it here for future debugs, just in case
	alt_duplex_concentration = 2.6868994 / (box_size*box_size*box_size) ;
	return duplex_concentration

def get_TM_from_H_S_duplexC(H,S,duplexC, strand_is_self_complementary = False):
	'''
	Get the duplex melting temperature from provided DeltaH, DeltaS for a given duplex concentration.
	
	The DeltaH and DeltaS should all be given in cal/mol (or cal/mol K), and T will be in degrees Celsius.
	If the strand is self-complementary, the third argument should be set to True.
	'''
	x = 4 - 3 * strand_is_self_complementary
	return H / (S + R_cal_mol_K * math.log( 2 * duplexC / x)) - 273.15
	
	
def get_TM(seq,boxsize = 16.4,salt_conc = 0.5,seq_comp = [],duplex_concentration = None, phosphate_correction = 0.):
	'''
	Computes melting temperature, DeltaH and DeltaS of a DNA duplex.

	seq: strand 1 sequence (in the 5'-3' order)
	boxsize: size of the simulation box in oxDNA units (this XOR duplex_concentration must be provided)
	salt_conc: molar concentration of [NaCl]
	seq_comp: strand 2 sequence (in the 3'-5' order). Defaults to the complementary sequence of seq.
	duplex_concentration: molar concentration of the DNA duplex (this XOR boxsize must be provided)
	excess_phosphate: number of phosphate pairs in excess of the number of base-pairs, important for the salt correction.
												in oxDNA 2 there is half a phosphate on each chain end, so this is 0 (default).
												that's also the case for phosphates only on the 5' or 3' end.
												a duplex with no phosphates on either end has -1 phosphate pairs.
												a duplex with phosphates on all ends has +1 phosphate pairs.
	'''
	
	if duplex_concentration == None:
		duplex_concentration = get_duplex_concentration_from_box(boxsize)

	length = len(seq);
	deltaH = 0;
	deltaS = 0;
	# dictionary for Delta Enthalpy and Delta Enthropy - written here in a way that's easy
  # to compare with the paper
	deltaHS = {}
# insert the values from Table 1, of Annu. Rev. Biophys. Biomol. Struct. 2004. 33:415-40
# "The Thermodynamics of DNA Structural Motifs", by SantaLucia and Hicks
	deltaHS.update({
	'AA/TT':( -7.6,  -21.3),
	'AT/TA':( -7.2,  -20.4),
	'TA/AT':( -7.2,  -21.3),
	'CA/GT':( -8.5,  -22.7),
	'GT/CA':( -8.4,  -22.4),
	'CT/GA':( -7.8,  -21.0),
	'GA/CT':( -8.2,  -22.2),
	'CG/GC':( -10.6, -27.2),
	'GC/CG':( -9.8,  -24.4),
	'GG/CC':( -8.0,  -19.9)
	})
	# insert the values from Table 5, of Biochemistry 1997, 36, 10581-10594
	# "Thermodynamics and NMR of Internal G,T Mismatches in DNA", by Hallawy and SantaLucia
	deltaHS.update({
	'AG/TT': ( 1.0,  0.9),
	'AT/TG': (-2.5, -8.3),
	'CG/GT': (-4.1, -11.7),
	'CT/GG': (-2.8, -8.0),
	'GG/CT': ( 3.3,  10.4),
	'GG/TT': ( 5.8,  16.3),
	'GT/CG': (-4.4, -12.3),
	'GT/TG': ( 4.1,  9.5),
	'TG/AT': (-0.1, -1.7),
	'TG/GT': (-1.4, -6.2),
	'TT/AG': (-1.3, -5.3)
	})
	# insert the values from Biochemistry 1998, 37, 9435-9444
	# "Nearest-Neighbor Thermodynamics of Internal A,C Mismatches in DNA: Sequence Dependence and pH Effects", by Allawy and SantaLucia
	# WARNING: they report "the contribution of single A,C
	#mismatches to duplex stability is strongly dependent on the solution pH and the nearest-neighbor context", which is at variance with what reported for the other mismatches. I'm using the pH 7.0 data, but the rest is also stored here for future use.
	deltaHS_ph7 = {
	'AA/TC' :( 2.3,  4.6),
	'AC/TA' :( 5.3,  14.6),
	'CA/GC' :( 1.9,  3.7),
	'CC/GA' :( 0.6, -0.6),
	'GA/CC' :( 5.2,  14.2),
	'GC/CA' :(-0.7, -3.8),
	'TA/AC' :( 3.4,  8.0),
	'TC/AA' :( 7.6,  20.2)
	}
	deltaHS_ph5 = {
	'AA/TC' :(-0.8, -3.8),
	'AC/TA' :(-6.3, -20.2),
	'CA/GC' :(-4.2, -13.6),
	'CC/GA' :(-1.3, -4.9),
	'GA/CC' :(-3.3, -10.3),
	'GC/CA' :(-4.9, -14.7),
	'TA/AC' :(-2.1, -7.6),
	'TC/AA' :( 2.2,  4.8)
	}
	deltaHS.update(deltaHS_ph7)
	# insert the values from Biochemistry 1998, 37, 2170-2179
	# "Nearest Neighbor Thermodynamic Parameters for Internal G,A Mismatches in DNA", by Allawy and SantaLucia
	deltaHS.update({
	'AA/TG': (-0.6, -2.3),
	'AG/TA': (-0.7, -2.3),
	'CA/GG': (-0.7, -2.3),
	'CG/GA': (-4.0, -13.2),
	'GA/CG': (-0.6, -1.0),
	'GG/CA': ( 0.5,  3.2),
	'TA/AG': ( 0.7,  0.7),
	'TG/AA': ( 3.0,  7.4)
	})
	# insert the values from Table 4, of Nucleic Acids Research, 1998, Vol. 26, No. 11
	# "Thermodynamics of internal C,T mismatches in DNA", by Allawy and SantaLucia
	deltaHS.update({
	'AC/TT': ( 0.7,  0.2),
	'AT/TC': (-1.2, -6.2),
	'CC/GT': (-0.8, -4.5),
	'CT/GC': (-1.5, -6.1),
	'GC/CT': ( 2.3,  5.4),
	'GT/CC': ( 5.2,  13.5),
	'TC/AT': ( 1.2,  0.7),
	'TT/AC': ( 1.0,  0.7)
	})
	# insert the values from Table 2, of Biochemistry, Vol. 38, No. 12, 1999
	# "Nearest-Neighbor Thermodynamics and NMR of DNA Sequences with Internal A,A, C,C, G,G, and T,T Mismatches", by Peyret, Seneviratne, Allawy and SantaLucia
	# notice that the error bars (not saved here) are pretty big. See original text for more info.
	deltaHS.update({
	'AA/TA' :( 1.2,  1.7),
	'CA/GA' :(-0.9, -4.2),
	'GA/CA' :(-2.9, -9.8),
	'TA/AA' :( 4.7,  12.9),
	'AC/TC' :( 0.0, -4.4),
	'CC/GC' :(-1.5, -7.2),
	'GC/CC' :( 3.6,  8.9),
	'TC/AC' :( 6.1,  16.4),
	'AG/TG' :(-3.1, -9.5),
	'CG/GG' :(-4.9, -15.3),
	'GG/CG' :(-6.0, -15.8),
	'TG/AG' :( 1.6,  3.6),
	'AT/TT' :(-2.7, -10.8),
	'CT/GT' :(-5.0, -15.8),
	'GT/CT' :(-2.2, -8.4),
	'TT/AT' :( 0.2, -1.5)
	})
	# add every complementary base-pair-step
	deltaHS_compl = {}
	for k,v in deltaHS.iteritems():
		if not k[::-1] in deltaHS.keys():
			deltaHS_compl[k[::-1]] = v
		elif k[::-1] != k:
			print 'ERROR: double addition of key',k
	deltaHS.update(deltaHS_compl)

	# check that the two sequences have the same length
	if seq_comp == []:
		seq_comp = complementary(seq)
	if len(seq) != len(seq_comp):
		print 'ERROR: function',__name__,'called with seq of length',len(seq),'and seq_comp of length',len(seq_comp),'but the two sequences should have the same length'
		sys.exit(-1)
	# check that there are no terminal mismatches, as the parameters here are only for internal mismatches
	exit = False
	if seq_comp[0]  != complementary(seq[0]):
		print 'ERROR: sequences in get_TM have terminal mismatches, which are not implemented in the program.'
		print "       strand 1 has",seq[0],"on the 5' end, but strand 2 has",seq_comp[0],"on the 3' end."
		exit = True
	if seq_comp[-1] != complementary(seq[-1]):
		print 'ERROR: sequences in get_TM have terminal mismatches, which are not implemented in the program.'
		print "       strand 1 has",seq[-1],"on the 3' end, but strand 2 has",seq_comp[-1],"on the 5' end."
		exit = True
	# check that there are no 2 consecutive mismatches
	for i in range(len(seq)):
		if seq[i] != complementary(seq_comp[i]) and seq[i+1] != complementary(seq_comp[i+1]):
			print "ERROR: the duplex contains two or more consecutive mismatches. Parameters for that (if they exist) are not implemented."
			exit = True
			break
	if exit: sys.exit(-1)
	
	for i in range (length-1) :
		my_id = seq[i:i+2] + '/'+ seq_comp[i:i+2]
		
		deltaH += deltaHS[my_id][0]
		deltaS += deltaHS[my_id][1]
		
	# add penalty for an AT sequence end
	n_terminal_AT = (seq[0] in ['A','T']) + (seq[-1] in ['A','T'])
	deltaH += 2.2 * n_terminal_AT
	deltaS += 6.9 * n_terminal_AT
	 

	deltaH += 0.2;
	deltaS += -5.6;
	# convert deltaH from kcal/mol to cal/mol
	deltaH *= 1000.0;

	divisor = 2.0;
 	# self-complementarity correction
	if(is_self_complementary(seq)):
		deltaS += -1.4;
		divisor = 0.5;
	# salt entropy correction
	salt_correction = 0.368 * (len(seq) + phosphate_correction) * math.log(salt_conc)
	deltaS += salt_correction
 


	Tm = deltaH / (deltaS + R_cal_mol_K * math.log(duplex_concentration/divisor));
	return Tm - 273.15, deltaH, deltaS

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Compute the SantaLucia melting temperature of a DNA "+ 
	"duplex of a given sequence. The sequence can be given either as read from the 5' end to the 3' end, "+
	"as is convention in the research community, or in the 3' to 5' order, as in the oxDNA program.")
	# argumets for strand 1 sequence
	seq_group = parser.add_mutually_exclusive_group(required=True)
	seq_group.add_argument('--seq_53', type = str, help = "the sequence of one of the strands, from the 5' end to the 3' end.")
	seq_group.add_argument('--seq_oxdna', type = str, help = "the sequence of one of the strands, from the 3' end to the 5' end (as is convention in the oxDNA program).")
	# arguments for strand 2 sequence
	comp_seq_group = parser.add_mutually_exclusive_group()
	comp_seq_group.add_argument('--comp_seq_53', type = str, help = "the sequence of the complementary strand, from the 5' end to the 3' end.")
	comp_seq_group.add_argument('--comp_seq_oxdna', type = str, help = "the sequence of the complementary strand, from the 3' end to the 5' end (as is convention in the oxDNA program).")
	# arguments for duplex concentration
	duplex_c_group = parser.add_mutually_exclusive_group(required = True)
	duplex_c_group.add_argument('--box', type = float, help = "the length of the edge of the cubic simulation box, in oxDNA units.")
	duplex_c_group.add_argument('--duplex_conc', type = float, help = "the molar concentration of the duplex.")
	# salt concentration
	parser.add_argument('--salt',type = float, help = "[NaCl] molar concentration",default = 0.5)
	# number of phosphate pairs in excess of the number of base-pairs (0 in oxDNA2)
	parser.add_argument('--extra_phosphate', type = int, help = "number of phosphate pairs in excess of the number of base-pairs in the duplex (0 for oxDNA2 or for a duplex with phosphates only on the 3' OR 5' end, 1 for duplex with phosphate at every chain end, -1 for a duplex for no phosphate at any chaing end)",default = 0)
	# whether to print only T_M
	parser.add_argument('-T','--only-T', action = 'store_true', help = 'just print the melting temperature as opposed to printing all the other useful informations (useful when calling the scrip a bazillion of times with different parameters)')
	# print help as needed
	if len(sys.argv) == 1:
		parser.print_help()
		parser.exit()

 	 
	# handle arguments
	args = parser.parse_args()
	salt = args.salt
	box = args.box
	duplex_concentration = args.duplex_conc
	phosphate_correction = args.extra_phosphate
	
	# primary sequence is always treated in the 5'-3' order internally
	# and the secondary sequence is always treated in the 3'-5' order internally
	if args.seq_53 != None:
		seq = args.seq_53.upper()
	else:
		seq = args.seq_oxdna[::-1].upper()

	if args.comp_seq_53 != None:
		comp_seq = args.comp_seq_53[::-1].upper()
	elif args.comp_seq_oxdna != None:
		comp_seq = args.comp_seq_oxdna.upper()
	else:
		comp_seq = complementary(seq)
	# get the melting temperature, DeltaH and DeltaS
	print_only_TM = args.only_T
	T_M, DeltaH, DeltaS = get_TM(seq, box, salt, comp_seq, duplex_concentration, phosphate_correction)
	if not print_only_TM:
		# pretty print the duplex sequence, to avoid mistakes
		print pretty_print_duplex(seq, comp_seq, seq_is_53 = args.seq_53 != None )
		# print the box size if duplex_conc is given, or duplex_conc if box is given
		if box == None:
			print 'box size = ',get_box_from_duplex_concentration(duplex_concentration),'oxDNA unit lengths'
		else:
			print 'duplex concentration = ',get_duplex_concentration_from_box(box),'M'

		print 'DeltaH\t'+str(DeltaH)+'\tcal mol^-1'
		print 'DeltaS\t'+str(DeltaS)+'\tcal mol^-1 K^-1'
		print 'T_M\t'+str(T_M)+'\tC'
	else:
		print "%.1f"%(T_M)

