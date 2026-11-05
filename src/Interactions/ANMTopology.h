/*
 * ANMTopology.h
 *
 * Topology reader shared by the interactions that support network-model proteins (AC, DNANM/DNACT, RNANM/RNACT).
 * It reads both topology formats into a single Topology object, which is then used to set up the particles.
 */

#ifndef ANMTOPOLOGY_H_
#define ANMTOPOLOGY_H_

#include "../defs.h"

#include <map>
#include <string>
#include <vector>

class BaseParticle;

namespace ANMTopology {

/// one line of a new-style topology
struct Strand {
	bool peptide = false;
	bool circular = false;
	std::string sequence;
};

/// one particle line of an old-style topology
struct OldParticle {
	int strand; // strand id as written in the file: < 0 for amino acids (in mixed topologies), > 0 for nucleotides
	bool peptide;
	char code; // one-letter base or amino acid code
	int n3, n5; // neighbour indices as written in the file (nside and cside for amino acids), -1 if absent
	int line; // line number in the file, for error messages
};

struct Topology {
	bool new_style = false;
	int N = 0;
	int N_strands = 0;
	int N_nucleotides = 0;
	int N_amino_acids = 0;
	int N_nucleic_strands = 0;
	bool protein_first = false; // true if the first particle is an amino acid

	std::vector<Strand> strands; // filled for new-style topologies
	std::vector<OldParticle> old_particles; // filled for old-style topologies
};

/// true if the last token on the first line of the file is "5->3"
bool is_new_topology(const std::string &filename);

/**
 * @brief Reads an old- or new-style topology.
 *
 * nucleic_type ("DNA" or "RNA") is the type of the nucleic acids supported by the interaction. Pass an empty string
 * if the interaction supports proteins only (AC).
 *
 * New-style topologies have one strand per line, a sequence followed by key=value pairs; strands with type=peptide are
 * proteins. The header may contain the optional "<N_nucleotides> <N_amino_acids> <N_nucleic_strands>" integers before
 * "5->3", which are checked against the strands.
 *
 * Old-style topologies have one particle per line:
 *   - mixed topologies (nucleic_type not empty) have the header "<N> <N_strands> <N_nucleotides> <N_amino_acids> <N_nucleic_strands>",
 *     amino acids have negative strand ids and nucleotides positive ones;
 *   - protein-only topologies (nucleic_type empty) have the header "<N> <N_strands>" and every line is an amino acid.
 * Amino acid lines read "strand aa nside cside", nucleotide lines "strand base n3 n5". Any further column (the bonded
 * neighbours written by older tools) is ignored: protein bonds always come from the parameter file.
 */
Topology read(const std::string &filename, const std::string &nucleic_type);

/**
 * @brief Sets types, strand ids and connectivity of the particles from a Topology returned by read().
 *
 * Amino acids must land on ANMParticle objects and nucleotides on nucleotide objects, so the allocation must have been
 * done accordingly (amino acids and nucleotides must be in two separate blocks). Nucleic acid strands get strand ids
 * 0..N_nucleic_strands-1 and peptide strands the following ids.
 *
 * Amino acid sequences are read from the N to the C terminus, so the backbone neighbours of residue i are i-1
 * ("nside") and i+1 ("cside"), the convention used by the old-style topologies shipped with oxDNA (e.g. examples/DNANM/KDPG).
 *
 * @param protein_n3n5 if true, set n3 = nside and n5 = cside for amino acids (used by the DNACT/RNACT models)
 */
void fill_particles(const Topology &topology, std::vector<BaseParticle *> &particles, bool protein_n3n5);

/**
 * @brief Registers every spring of the parameter file as a bond between the two (protein) particles it connects.
 *
 * Pairs that are already bonded are skipped.
 */
void add_spring_bonds(const std::map<std::pair<int, int>, double> &springs, std::vector<BaseParticle *> &particles);

}

#endif /* ANMTOPOLOGY_H_ */
