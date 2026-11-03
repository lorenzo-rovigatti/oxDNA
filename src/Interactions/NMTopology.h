/*
 * NMTopology.h
 *
 * Shared reader for new-style (5->3) topology files used by the interactions that support
 * network-model proteins (AC, DNANM/DNACT, RNANM/RNACT).
 */

#ifndef NMTOPOLOGY_H_
#define NMTOPOLOGY_H_

#include "../defs.h"

#include <map>
#include <string>
#include <vector>

class BaseParticle;

namespace NMTopology {

struct Strand {
	bool peptide = false;
	bool circular = false;
	std::string sequence;
};

/// true if the last token on the first line of the file is "5->3"
bool is_new_topology(const std::string &filename);

/**
 * @brief Parses a new-style topology file into strands.
 *
 * Each line is one strand: a sequence followed by key=value pairs. Strands with type=peptide are proteins,
 * all others are nucleic acids, whose type (if given) must be nucleic_type ("DNA" or "RNA"). Pass an empty
 * nucleic_type if the interaction supports proteins only.
 *
 * The particle counts are returned through the last three arguments. If the header contains the optional
 * "<N_nucleotides> <N_amino_acids> <N_nucleic_strands>" integers, they are checked against the strands.
 */
std::vector<Strand> read(const std::string &filename, const std::string &nucleic_type, int &N_nucleotides, int &N_amino_acids, int &N_nucleic_strands);

/**
 * @brief Sets types, strand ids and connectivity of the particles from the strands returned by read().
 *
 * Particles are assigned in line order, and peptide strands must land on ANMParticle objects (so allocation must
 * have been done accordingly). Nucleic acid strands get strand ids 0..N_nucleic_strands-1 in line order, peptide
 * strands get the following ids, matching the numbering used for old-style topologies.
 *
 * Nucleotides are connected 5'->3' along the line, as in DNAInteraction. Amino acid sequences are read from the N to
 * the C terminus, so the backbone neighbours of residue i are i-1 ("nside") and i+1 ("cside"), the convention used by
 * the old-style topologies shipped with oxDNA (e.g. examples/DNANM/KDPG). The interactions do not use the direction.
 *
 * @param protein_n3n5 if true, set n3 = nside and n5 = cside for amino acids (used by the DNACT/RNACT models)
 * @param protein_backbone_bonds if true, register backbone neighbours as bonded (used by AC and RNANM)
 */
void fill_particles(const std::vector<Strand> &strands, std::vector<BaseParticle *> &particles, int N_nucleic_strands, bool protein_n3n5, bool protein_backbone_bonds);

/**
 * @brief Registers every spring of the parameter file as a bond between the two (protein) particles it connects.
 *
 * Old-style topologies list these bonds on each amino acid line. New-style topologies do not, so the parameter
 * file is the only source. Pairs that are already bonded are skipped.
 */
void add_spring_bonds(const std::map<std::pair<int, int>, double> &springs, std::vector<BaseParticle *> &particles);

}

#endif /* NMTOPOLOGY_H_ */
