/*
 * NMTopology.cpp
 */

#include "NMTopology.h"

#include "../Particles/BaseParticle.h"
#include "../Particles/ANMParticle.h"
#include "../Utilities/TopologyParser.h"
#include "../Utilities/Utils.h"

#include <fstream>

namespace NMTopology {

bool is_new_topology(const std::string &filename) {
	std::ifstream topology(filename, std::ios::in);
	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", filename.c_str());
	}
	std::string line;
	std::getline(topology, line);
	auto split = Utils::split(line, ' ');
	return split.size() > 0 && split.back() == "5->3";
}

std::vector<Strand> read(const std::string &filename, const std::string &nucleic_type, int &N_nucleotides, int &N_amino_acids, int &N_nucleic_strands) {
	TopologyParser parser(filename);
	parser.parse_new_topology();

	std::vector<Strand> strands;
	N_nucleotides = N_amino_acids = N_nucleic_strands = 0;
	int ns = 0;
	for(auto &seq_input : parser.lines()) {
		Strand strand;
		getInputString(&seq_input, "specs", strand.sequence, 1);
		getInputBool(&seq_input, "circular", &strand.circular, 0);

		std::string type;
		if(getInputString(&seq_input, "type", type, 0) == KEY_FOUND && type == "peptide") {
			strand.peptide = true;
			N_amino_acids += strand.sequence.size();
		}
		else {
			if(nucleic_type.empty()) {
				throw oxDNAException("topology file, strand %d (line %d): this interaction only supports type=peptide strands", ns, ns + 2);
			}
			if(!type.empty() && type != nucleic_type) {
				throw oxDNAException("topology file, strand %d (line %d): this interaction only supports type=peptide and type=%s strands, found type=%s", ns, ns + 2, nucleic_type.c_str(), type.c_str());
			}
			N_nucleotides += strand.sequence.size();
			N_nucleic_strands++;
		}

		strands.push_back(strand);
		ns++;
	}

	if(N_nucleotides + N_amino_acids != parser.N()) {
		throw oxDNAException("The topology header says there are %d particles, but the strands contain %d", parser.N(), N_nucleotides + N_amino_acids);
	}

	auto extras = parser.header_extras();
	if(extras.size() == 3) {
		if(extras[0] != N_nucleotides || extras[1] != N_amino_acids || extras[2] != N_nucleic_strands) {
			throw oxDNAException("The topology header lists %d nucleotides, %d amino acids and %d nucleic acid strands, but the strands contain %d, %d and %d",
					extras[0], extras[1], extras[2], N_nucleotides, N_amino_acids, N_nucleic_strands);
		}
	}
	else if(extras.size() != 0) {
		throw oxDNAException("Invalid topology header: expected '<N> <N_strands> 5->3' or '<N> <N_strands> <N_nucleotides> <N_amino_acids> <N_nucleic_strands> 5->3'");
	}

	return strands;
}

void fill_particles(const std::vector<Strand> &strands, std::vector<BaseParticle *> &particles, int N_nucleic_strands, bool protein_n3n5, bool protein_backbone_bonds) {
	// bonds are stored by index, so all indices have to be set before any bond is added
	for(uint i = 0; i < particles.size(); i++) {
		particles[i]->index = i;
	}

	int idx = 0;
	int nucleic_id = 0;
	int peptide_id = N_nucleic_strands;
	for(uint ns = 0; ns < strands.size(); ns++) {
		const Strand &strand = strands[ns];
		int N_in_strand = strand.sequence.size();
		int first = idx;

		for(int i = 0; i < N_in_strand; i++, idx++) {
			BaseParticle *p = particles[idx];
			char c = strand.sequence[i];
			bool is_anm = dynamic_cast<ANMParticle *>(p) != nullptr;
			p->index = idx;
			p->n3 = p->n5 = P_VIRTUAL;

			if(strand.peptide) {
				if(!is_anm) {
					throw oxDNAException("topology file, strand %d (line %d): amino acids and nucleotides must be listed in two separate blocks (all proteins first or all nucleic acids first)", ns, ns + 2);
				}
				p->type = p->btype = Utils::decode_aa(c);
				if(p->type == A_INVALID) {
					throw oxDNAException("topology file, strand %d (line %d): invalid amino acid '%c'", ns, ns + 2, c);
				}
				p->strand_id = peptide_id;
			}
			else {
				if(is_anm) {
					throw oxDNAException("topology file, strand %d (line %d): amino acids and nucleotides must be listed in two separate blocks (all proteins first or all nucleic acids first)", ns, ns + 2);
				}
				p->type = p->btype = Utils::decode_base(c);
				if(p->type == P_INVALID) {
					throw oxDNAException("topology file, strand %d (line %d): invalid base '%c'", ns, ns + 2, c);
				}
				p->strand_id = nucleic_id;

				// same connectivity as DNAInteraction::read_topology for new-style topologies
				if(i > 0) {
					p->n5 = particles[idx - 1];
					p->affected.push_back(ParticlePair(p->n5, p));
				}
				if(i < N_in_strand - 1) {
					p->n3 = particles[idx + 1];
					p->affected.push_back(ParticlePair(p->n3, p));
				}
				else if(strand.circular && N_in_strand > 1) {
					BaseParticle *first_p = particles[first];
					p->n3 = first_p;
					p->affected.push_back(ParticlePair(first_p, p));
					first_p->n5 = p;
					first_p->affected.push_back(ParticlePair(p, first_p));
				}
			}
		}

		if(strand.peptide) {
			// backbone: sequences are written from the N to the C terminus, so residue i has nside = i - 1 and cside = i + 1,
			// the convention used by the old-style topologies shipped with oxDNA (e.g. examples/DNANM/KDPG)
			int last = idx - 1;
			for(int i = first; i <= last; i++) {
				BaseParticle *nside = (i > first) ? particles[i - 1] : ((strand.circular && last > first) ? particles[last] : P_VIRTUAL);
				BaseParticle *cside = (i < last) ? particles[i + 1] : ((strand.circular && last > first) ? particles[first] : P_VIRTUAL);
				auto *p = static_cast<ANMParticle *>(particles[i]);
				if(protein_n3n5) {
					p->n3 = nside;
					p->n5 = cside;
				}
				if(protein_backbone_bonds && cside != P_VIRTUAL) {
					p->add_bonded_neighbor(cside);
				}
			}
			peptide_id++;
		}
		else {
			nucleic_id++;
		}
	}
}

void add_spring_bonds(const std::map<std::pair<int, int>, double> &springs, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	for(auto &spring : springs) {
		int i = spring.first.first;
		int j = spring.first.second;
		if(i < 0 || i >= N || j < 0 || j >= N || i == j) {
			throw oxDNAException("Invalid spring %d-%d in the parameter file", i, j);
		}
		auto *p = dynamic_cast<ANMParticle *>(particles[i]);
		if(p == nullptr || dynamic_cast<ANMParticle *>(particles[j]) == nullptr) {
			throw oxDNAException("The parameter file contains a spring between particles %d and %d, but only amino acids can be connected by springs", i, j);
		}
		p->add_bonded_neighbor(particles[j]); // does nothing if the two particles are already bonded
	}
}

}
