/*
 * ANMTopology.cpp
 *  Created Nov 2026
 *  Author: erik
 *  Implements topology readers for the ANM models
 */

#include "ANMTopology.h"

#include "../Particles/BaseParticle.h"
#include "../Particles/ANMParticle.h"
#include "TopologyParser.h"
#include "Utils.h"

#include <fstream>
#include <set>
#include <sstream>

namespace ANMTopology {

namespace {

const char *INTERLEAVED_ERROR = "topology file: amino acids and nucleotides must be listed in two separate blocks (all proteins first or all nucleic acids first)";

Topology read_new(const std::string &filename, const std::string &nucleic_type) {
	TopologyParser parser(filename);
	parser.parse_new_topology();

	Topology top;
	top.new_style = true;
	int ns = 0;
	for(auto &seq_input : parser.lines()) {
		Strand strand;
		getInputString(&seq_input, "specs", strand.sequence, 1);
		getInputBool(&seq_input, "circular", &strand.circular, 0);

		std::string type;
		if(getInputString(&seq_input, "type", type, 0) == KEY_FOUND && type == "peptide") {
			strand.peptide = true;
			top.N_amino_acids += strand.sequence.size();
		}
		else {
			if(nucleic_type.empty()) {
				throw oxDNAException("topology file, strand %d (line %d): this interaction only supports type=peptide strands", ns, ns + 2);
			}
			if(!type.empty() && type != nucleic_type) {
				throw oxDNAException("topology file, strand %d (line %d): this interaction only supports type=peptide and type=%s strands, found type=%s", ns, ns + 2, nucleic_type.c_str(), type.c_str());
			}
			top.N_nucleotides += strand.sequence.size();
			top.N_nucleic_strands++;
		}

		top.strands.push_back(strand);
		ns++;
	}

	top.N = parser.N();
	top.N_strands = top.strands.size();
	top.protein_first = top.strands.size() > 0 && top.strands[0].peptide;

	if(top.N_nucleotides + top.N_amino_acids != top.N) {
		throw oxDNAException("The topology header says there are %d particles, but the strands contain %d", top.N, top.N_nucleotides + top.N_amino_acids);
	}

	auto extras = parser.header_extras();
	if(extras.size() == 3) {
		if(extras[0] != top.N_nucleotides || extras[1] != top.N_amino_acids || extras[2] != top.N_nucleic_strands) {
			throw oxDNAException("The topology header lists %d nucleotides, %d amino acids and %d nucleic acid strands, but the strands contain %d, %d and %d",
					extras[0], extras[1], extras[2], top.N_nucleotides, top.N_amino_acids, top.N_nucleic_strands);
		}
	}
	else if(extras.size() != 0) {
		throw oxDNAException("Invalid topology header: expected '<N> <N_strands> 5->3' or '<N> <N_strands> <N_nucleotides> <N_amino_acids> <N_nucleic_strands> 5->3'");
	}

	return top;
}

Topology read_old(const std::string &filename, const std::string &nucleic_type) {
	bool mixed = !nucleic_type.empty();

	std::ifstream file(filename, std::ios::in);
	if(!file.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", filename.c_str());
	}

	Topology top;
	std::string line;
	std::getline(file, line);
	std::stringstream head(line);
	head >> top.N >> top.N_strands;
	int header_nucleotides = 0, header_amino_acids = 0, header_nucleic_strands = 0;
	if(mixed) {
		head >> header_nucleotides >> header_amino_acids >> header_nucleic_strands;
		bool invalid = head.fail() || top.N < 0 || top.N_strands < 0 || top.N_strands > top.N;
		invalid = invalid || header_nucleotides < 0 || header_nucleotides > top.N || header_amino_acids < 0 || header_amino_acids > top.N;
		invalid = invalid || header_nucleic_strands < 0 || header_nucleic_strands > top.N;
		if(invalid) {
			throw oxDNAException("Invalid topology header: old-style topologies for this interaction must start with '<N> <N_strands> <N_nucleotides> <N_amino_acids> <N_nucleic_strands>'");
		}
	}
	else if(head.fail() || top.N < 0 || top.N_strands < 0 || top.N_strands > top.N) {
		throw oxDNAException("Invalid topology header: old-style protein topologies must start with '<N> <N_strands>'");
	}

	std::set<int> nucleic_strand_ids;
	int line_number = 1;
	int last_strand = 0;
	while(std::getline(file, line)) {
		line_number++;
		if(line.find_first_not_of(" \t\r") == std::string::npos || line[0] == '#') {
			continue;
		}
		if((int) top.old_particles.size() == top.N) {
			throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", top.N);
		}

		// anything after the fourth column (e.g. the list of bonded neighbours of amino acids) is ignored
		std::stringstream ss(line);
		OldParticle p;
		std::string code;
		ss >> p.strand >> code >> p.n3 >> p.n5;
		if(ss.fail()) {
			throw oxDNAException("topology file, line %d: invalid syntax, expected '<strand> <base or amino acid> <n3 or nside> <n5 or cside>'", line_number);
		}
		if(p.strand == 0) {
			throw oxDNAException("topology file, line %d: strand ids cannot be 0", line_number);
		}
		p.peptide = !mixed || p.strand < 0;
		if(code.size() != 1) {
			if(p.peptide) {
				throw oxDNAException("topology file, line %d: invalid amino acid '%s'", line_number, code.c_str());
			}
			throw oxDNAException("topology file, line %d: only %s base characters are permitted, found '%s'", line_number, nucleic_type.c_str(), code.c_str());
		}
		p.code = code[0];
		if(p.n3 >= top.N || p.n5 >= top.N) {
			throw oxDNAException("topology file, line %d: invalid neighbour index", line_number);
		}

		if(p.peptide) {
			top.N_amino_acids++;
		}
		else {
			top.N_nucleotides++;
			nucleic_strand_ids.insert(p.strand);
		}
		p.line = line_number;
		last_strand = p.strand;
		top.old_particles.push_back(p);
	}

	if((int) top.old_particles.size() < top.N) {
		throw oxDNAException("Not enough particles found in the topology file (should be %d). Aborting", top.N);
	}

	top.N_nucleic_strands = nucleic_strand_ids.size();
	top.protein_first = top.N > 0 && top.old_particles[0].peptide;

	if(mixed) {
		if(header_nucleotides != top.N_nucleotides || header_amino_acids != top.N_amino_acids || header_nucleic_strands != top.N_nucleic_strands) {
			throw oxDNAException("The topology header lists %d nucleotides, %d amino acids and %d nucleic acid strands, but the file contains %d, %d and %d",
					header_nucleotides, header_amino_acids, header_nucleic_strands, top.N_nucleotides, top.N_amino_acids, top.N_nucleic_strands);
		}
	}
	else if(std::abs(last_strand) != top.N_strands) {
		throw oxDNAException("Mismatching strand number in header of Topology file");
	}

	return top;
}

void fill_particles_new(const Topology &top, std::vector<BaseParticle *> &particles, bool protein_n3n5) {
	int idx = 0;
	int nucleic_id = 0;
	int peptide_id = top.N_nucleic_strands;
	for(uint ns = 0; ns < top.strands.size(); ns++) {
		const Strand &strand = top.strands[ns];
		int N_in_strand = strand.sequence.size();
		int first = idx;

		for(int i = 0; i < N_in_strand; i++, idx++) {
			BaseParticle *p = particles[idx];
			char c = strand.sequence[i];
			bool is_anm = dynamic_cast<ANMParticle *>(p) != nullptr;
			p->n3 = p->n5 = P_VIRTUAL;

			if(strand.peptide) {
				if(!is_anm) {
					throw oxDNAException(INTERLEAVED_ERROR);
				}
				p->type = p->btype = Utils::decode_aa(c);
				if(p->type == A_INVALID) {
					throw oxDNAException("topology file, strand %d (line %d): invalid amino acid '%c'", ns, ns + 2, c);
				}
				p->strand_id = peptide_id;
			}
			else {
				if(is_anm) {
					throw oxDNAException(INTERLEAVED_ERROR);
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
			// sequences are written from the N to the C terminus, so residue i has nside = i - 1 and cside = i + 1
			int last = idx - 1;
			if(protein_n3n5) {
				for(int i = first; i <= last; i++) {
					bool closed = strand.circular && last > first;
					particles[i]->n3 = (i > first) ? particles[i - 1] : (closed ? particles[last] : P_VIRTUAL);
					particles[i]->n5 = (i < last) ? particles[i + 1] : (closed ? particles[first] : P_VIRTUAL);
				}
			}
			peptide_id++;
		}
		else {
			nucleic_id++;
		}
	}
}

void fill_particles_old(const Topology &top, std::vector<BaseParticle *> &particles, bool protein_n3n5) {
	auto neighbour = [&particles](int idx) {
		return (idx < 0) ? P_VIRTUAL : particles[idx];
	};

	for(uint i = 0; i < particles.size(); i++) {
		const OldParticle &rec = top.old_particles[i];
		BaseParticle *p = particles[i];
		bool is_anm = dynamic_cast<ANMParticle *>(p) != nullptr;
		p->n3 = p->n5 = P_VIRTUAL;

		if(rec.peptide) {
			if(!is_anm) {
				throw oxDNAException(INTERLEAVED_ERROR);
			}
			p->type = p->btype = Utils::decode_aa(rec.code);
			if(p->type == A_INVALID) {
				throw oxDNAException("topology file, line %d: invalid amino acid '%c'", rec.line, rec.code);
			}
			// protein strand ids are negative and start from -1 in mixed topologies, positive and start from 1 in protein-only ones
			p->strand_id = std::abs(rec.strand) + top.N_nucleic_strands - 1;
			if(protein_n3n5) {
				p->n3 = neighbour(rec.n3);
				p->n5 = neighbour(rec.n5);
			}
		}
		else {
			if(is_anm) {
				throw oxDNAException(INTERLEAVED_ERROR);
			}
			p->type = p->btype = Utils::decode_base(rec.code);
			if(p->type == P_INVALID) {
				throw oxDNAException("topology file, line %d: invalid base '%c'", rec.line, rec.code);
			}
			p->strand_id = rec.strand - 1;
			p->n3 = neighbour(rec.n3);
			p->n5 = neighbour(rec.n5);
			if(p->n3 != P_VIRTUAL) p->affected.emplace_back(ParticlePair(p->n3, p));
			if(p->n5 != P_VIRTUAL) p->affected.emplace_back(ParticlePair(p, p->n5));
		}
	}
}

}

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

Topology read(const std::string &filename, const std::string &nucleic_type) {
	if(is_new_topology(filename)) {
		return read_new(filename, nucleic_type);
	}
	return read_old(filename, nucleic_type);
}

void fill_particles(const Topology &top, std::vector<BaseParticle *> &particles, bool protein_n3n5) {
	if((int) particles.size() != top.N) {
		throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match. Aborting");
	}
	// bonds are stored by index, so all indices have to be set before any bond is added
	for(uint i = 0; i < particles.size(); i++) {
		particles[i]->index = i;
	}

	if(top.new_style) {
		fill_particles_new(top, particles, protein_n3n5);
	}
	else {
		fill_particles_old(top, particles, protein_n3n5);
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
