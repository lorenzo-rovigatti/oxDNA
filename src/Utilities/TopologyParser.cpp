/*
 * TopologyParser.cpp
 *
 *  Created on: Mar 1, 2023
 *      Author: lorenzo
 */

#include "TopologyParser.h"

#include "../Particles/BaseParticle.h"

#include <fstream>

TopologyParser::TopologyParser(std::string filename) :
				_filename(filename) {

}

TopologyParser::~TopologyParser() {

}

void TopologyParser::parse(std::vector<BaseParticle *> &particles) {
	std::ifstream topology;
	topology.open(_filename, std::ios::in);

	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", _filename.c_str());
	}

	std::string line;
	std::getline(topology, line);

	auto split = Utils::split(line, ' ');
	int N_from_topology = 0;
	if(split.size() == 2) {
		_N = std::atoi(split[0].c_str());
		_N_strands = std::atoi(split[1].c_str());
		N_from_topology = _parse_old_topology(topology, particles);
	}
	else if(split.size() == 3 && split[2] == "5->3") {
		_N = std::atoi(split[0].c_str());
		_N_strands = std::atoi(split[1].c_str());
		N_from_topology = _parse_new_topology(topology, particles);
	}
	else {
		throw oxDNAException("The first line of the topology file should contain two integers and a \"5->3\" token "
				"for new-style topology or two integers for old-style topology");
	}

	topology.close();

	if(N_from_topology != (int) particles.size()) {
		throw oxDNAException("The number of particles specified in the header of the topology file (%d) "
				"and the number of particles found in the topology (%d) don't match. Aborting",  N_from_topology, particles.size());
	}
}

int TopologyParser::_parse_old_topology(std::ifstream &topology, std::vector<BaseParticle *> &particles) {
	int N_from_conf = particles.size();
	char base[256];
	std::string line;
	int strand, i = 0;
	while(topology.good()) {
		std::getline(topology, line);
		if(line.size() == 0 || line[0] == '#') {
			continue;
		}
		if(i == N_from_conf) {
			throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", N_from_conf);
		}

		int tmpn3, tmpn5;
		int res = sscanf(line.c_str(), "%d %s %d %d", &strand, base, &tmpn3, &tmpn5);

		if(res < 4) {
			throw oxDNAException("Line %d of the topology file has an invalid syntax", i + 2);
		}

		BaseParticle *p = particles[i];

		if(tmpn3 < 0) {
			p->n3 = P_VIRTUAL;
		}
		else {
			p->n3 = particles[tmpn3];
		}
		if(tmpn5 < 0) {
			p->n5 = P_VIRTUAL;
		}
		else {
			p->n5 = particles[tmpn5];
		}

		// store the strand id
		// for a design inconsistency, in the topology file
		// strand ids start from 1, not from 0
		p->strand_id = strand - 1;

		// the base can be either a char or an integer
		if(strlen(base) == 1) {
			p->type = Utils::decode_base(base[0]);
			p->btype = Utils::decode_base(base[0]);
		}
		else {
			if(atoi(base) > 0) {
				p->type = atoi(base) % 4;
			}
			else {
				p->type = 3 - ((3 - atoi(base)) % 4);
			}
			p->btype = atoi(base);
		}

		if(p->type == P_INVALID) {
			throw oxDNAException("Particle #%d in strand #%d contains a non valid base '%c'. Aborting", i, strand, base);
		}
		p->index = i;
		i++;

		// here we fill the affected vector
		if(p->n3 != P_VIRTUAL) {
			p->affected.push_back(ParticlePair(p->n3, p));
		}

		if(p->n5 != P_VIRTUAL) {
			p->affected.push_back(ParticlePair(p, p->n5));
		}
	}

	return i;
}

std::vector<int> TopologyParser::_btypes_from_sequence(const std::string &sequence) const {
	std::vector<int> result;
	bool open_parenthesis = false;
	std::string parenthesis_token;

	for(char c : sequence) {
		if(c == '(') {
			open_parenthesis = true;
			parenthesis_token.clear();
			continue;
		}
		else if(c == ')') {
			if(!open_parenthesis) {
				throw oxDNAException("unbalanced parenthesis");
			}
			open_parenthesis = false;
			int btype = std::atoi(parenthesis_token.c_str());
			result.push_back(btype);
		}
		else {
			if(open_parenthesis) {
				parenthesis_token.push_back(c);
			}
			else {
				int btype = Utils::decode_base(c);
				if(btype == P_INVALID) {
					throw oxDNAException("invalid base %c", c);
				}
				result.push_back(btype);
			}
		}
	}

	if(open_parenthesis) {
		throw oxDNAException("missing closing parenthesis");
	}

	return result;
}

int TopologyParser::_parse_new_topology(std::ifstream &topology, std::vector<BaseParticle *> &particles) {
	int N_from_conf = particles.size();
	int current_idx = 0;
	for(int ns = 0; ns < _N_strands; ns++) {
		std::string line;
		std::getline(topology, line);
		auto split = Utils::split(line, ' ');

		std::vector<int> btypes;
		try {
			btypes = _btypes_from_sequence(split[0]);
		}
		catch(oxDNAException &e) {
			throw oxDNAException("topology file, strand %d (line %d): %s", ns, ns + 1, e.what());
		}

		int N_in_strand = btypes.size();
		for(int i = 0; i < N_in_strand; i++, current_idx++) {
			if(current_idx == N_from_conf) {
				throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", N_from_conf);
			}

			BaseParticle *p = particles[current_idx];
			p->strand_id = ns;
			p->btype = btypes[i];
			p->type = (p->btype < 0) ? 3 - ((3 - p->btype) % 4) : p->btype % 4;
			p->n3 = p->n5 = P_VIRTUAL;
			if(i != 0) {
				p->n5 = particles[current_idx - 1];
				p->affected.push_back(ParticlePair(p->n5, p));
			}
			if(i != N_in_strand - 1) {
				p->n3 = particles[current_idx + 1];
				p->affected.push_back(ParticlePair(p->n3, p));
			}
		}
	}

	return current_idx;
}
