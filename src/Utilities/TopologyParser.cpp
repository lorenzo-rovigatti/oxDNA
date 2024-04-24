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
	std::ifstream topology(_filename, std::ios::in);

	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", _filename.c_str());
	}

	std::string line;
	std::getline(topology, line);

	auto split = Utils::split(line, ' ');
	if(split.size() == 2) {
		_N = std::atoi(split[0].c_str());
		_N_strands = std::atoi(split[1].c_str());
		_is_new_topology = false;
	}
	else if(split.size() == 3 && split[2] == "5->3") {
		_N = std::atoi(split[0].c_str());
		_N_strands = std::atoi(split[1].c_str());
		_is_new_topology = true;
	}
	else {
		throw oxDNAException("The first line of the topology file should contain two integers and a \"5->3\" token "
				"for new-style topology or two integers for old-style topology");
	}
	topology.close();
}

TopologyParser::~TopologyParser() {

}

bool TopologyParser::is_new_topology() {
	return _is_new_topology;
}

int TopologyParser::parse_old_topology(std::vector<BaseParticle *> &particles) {
	std::ifstream topology(_filename, std::ios::in);
	std::string line;
	std::getline(topology, line);

	char base[256];
	int strand, i = 0;
	while(topology.good()) {
		std::getline(topology, line);
		if(line.size() == 0 || line[0] == '#') {
			continue;
		}
		if(i == N()) {
			throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", N());
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

	topology.close();

	return i;
}

void TopologyParser::parse_new_topology() {
	std::ifstream topology(_filename, std::ios::in);
	std::string line;
	std::getline(topology, line);

	for(int ns = 0; ns < _N_strands; ns++) {
		std::getline(topology, line);

		if(!topology.good()) {
			throw oxDNAException("Not enough lines found in the topology file: %d instead of %d", ns, _N_strands);
		}

		auto split = Utils::split(line, ' ');

		std::string input_source = "specs=" + split[0] + "\n";
		for(uint i = 1; i < split.size(); i++) {
			input_source += split[i] + "\n";
		}
		input_file seq_input;
		seq_input.init_from_string(input_source);
		_lines.push_back(seq_input);
	}

	topology.close();
}
