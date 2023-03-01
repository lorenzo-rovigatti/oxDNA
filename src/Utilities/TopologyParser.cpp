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

	char line[512];
	topology.getline(line, 512);

	if(sscanf(line, "%d %d\n", &_N, &_N_strands) == 2) {
		_parse_old_topology(topology, particles);
	}
	else {
		_parse_new_topology(topology, particles);
	}

	topology.close();

	if(N() != (int) particles.size()) {
		throw oxDNAException("Number of lines in the configuration file (%d) and\nnumber of particles in the topology file (%d) don't match. Aborting", particles.size(), N());
	}
}

void TopologyParser::_parse_old_topology(std::ifstream &topology, std::vector<BaseParticle *> &particles) {
	int N_from_conf = particles.size();
	char base[256], line[512];;
	int strand, i = 0;
	while(topology.good()) {
		topology.getline(line, 512);
		if(strlen(line) == 0 || line[0] == '#') {
			continue;
		}
		if(i == N_from_conf) {
			throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", N_from_conf);
		}

		int tmpn3, tmpn5;
		int res = sscanf(line, "%d %s %d %d", &strand, base, &tmpn3, &tmpn5);

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

	if(i < N_from_conf) {
		throw oxDNAException("Not enough particles found in the topology file (%d found, should be %d). Aborting", i, N_from_conf);
	}
}

void TopologyParser::_parse_new_topology(std::ifstream &topology, std::vector<BaseParticle *> &particles) {

}
