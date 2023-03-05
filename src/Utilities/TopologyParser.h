/*
 * TopologyParser.h
 *
 *  Created on: Mar 1, 2023
 *      Author: lorenzo
 */

#ifndef SRC_UTILITIES_TOPOLOGYPARSER_H_
#define SRC_UTILITIES_TOPOLOGYPARSER_H_

#include "../defs.h"

class BaseParticle;

class TopologyParser {
public:
	TopologyParser(std::string filename);
	virtual ~TopologyParser();
	TopologyParser(const TopologyParser &other) = delete;
	TopologyParser(TopologyParser &&other) = delete;

	int N() {
		return _N;
	}

	int N_strands() {
		return _N_strands;
	}

	void parse(std::vector<BaseParticle *> &particles);

protected:
	int _N = -1;
	int _N_strands = -1;
	std::string _filename;

	std::vector<int> _btypes_from_sequence(const std::string &sequence) const;

	// these two methods return the number of initialised particles
	int _parse_old_topology(std::ifstream &topology, std::vector<BaseParticle *> &particles);
	int _parse_new_topology(std::ifstream &topology, std::vector<BaseParticle *> &particles);
};

#endif /* SRC_UTILITIES_TOPOLOGYPARSER_H_ */
