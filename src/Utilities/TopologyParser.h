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

	bool is_new_topology();

	int parse_old_topology(std::vector<BaseParticle *> &particles);
	void parse_new_topology();

	std::vector<input_file> lines() {
		return _lines;
	}

protected:
	int _N = -1;
	int _N_strands = -1;
	std::string _filename;
	bool _is_new_topology;
	std::vector<input_file> _lines;
};

#endif /* SRC_UTILITIES_TOPOLOGYPARSER_H_ */
