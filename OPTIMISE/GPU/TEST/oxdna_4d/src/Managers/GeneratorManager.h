/*
 * GeneratorManager.h
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#ifndef GENERATORMANAGER_H_
#define GENERATORMANAGER_H_

#include "../defs.h"
#include "../Backends/AnalysisBackend.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

/**
 * @brief Manages the generation of an initial configuration.
 */
class GeneratorManager {
protected:
	input_file _input;
	char _output_conf[256];
	char _trajectory[256];

	InteractionPtr _interaction;

	bool _use_density;
	double _box_side;
	double _box_side_x, _box_side_y, _box_side_z;
	double _density;
	int _N;
	std::vector<BaseParticle *> _particles;
	std::vector<std::shared_ptr<Molecule>> _molecules;

	bool _external_forces;
	std::string _external_filename;

	std::shared_ptr<BaseBox> _mybox;

public:
	GeneratorManager(input_file input, char *third_argument);
	virtual ~GeneratorManager();

	void load_options();
	void init();
	void generate();
};

#endif /* ANALYSISMANAGER_H_ */
