/*
 * GeneratorManager.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#include "GeneratorManager.h"
#include "../Interactions/InteractionFactory.h"
#include "../Forces/ForceFactory.h"
#include "../PluginManagement/PluginManager.h"
#include "../Boxes/BoxFactory.h"

GeneratorManager::GeneratorManager(input_file input, char *third_argument) :
				_input(input) {
	_use_density = false;
	_density = -1.;
	_box_side = -1.;
	_box_side_x = _box_side_y = _box_side_z = -1.;

	// first we look for an x in argv[2]
	std::string tmpstr(third_argument);
	size_t found_x = tmpstr.find('x');
	if(found_x != std::string::npos) {
		std::vector<std::string> sides_str = Utils::split(tmpstr, 'x');
		_box_side_x = atof(sides_str[0].c_str());
		_box_side_y = atof(sides_str[1].c_str());
		_box_side_z = atof(sides_str[2].c_str());

		_box_side = _box_side_x;
		if(_box_side_y < _box_side) _box_side = _box_side_y;
		if(_box_side_z < _box_side) _box_side = _box_side_z;

		OX_LOG(Logger::LOG_INFO, "Detected non cubic box; box sides %g, %g, %g", _box_side_x, _box_side_y, _box_side_z);
		_use_density = false;
	}
	else {
		number argv2 = atof(tmpstr.c_str());
		if(argv2 <= 0.) {
			throw oxDNAException("Refusing to run with box_side/density = %g (converted from \"%s\")", argv2, tmpstr.c_str());
		}
		if(argv2 < 2.) {
			_use_density = true;
			_density = argv2;
			OX_LOG(Logger::LOG_INFO, "Generating configuration with density %g", _density);
		}
		else {
			_use_density = false;
			_box_side = argv2;
			_box_side_x = _box_side_y = _box_side_z = _box_side;
			OX_LOG(Logger::LOG_INFO, "Generating configuration with side %g", _box_side);
		}
	}

	_N = 0;
	_interaction = NULL;

	_external_forces = false;
	_external_filename = std::string("");

	ConfigInfo::init(&_particles, &_molecules);
}

GeneratorManager::~GeneratorManager() {
	for(auto particle : _particles) {
		delete particle;
	}
}

void GeneratorManager::load_options() {
	CONFIG_INFO->sim_input = &_input;

	getInputString(&_input, "trajectory_file", _trajectory, 1);
	getInputString(&_input, "conf_file", _output_conf, 1);

	// seed;
	int seed;
	if(getInputInt(&_input, "seed", &seed, 0) == KEY_NOT_FOUND) {
		seed = time(NULL);
		int rand_seed = 0;
		FILE *f = fopen("/dev/urandom", "rb");
		if(f == NULL) {
			OX_LOG(Logger::LOG_INFO, "Can't open /dev/urandom, using system time as a random seed");
		}
		else {
			if(fread((void *) &rand_seed, sizeof(rand_seed), 1, f) != 0) seed += rand_seed;
			else OX_LOG(Logger::LOG_INFO, "Can't read from /dev/urandom, using system time as a random seed");
			fclose(f);
		}
	}
	OX_LOG(Logger::LOG_INFO, "Setting the random number generator with seed = %d", seed);
	srand48((long int) seed);

	Logger::instance()->get_settings(_input);

	PluginManager::instance()->init(_input);

	_mybox = std::shared_ptr<BaseBox>(BoxFactory::make_box(_input));

	_interaction = InteractionFactory::make_interaction(_input);
	_interaction->get_settings(_input);

	getInputBool(&_input, "external_forces", &_external_forces, 0);
}

void GeneratorManager::init() {
	_interaction->init();

	_N = _interaction->get_N_from_topology();
	_particles.resize(_N);
	int N_strands;
	_interaction->read_topology(&N_strands, _particles);

	// initialise the molecules
	_molecules.resize(N_strands, std::make_shared<Molecule>());
	for(auto p : _particles) {
		int mol_id = p->strand_id;
		_molecules[mol_id]->add_particle(p);
	}

	if(_use_density) {
		_box_side = pow(_N / _density, 1. / 3.);
		_box_side_x = _box_side_y = _box_side_z = _box_side;
		OX_LOG(Logger::LOG_INFO, "Generating cubic configuration with box_side %g", _box_side);
	}
	else {
		_density = _N / (_box_side_x * _box_side_y * _box_side_z);
		OX_LOG(Logger::LOG_INFO, "Generating configuration with density %g (%d particles, box sides %g %g %g)", _density, _N, _box_side_x, _box_side_y, _box_side_z);
	}

	// setting the box of the interaction
	_mybox->init(_box_side_x, _box_side_y, _box_side_z);
	_interaction->set_box(_mybox.get());

	CONFIG_INFO->interaction = _interaction.get();
	CONFIG_INFO->box = _mybox.get();

	// initialise external forces
	ForceFactory::instance()->make_forces(_particles, _mybox.get());
}

void GeneratorManager::generate() {
	_interaction->generate_random_configuration(_particles);

	ofstream conf_output(_output_conf);
	conf_output.precision(15);

	conf_output << "t = 0" << endl;
	conf_output << "b = " << _mybox->box_sides().x << " " << _mybox->box_sides().y << " " << _mybox->box_sides().z << endl;
	conf_output << "E = 0 0 0 " << endl;

	for(int i = 0; i < _N; i++) {
		BaseParticle *p = _particles[i];
		LR_vector mypos = _mybox->get_abs_pos(p);
		LR_matrix oT = p->orientation.get_transpose();
		conf_output << mypos.x << " " << mypos.y << " " << mypos.z << " ";
		conf_output << oT.v1.x << " " << oT.v1.y << " " << oT.v1.z << " ";
		conf_output << oT.v3.x << " " << oT.v3.y << " " << oT.v3.z << " ";
		conf_output << p->vel.x << " " << p->vel.y << " " << p->vel.z << " ";
		conf_output << p->L.x << " " << p->L.y << " " << p->L.z << endl;
	}
	conf_output.close();
}
