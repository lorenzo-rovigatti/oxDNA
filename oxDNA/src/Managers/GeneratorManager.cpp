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

GeneratorManager::GeneratorManager(int argc, char *argv[]) {
	_use_density = false;
	_density = -1.;
	_box_side = -1.;
	_box_side_x = _box_side_y = _box_side_z = -1.;

	// first we look for an x in argv[2]
	std::string tmpstr = argv[2];
	size_t found_x = tmpstr.find('x');
	if (found_x != std::string::npos) {
		std::vector< std::string > sides_str = Utils::split(tmpstr, 'x');
		_box_side_x = atof(sides_str[0].c_str());
		_box_side_y = atof(sides_str[1].c_str());
		_box_side_z = atof(sides_str[2].c_str());

		_box_side = _box_side_x;
		if (_box_side_y < _box_side) _box_side = _box_side_y;
		if (_box_side_z < _box_side) _box_side = _box_side_z;
		
		OX_LOG (Logger::LOG_INFO, "Detected non cubic box; box sides %g, %g, %g", _box_side_x, _box_side_y, _box_side_z);
		_use_density = false;
	}
	else {
		double argv2 = atof (argv[2]);
		if (argv2 <= 0.) throw oxDNAException ("Refusing to run with box_side/density = %g (converted from \"%s\")", argv2, argv[2]);
		if (argv2 < 2.) {
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

	loadInputFile(&_input, argv[1]);
	if(_input.state == ERROR) throw oxDNAException("Caught an error while opening the input file");
	argc -= 3;
	if(argc > 0) addCommandLineArguments(&_input, argc, argv+3);

	_N = 0;
	_particles = NULL;
	_interaction = NULL;
	_init_completed = false;

	_external_forces = false;
	_external_filename = std::string("");
}

GeneratorManager::~GeneratorManager() {
	cleanInputFile(&_input);
	if(_init_completed) {
		delete _interaction;
		for(int i = 0; i < _N; i++) delete _particles[i];
		delete[] _particles;
	}

	if (_external_forces) ForceFactory<double>::instance()->clear();
}

void GeneratorManager::load_options() {
	getInputString(&_input, "trajectory_file", _trajectory, 1);
	getInputString(&_input, "conf_file", _output_conf, 1);
	
	// read wether to use external forces
	getInputBool(&_input, "external_forces", &_external_forces, 0);
	if (_external_forces) getInputString(&_input, "external_forces_file", _external_filename, 0);

	// seed;
	int seed;
	if (getInputInt(&_input, "seed", &seed, 0) == KEY_NOT_FOUND) {
		seed = time(NULL);
		int rand_seed = 0;
		FILE *f = fopen("/dev/urandom", "rb");
		if (f == NULL) {
			OX_LOG(Logger::LOG_INFO, "Can't open /dev/urandom, using system time as a random seed");
		}
		else {
			if(fread((void *) &rand_seed, sizeof(rand_seed), 1, f) != 0) seed += rand_seed;
			else OX_LOG(Logger::LOG_INFO, "Can't read from /dev/urandom, using system time as a random seed");
			fclose(f);
		}
	}
	OX_LOG(Logger::LOG_INFO, "Setting the random number generator with seed = %d", seed);
	srand48((long int)seed);

	Logger::instance()->get_settings(_input);

	PluginManager *pm = PluginManager::instance();
	pm->init(_input);

	// added by FR
	_mybox = BoxFactory::make_box<double>(_input);

	_interaction = InteractionFactory::make_interaction<double>(_input);
	_interaction->get_settings(_input);
}

void GeneratorManager::init() {
	_interaction->init();

	_N = _interaction->get_N_from_topology();
	_particles = new BaseParticle<double>*[_N];
	int N_strands;
	_interaction->read_topology(_N, &N_strands, _particles);

	if (_use_density){
		 _box_side = pow(_N / _density, 1. / 3.);
		 _box_side_x = _box_side_y = _box_side_z = _box_side;
		OX_LOG(Logger::LOG_INFO, "Generating cubic configuration with box_side %g", _box_side);
	}
	else {
		_density = _N / (_box_side_x * _box_side_y * _box_side_z);
		OX_LOG(Logger::LOG_INFO, "Generating configuration with density %g (%d particles, box sides %g %g %g)", _density, _N, _box_side_x, _box_side_y, _box_side_z);
	}
	
	// setting the box side for the interaction
	_interaction->set_box_side(_box_side);

	_mybox->init(_box_side_x, _box_side_y, _box_side_z);
	_interaction->set_box(_mybox);

	// initializing external forces
	if (_external_forces) { 
		ForceFactory<double>::instance()->read_external_forces(_external_filename, _particles, this->_N, false, &_box_side);
	}
 
	_init_completed = true;
}

void GeneratorManager::generate() {
	_interaction->generate_random_configuration(_particles, _N, _box_side);

	ofstream conf_output(_output_conf);
	conf_output.precision(15);

	conf_output << "t = 0" << endl;
	conf_output << "b = " << _box_side_x << " " << _box_side_y << " " << _box_side_z << endl;
	conf_output << "E = 0 0 0 " << endl;

	for(int i = 0; i < _N; i++) {
		BaseParticle<double> *p = _particles[i];
		LR_vector<double> mypos = p->get_abs_pos(_box_side_x, _box_side_y, _box_side_z);
		LR_matrix<double> oT = p->orientation.get_transpose();
		conf_output << mypos.x << " " << mypos.y << " " << mypos.z << " ";
		conf_output << oT.v1.x << " " << oT.v1.y << " " << oT.v1.z << " ";
		conf_output << oT.v3.x << " " << oT.v3.y << " " << oT.v3.z << " ";
		conf_output << p->vel.x << " " << p->vel.y << " " << p->vel.z << " ";
		conf_output << p->L.x << " " << p->L.y << " " << p->L.z << endl;
	}
	conf_output.close();
}
