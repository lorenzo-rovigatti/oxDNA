/*
 * SimBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include <sstream>

#include "SimBackend.h"
#include "../Utilities/Utils.h"
#include "../Interactions/InteractionFactory.h"
#include "../Forces/ForceFactory.h"
#include "../Lists/ListFactory.h"
#include "../PluginManagement/PluginManager.h"

template<typename number>
SimBackend<number>::SimBackend() {
	// we need to initialize everything so that we can check what we can
	// and what we can't delete[] in the destructor
	_timer_msgs_number = 0;
	_enable_fix_diffusion = true;
	_print_timings = false;
	_external_forces = false;
	_N_updates = 0;
	_confs_to_skip = 0;
	_particles = NULL;
	_sim_type = -1;
	_is_CUDA_sim = false;
	_interaction = NULL;
	_N_strands = -1;
	_N = 0;
	// here to avoid valgrind's "Conditional jump or move depends on
	// uninitialised value"
	_timer.state = TIMER_NOT_INIT;
	_start_step_from_file = (llint) 0;
	_U = (number) 0.f;
	_K = (number) 0.f;
	_U_hydr = (number) 0.f;
	_backend_info = std::string("");
	_initial_conf_is_binary = false;
	_reseed = false;
	_back_in_box = false;
	_custom_conf_name = false;
	_obs_output_last_conf_bin = NULL;
	_P = (number) 0.;
	_lists = NULL;
}

template<typename number>
SimBackend<number>::~SimBackend() {
	if(_particles != NULL) {
		for(int i = 0; i < _N; i++) delete _particles[i];
		delete[] _particles;
	}
	if(_interaction != NULL) delete _interaction;

	ForceFactory<number>::instance()->clear();

	typename vector<ObservableOutput<number> *>::iterator it;
	for(it = _obs_outputs.begin(); it != _obs_outputs.end(); it++) delete *it;

	if(_timer.state == TIMER_READY) {
		prepare_timer_results(&_timer);

		OX_LOG(Logger::LOG_INFO, "Timings informations:\n");
		print_times(&_timer, Logger::instance()->get_log_stream());
		OX_LOG(Logger::LOG_NOTHING, "");

		if(_print_timings == true) {
			FILE *timings_file = fopen(_timings_filename, "a");
			fprintf(timings_file, "%d %lf\n", _N, _timer.timings[0]);
			fclose(timings_file);
		}

		destroy_timer(&_timer);
	}

	// destroy lists;
	if (_lists != NULL) delete _lists;

	PluginManager::clear();
}

template<typename number>
void SimBackend<number>::get_settings(input_file &inp) {
	int tmp;

	// initialise the plugin manager with the input file
	PluginManager *pm = PluginManager::instance();
	pm->init(inp);

	_interaction = InteractionFactory::make_interaction<number>(inp);
	_interaction->get_settings(inp);

	_lists = ListFactory::make_list<number>(inp, _N, _box_side);
	_lists->get_settings(inp);

	int tmpi;
	if (getInputBoolAsInt(&inp, "binary_initial_conf", &tmpi, 0) == KEY_FOUND) _initial_conf_is_binary = (tmpi > 0);
	if (getInputBoolAsInt(&inp, "fix_diffusion", &tmpi, 0) == KEY_FOUND) _enable_fix_diffusion = (tmpi > 0);

	// we only reseed the RNG if:
	// a) we have a binary conf
	// b) no seed was specified in the input file
	// c) restart_step_counter is set to one (true) -- LR: looking at the code, it looks like it is quite the contrary
	getInputBoolAsInt (&inp, "restart_step_counter", &tmpi, 1);
	_reseed = (getInputInt (&inp, "seed", &tmpi, 0) == KEY_NOT_FOUND || tmpi == 0);

	getInputInt(&inp, "confs_to_skip", &_confs_to_skip, 0);

	int val = getInputBoolAsInt(&inp, "external_forces", &tmp, 0);
	if(val == KEY_FOUND) {
		_external_forces = (tmp != 0);
		if (_external_forces) {
			getInputString (&inp, "external_forces_file", _external_filename, 1);
		}
	}
	else if(val == KEY_INVALID) throw oxDNAException("external_forces must be either 0 (false, no) or 1 (true, yes)");

	if(getInputBoolAsInt(&inp, "print_timings", &tmp, 0) == KEY_FOUND) {
		_print_timings = (tmp != 0);
		if(_print_timings) {
			getInputString(&inp, "timings_filename", _timings_filename, 1);
			FILE *timings_file = fopen(_timings_filename, "a");
			if(timings_file == NULL) throw oxDNAException("Timings file '%s' is not writable\n", _timings_filename);
			fclose(timings_file);
		}
	}

	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	_T = Utils::get_temperature<number>(raw_T);

	// here we fill the _obs_outputs vector
	int i = 1;
	bool found = true;
	while(found) {
		stringstream ss;
		ss << "data_output_" << i;
		string obs_string;
		if(getInputString(&inp, ss.str().c_str(), obs_string, 0) == KEY_FOUND) {
			ObservableOutput<number> *new_obs_out = new ObservableOutput<number>(obs_string, inp);
			_obs_outputs.push_back(new_obs_out);
		}
		else found = false;

		i++;
	}

	if(getInputBoolAsInt(&inp, "back_in_box", &tmp, 0) == KEY_FOUND){
		_back_in_box = (tmp != 0);
		if (_back_in_box) OX_LOG(Logger::LOG_INFO, "ascii configuration files will have the particles put back in the box");
	}

	// just to print the message
	if(getInputBoolAsInt(&inp, "major_minor_grooving", &tmp, 0) == KEY_FOUND){
		if (tmp != 0) OX_LOG(Logger::LOG_INFO, "Using different widths for major and minor grooves");
	}

	// we build the default stream of observables for trajectory and last configuration
	char traj_file[512];

	// Trajectory
	getInputString(&inp, "trajectory_file", traj_file, 1);
	std::string fake = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n}\n", traj_file);
	_obs_output_trajectory = new ObservableOutput<number>(fake, inp);
	_obs_output_trajectory->add_observable("type = configuration");
	_obs_outputs.push_back(_obs_output_trajectory);

	// Last configuration
	char lastconf_file[512] = "last_conf.dat";
	getInputString(&inp, "lastconf_file", lastconf_file, 0);
	fake = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n\tonly_last = 1\n}\n", lastconf_file);
	_obs_output_last_conf = new ObservableOutput<number>(fake, inp);
	_obs_output_last_conf->add_observable("type = configuration");
	_obs_outputs.push_back(_obs_output_last_conf);

	// Last configuration in binary, optional
	char lastconf_file_bin[512];
	if((getInputString(&inp, "lastconf_file_bin", lastconf_file_bin, 0) == KEY_FOUND)) {
		fake = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n\tonly_last = 1\n\tbinary = 1\n}\n", lastconf_file_bin);
		_obs_output_last_conf_bin = new ObservableOutput<number>(fake, inp);
		_obs_output_last_conf_bin->add_observable("type = binary_configuration");
		_obs_outputs.push_back(_obs_output_last_conf_bin);
	}

	// Reduced configurations, optional
	llint reduced_conf_every;
	if(getInputLLInt(&inp, "print_reduced_conf_every", &reduced_conf_every, 0) == KEY_FOUND && reduced_conf_every > 0) {
		getInputString(&inp, "reduced_conf_output_dir", _reduced_conf_output_dir, 1);
		fake = Utils::sformat("{\n\tname = reduced_conf.dat\n\tprint_every = %lld\n\tonly_last = 1\n}\n", reduced_conf_every);
		_obs_output_reduced_conf = new ObservableOutput<number>(fake, inp);
		_obs_output_reduced_conf->add_observable("type = configuration\nreduced = true");
		_obs_outputs.push_back(_obs_output_reduced_conf);
	}
}


template<typename number>
void SimBackend<number>::init(char conf_filename[256]) {
	init_timer(&_timer, _timer_msgs_number, _timer_msgs);

	_conf_input.open(conf_filename);
	if(_conf_input.good() == false) throw oxDNAException("Can't read configuration file '%s'", conf_filename);

	_interaction->init();

	_rcut = _interaction->get_rcut();
	_sqr_rcut = SQR(_rcut);

	// check number of particles
	_N = _interaction->get_N_from_topology();
	_particles = new BaseParticle<number>*[_N];
	_interaction->read_topology(_N, &_N_strands, _particles);

	_conf_input.seekg(0);

	// we need to skip a certain number of lines, depending on how many
	// particles we have and how many configurations we want to skip
	if(_confs_to_skip > 0) {
		OX_LOG(Logger::LOG_INFO, "Skipping %d configuration(s)", _confs_to_skip);
		int i;
		for(i = 0; i < _confs_to_skip && _conf_input.good(); i++) {
			bool check = false;
			if (_initial_conf_is_binary) check = _read_next_configuration(true);
			else check = _read_next_configuration();
			if(!check) throw oxDNAException("Skipping %d configuration(s) is not possible, as the initial configuration file only contains %d.", _confs_to_skip, i);
		}
	}

	bool check = false;
	if (_initial_conf_is_binary) check = _read_next_configuration(true);
	else check = _read_next_configuration();
	if (!check) throw oxDNAException("Could not read the initial configuration, aborting");

	_start_step_from_file = _read_conf_step;

	if (_external_forces) ForceFactory<number>::instance()->read_external_forces(std::string(_external_filename), _particles, this->_N, _is_CUDA_sim, &_box_side);

	this->_U = (number) 0;
	this->_K = (number) 0;
	this->_U_stack = (number) 0;

	_interaction->set_box_side(_box_side);

	_lists->init(_particles, _rcut);

	// initializes the observable output machinery. This part has to follow
	// read_topology() since _particles has to be initialized
	ConfigInfo<number> conf_info(_particles, &_box_side, _interaction, &_N, &_backend_info, _lists);
	typename vector<ObservableOutput<number> *>::iterator it;
	for(it = _obs_outputs.begin(); it != _obs_outputs.end(); it++) (*it)->init(conf_info);

	OX_LOG(Logger::LOG_INFO, "N: %d", _N);
}

// apparently this is the right syntax to have a templated method in a templated classes. Oh c++...
template<typename number>
template<typename number_n>
LR_vector<number_n> SimBackend<number>::_read_next_vector(bool binary) {
	LR_vector<number_n> res;
	if(binary) {
		double tmpf;
		_conf_input.read ((char *)&tmpf, sizeof(double));
		res.x = tmpf;
		_conf_input.read ((char *)&tmpf, sizeof(double));
		res.y = tmpf;
		_conf_input.read ((char *)&tmpf, sizeof(double));
		res.z = tmpf;
	}
	else _conf_input >> res.x >> res.y >> res.z;

	return res;
}

template<typename number>
bool SimBackend<number>::_read_next_configuration(bool binary) {
	double box_side;
	// parse headers. Binary and ascii configurations have different headers, and hence
	// we have to separate the two procedures
	if(binary) {
		OX_LOG(Logger::LOG_INFO, "Reading binary configuration");

		// first bytes: step
		_conf_input.read ((char *)&_read_conf_step, sizeof(llint));

		unsigned short rndseed[3];
		_conf_input.read ((char *)rndseed, 3 * sizeof(unsigned short));
		// we only use the seed if:
		// a) seed was not specified;
		// b) restart_step_counter  == 0;
		if (_reseed) {
			OX_LOG(Logger::LOG_INFO,"Overriding seed: restoring seed: %hu %hu %hu from binary conf", rndseed[0], rndseed[1], rndseed[2]);
			seed48 (rndseed);
		}

		double tmpf;
		_conf_input.read ((char *)&box_side, sizeof(double));
		_conf_input.read ((char *)&tmpf, sizeof(double));
		_conf_input.read ((char *)&tmpf, sizeof(double));

		// energies and such; discarded
		_conf_input.read ((char *)&tmpf, sizeof(double));
		_conf_input.read ((char *)&tmpf, sizeof(double));
		_conf_input.read ((char *)&tmpf, sizeof(double));
	}
	else {
		char line[512];
		_conf_input.getline(line, 512);
		if(_conf_input.eof()) return false;
		int res = sscanf(line, "t = %lld", &_read_conf_step);
		_conf_input.getline(line, 512);

		res += sscanf(line, "b = %lf %*f %*f", &box_side);
		if(res != 2) throw oxDNAException("Malformed headers found in the initial configuration file");

		_conf_input.getline(line, 512);
	}

	_box_side = (number) box_side;

	// the following part is in double precision since
	// we want to be able to restart from confs that have
	// large numbers in the conf file and use float precision
	// later
	int *nins, k, i;
	LR_vector<double> *tmp_poss, *scdm;

	tmp_poss = new LR_vector<double>[_N];
	nins = new int[_N_strands];
	scdm = new LR_vector<double> [_N_strands];
	for (k = 0; k < _N_strands; k ++) {
		nins[k] = 0;
		scdm[k] = LR_vector<double> ((double)0., (double)0., (double)0.);
	}

	i = 0;
	while(!_conf_input.eof() && i < _N) {
		BaseParticle<number> *p = this->_particles[i];

		tmp_poss[i] = _read_next_vector<double>(binary);
		k = p->strand_id;
		scdm[k] += tmp_poss[i];
		nins[k] ++;

		if (!binary) {
			p->orientation.v1 = _read_next_vector<number>(binary);
			p->orientation.v3 = _read_next_vector<number>(binary);
			// get v2 from v1 and v3
			p->orientation.v1.normalize();
			p->orientation.v3.normalize();
			p->orientation.v1 -= p->orientation.v3 * (p->orientation.v1*p->orientation.v3);
			p->orientation.v1.normalize();
			p->orientation.v2 = p->orientation.v3.cross(p->orientation.v1);
		}
		else {
			int x, y, z;
			_conf_input.read ((char *)&x, sizeof(int));
			_conf_input.read ((char *)&y, sizeof(int));
			_conf_input.read ((char *)&z, sizeof(int));
			p->set_pos_shift (x, y, z);

			p->orientation.v1 = _read_next_vector<number>(binary);
			p->orientation.v2 = _read_next_vector<number>(binary);
			p->orientation.v3 = _read_next_vector<number>(binary);
		}
		p->orientation.transpone();

		p->vel = _read_next_vector<number>(binary);
		p->L = _read_next_vector<number>(binary);

		p->init();
		p->set_positions();
		p->orientationT = p->orientation.get_transpose();

		i++;
	}

	// this is needed because, if reading from an ascii trajectory, at this stage the _conf_input pointer points to a \n
	if(!binary && !_conf_input.eof()) {
		char line[512];
		_conf_input.getline(line, 512);
	}
	
	// discarding the final '\n' in the binary file...
	if (binary && !_conf_input.eof()) {
		char tmpc;
		_conf_input.read ((char *)&tmpc, sizeof(char));
	}

	if(i != this->_N) {
		if(_confs_to_skip > 0) throw oxDNAException("Wrong number of particles (%d) found in configuration. Maybe you skipped too many configurations?", i);
		else throw oxDNAException("The number of lines found in configuration file (%d) doesn't match the parsed number of particles (%d)", i, this->_N);
	}

	for(k = 0; k < this->_N_strands; k ++) scdm[k] /= (double) nins[k];

	for (i = 0; i < this->_N; i ++) {
		BaseParticle<number> *p = this->_particles[i];
		k = p->strand_id;

		LR_vector<double> p_pos = tmp_poss[i];
		if(_enable_fix_diffusion && !binary) {
			// we need to manually set the particle shift so that the particle absolute position
			// is the right one
			LR_vector<number> scdm_number(scdm[k].x, scdm[k].y, scdm[k].z);
			p->shift(scdm_number, box_side);
			p_pos.x -= box_side * floor(scdm[k].x / box_side + 0.001);
			p_pos.y -= box_side * floor(scdm[k].y / box_side + 0.001);
			p_pos.z -= box_side * floor(scdm[k].z / box_side + 0.001);
		}
		p->pos = LR_vector<number>(p_pos.x, p_pos.y, p_pos.z);
	}

	_interaction->check_input_sanity(_particles, _N);

	delete[] tmp_poss;
	delete[] scdm;
	delete[] nins;

	return true;
}

template<typename number>
void SimBackend<number>::_print_ready_observables(llint curr_step) {
	typename vector<ObservableOutput<number> *>::iterator it;
	for(it = _obs_outputs.begin(); it != _obs_outputs.end(); it++) {
		if((*it)->is_ready(curr_step)) (*it)->print_output(curr_step);
	}
}

template<typename number>
void SimBackend<number>::print_observables(llint curr_step) {
	bool someone_ready = false;
	typename vector<ObservableOutput<number> *>::iterator it;
	for(it = _obs_outputs.begin(); it != _obs_outputs.end(); it++) {
		if((*it)->is_ready(curr_step)) someone_ready = true;
	}
	if(someone_ready) _print_ready_observables(curr_step);
	_backend_info = std::string ("");
}

template<typename number>
void SimBackend<number>::fix_diffusion() {
	if(!_enable_fix_diffusion) return;

	LR_vector<number> *scdm = new LR_vector<number>[this->_N_strands];
	int *ninstrand = new int [this->_N_strands];
	for(int k = 0; k < this->_N_strands; k ++) {
		scdm[k].x = scdm[k].y = scdm[k].z = (number) 0.f;
		ninstrand[k] = 0;
	}

	// compute com for each strand;
	for (int i = 0; i < this->_N; i ++) {
		BaseParticle<number> *p = this->_particles[i];
		scdm[p->strand_id] += p->pos;
		ninstrand[p->strand_id]++;
	}

	for (int k = 0; k < this->_N_strands; k ++) scdm[k] = scdm[k] / (number) ninstrand[k];

	// change particle position and fix orientation matrix;
	for (int i = 0; i < this->_N; i ++) {
		BaseParticle<number> *p = this->_particles[i];
		p->shift(scdm[p->strand_id], this->_box_side);
		p->orientation.orthonormalize();
		p->orientationT = p->orientation.get_transpose();
		p->set_positions();
	}

	delete[] ninstrand;
	delete[] scdm;
}

template<typename number>
void SimBackend<number>::print_conf(llint curr_step, bool reduced, bool only_last) {
	if(reduced) {
		char conf_name[512];
		sprintf(conf_name, "%s/reduced_conf%lld.dat", _reduced_conf_output_dir, curr_step);
		_obs_output_reduced_conf->change_output_file(conf_name);
		_obs_output_reduced_conf->print_output(curr_step);
	}
	else {
		if(!only_last) _obs_output_trajectory->print_output(curr_step);
		_obs_output_last_conf->print_output(curr_step);
		if(_obs_output_last_conf_bin != NULL) _obs_output_last_conf_bin->print_output(curr_step);
	}
}

template class SimBackend<float>;
template class SimBackend<double>;

