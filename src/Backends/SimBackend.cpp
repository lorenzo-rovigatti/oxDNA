/*
 * SimBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include <sstream>
#include <fstream>

#include "SimBackend.h"
#include "../Utilities/Utils.h"
#include "../Utilities/ConfigInfo.h"
#include "../Interactions/InteractionFactory.h"
#include "../Observables/ObservableFactory.h"
#include "../Forces/ForceFactory.h"
#include "../Lists/ListFactory.h"
#include "../Boxes/BoxFactory.h"
#include "../PluginManagement/PluginManager.h"
#include "../Particles/BaseParticle.h"
#include "../Utilities/Timings.h"

SimBackend::SimBackend() {
	// we need to initialize everything so that we can check what we can
	// and what we can't delete[] in the destructor
	_enable_fix_diffusion = true;
	_external_forces = false;
	_N_updates = 0;
	_confs_to_skip = 0;
	_bytes_to_skip = 0;
	_interaction = nullptr;
	_N_strands = -1;
	start_step_from_file = (llint) 0;
	_U = (number) 0.f;
	_backend_info = std::string("");
	_initial_conf_is_binary = false;
	_reseed = false;
	_back_in_box = false;
	_custom_conf_name = false;
	_read_conf_step = 0;
	_obs_output_last_conf_bin = nullptr;
	_P = (number) 0.;
	_lists = nullptr;
	_box = nullptr;
	_max_io = 1.e30;
	_read_conf_step = -1;
	_conf_interval = -1;
	_obs_output_trajectory = _obs_output_stdout = _obs_output_file = _obs_output_reduced_conf = _obs_output_last_conf = _obs_output_checkpoints = _obs_output_last_checkpoint = nullptr;
	_mytimer = nullptr;
	_restart_step_counter = false;
	_rcut = -1.;
	_sqr_rcut = -1.;
	_T = -1.;

	ConfigInfo::init(&_particles, &_molecules);
	_config_info = ConfigInfo::instance().get();
	_config_info->subscribe("T_updated", [this]() {
		this->_on_T_update();
	});
}

SimBackend::~SimBackend() {
	for(auto particle : _particles) {
		delete particle;
	}

	// here we print the input output information
	llint total_file = 0;
	llint total_stderr = 0;
	OX_LOG(Logger::LOG_INFO, "Aggregated I/O statistics (set debug=1 for file-wise information)");
	for(auto const &element : _obs_outputs) {
		llint now = element->get_bytes_written();
		auto fname = element->get_output_name();
		if(!strcmp(fname.c_str(), "stderr") || !strcmp(fname.c_str(), "stdout")) {
			total_stderr += now;
		}
		else {
			total_file += now;
		}
		std::string mybytes = Utils::bytes_to_human(now);
		OX_DEBUG("  on %s: %s", fname.c_str(), mybytes.c_str());
	}
	OX_LOG(Logger::LOG_NOTHING, "\t%s written to files" , Utils::bytes_to_human(total_file).c_str());
	OX_LOG(Logger::LOG_NOTHING, "\t%s written to stdout/stderr" , Utils::bytes_to_human(total_stderr).c_str());

	if(_mytimer != nullptr) {
		double time_passed = (double) _mytimer->get_time() / (double) CLOCKS_PER_SEC;
		if(time_passed > 0.) {
			OX_LOG(Logger::LOG_NOTHING, "\tFor a total of %8.3lg MB/s\n", (total_file + total_stderr) / ((1024.*1024.) * time_passed));
		}
	}

	// we explicitly clean up ConfigInfo so that observables are destructed before PluginManager unloads
	// any dynamically-linked library
	ConfigInfo::clear();
}

void SimBackend::get_settings(input_file &inp) {
	int tmp;

	_config_info->sim_input = &inp;

	// initialise the plugin manager with the input file
	PluginManager::instance()->init(inp);

	// initialise the timer
	_mytimer = TimingManager::instance()->new_timer(std::string("SimBackend"));
	_obs_timer = TimingManager::instance()->new_timer(string("Observables"), string("SimBackend"));

	_interaction = InteractionFactory::make_interaction(inp);
	_interaction->get_settings(inp);

	_box = BoxFactory::make_box(inp);
	_box->get_settings(inp);

	_lists = ListFactory::make_list(inp, _particles, _box.get());
	_lists->get_settings(inp);

	getInputBool(&inp, "restart_step_counter", &_restart_step_counter, 0);

	// reload configuration
	std::string reload_from;
	if(getInputString(&inp, "reload_from", reload_from, 0) == KEY_FOUND) {
		// set default variables in this case
		_initial_conf_is_binary = true;

		// check that conf_file is not specified
		std::string tmpstring;
		if(getInputString(&inp, "conf_file", tmpstring, 0) == KEY_FOUND)
			throw oxDNAException("Input file error: \"conf_file\" cannot be specified if \"reload_from\" is specified");

		// check that restart_step_counter is set to 0
		if(_restart_step_counter)
			throw oxDNAException("Input file error: \"restart_step_counter\" must be set to false if \"reload_from\" is specified");

		int my_seed;
		if(getInputInt(&inp, "seed", &my_seed, 0) == KEY_FOUND)
			throw oxDNAException("Input file error: \"seed\" must not be specified if \"reload_from\" is specified");

		_conf_filename = std::string(reload_from);
	}
	else {
		getInputString(&inp, "conf_file", _conf_filename, 1);
	}

	getInputBool(&inp, "binary_initial_conf", &_initial_conf_is_binary, 0);
	if(_initial_conf_is_binary) {
		OX_LOG(Logger::LOG_INFO, "Reading binary configuration");
	}

	getInputBool(&inp, "fix_diffusion", &_enable_fix_diffusion, 0);

	// we only reseed the RNG if:
	// a) we have a binary conf
	// b) no seed was specified in the input file
	// c) restart_step_counter is set to one (true) -- LR: looking at the code, it looks like it is quite the contrary
	int tmpi;
	getInputBoolAsInt(&inp, "restart_step_counter", &tmpi, 1);
	_reseed = (getInputInt(&inp, "seed", &tmpi, 0) == KEY_NOT_FOUND || tmpi == 0);

	getInputInt(&inp, "confs_to_skip", &_confs_to_skip, 0);
	if(_confs_to_skip == 0) {
		getInputLLInt(&inp, "bytes_to_skip", &_bytes_to_skip, 0);
	}

	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	_config_info->update_temperature(Utils::get_temperature(raw_T));

	// here we fill the _obs_outputs vector
	auto new_outputs = ObservableFactory::make_observables();
	for(auto new_output : new_outputs) {
		add_output(new_output);
	}

	getInputBool(&inp, "back_in_box", &_back_in_box, 0);
	if(_back_in_box) {
		OX_LOG(Logger::LOG_INFO, "ascii configuration files will have the particles put back in the box");
	}

	// just to print the message
	if(getInputBoolAsInt(&inp, "major_minor_grooving", &tmp, 0) == KEY_FOUND) {
		if(tmp != 0) {
			OX_LOG(Logger::LOG_INFO, "Using different widths for major and minor grooves");
		}
	}

	// we build the default stream of observables for trajectory and last configuration
	bool traj_print_momenta = true;
	getInputBool(&inp, "trajectory_print_momenta", &traj_print_momenta, 0);
	std::string traj_file;
	// Trajectory
	getInputString(&inp, "trajectory_file", traj_file, 1);
	std::string output_inp_text = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n}\n", traj_file.c_str());
	_obs_output_trajectory = std::make_shared<ObservableOutput>(output_inp_text);

	std::string obs_text = Utils::sformat("type = configuration\nprint_momenta = %d", traj_print_momenta);
	_obs_output_trajectory->add_observable(obs_text);
	add_output(_obs_output_trajectory);

	// Last configuration
	std::string lastconf_file = "last_conf.dat";
	getInputString(&inp, "lastconf_file", lastconf_file, 0);
	output_inp_text = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n\tonly_last = 1\n}\n", lastconf_file.c_str());
	_obs_output_last_conf = std::make_shared<ObservableOutput>(output_inp_text);
	_obs_output_last_conf->add_observable("type = configuration\nid = last_conf");
	add_output(_obs_output_last_conf);

	// Last configuration in binary, optional
	std::string lastconf_file_bin;
	if((getInputString(&inp, "lastconf_file_bin", lastconf_file_bin, 0) == KEY_FOUND)) {
		output_inp_text = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n\tonly_last = 1\n\tbinary = 1\n}\n", lastconf_file_bin.c_str());
		_obs_output_last_conf_bin = std::make_shared<ObservableOutput>(output_inp_text);
		_obs_output_last_conf_bin->add_observable("type = binary_configuration");
		add_output(_obs_output_last_conf_bin);
	}

	// Reduced configurations, optional
	llint reduced_conf_every;
	if(getInputLLInt(&inp, "print_reduced_conf_every", &reduced_conf_every, 0) == KEY_FOUND && reduced_conf_every > 0) {
		getInputString(&inp, "reduced_conf_output_dir", _reduced_conf_output_dir, 1);
		output_inp_text = Utils::sformat("{\n\tname = reduced_conf.dat\n\tprint_every = %lld\n\tonly_last = 1\n}\n", reduced_conf_every);
		_obs_output_reduced_conf = std::make_shared<ObservableOutput>(output_inp_text);
		_obs_output_reduced_conf->add_observable("type = configuration\nreduced = true");
		add_output(_obs_output_reduced_conf);
	}

	// checkpoints trajectory, optional
	llint checkpoint_every;
	if(getInputLLInt(&inp, "checkpoint_every", &checkpoint_every, 0) == KEY_FOUND && checkpoint_every > 0) {
		int tmp1 = getInputString(&inp, "checkpoint_trajectory", _checkpoint_traj, 0);
		if(tmp1 == KEY_FOUND) {
			output_inp_text = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n\tonly_last = false\n}\n", _checkpoint_traj.c_str(), checkpoint_every);
			_obs_output_checkpoints = std::make_shared<ObservableOutput>(output_inp_text);
			_obs_output_checkpoints->add_observable("type = checkpoint");
			add_output(_obs_output_checkpoints);
			OX_LOG(Logger::LOG_INFO, "Setting up a trajectory of checkpoints to file %s every %lld steps",_checkpoint_traj.c_str(), checkpoint_every);
		}

		int tmp2 = getInputString(&inp, "checkpoint_file", _checkpoint_file, 0);
		if(tmp2 == KEY_FOUND) {
			output_inp_text = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n\tonly_last = true\n}\n", _checkpoint_file.c_str(), checkpoint_every);
			_obs_output_last_checkpoint = std::make_shared<ObservableOutput>(output_inp_text);
			_obs_output_last_checkpoint->add_observable("type = checkpoint");
			add_output(_obs_output_last_checkpoint);
			OX_LOG(Logger::LOG_INFO, "Setting up last checkpoint to file %s every %lld steps",_checkpoint_file.c_str(), checkpoint_every);
		}

		if(tmp1 != KEY_FOUND && tmp2 != KEY_FOUND)
			throw oxDNAException(
					"Input file error: At least one of \"checkpoint_file\" or \"checkpoint_trajectory\" must be specified if \"checkpoint_every\" is specified.");
	}

	// set the max IO
	if(getInputNumber(&inp, "max_io", &_max_io, 0) == KEY_FOUND) {
		if(_max_io < 0)
			throw oxDNAException("Cannot run with a negative I/O limit. Set the max_io key to something > 0");
		else
			OX_LOG(Logger::LOG_INFO, "Setting the maximum IO limit to %g MB/s", _max_io);
	}
	else {
		_max_io = 1.; // default value for a simulation is 1 MB/s;
	}

	getInputBool(&inp, "external_forces", &_external_forces, 0);

}

void SimBackend::init() {
	_conf_input.open(_conf_filename.c_str());
	if(_conf_input.good() == false) {
		throw oxDNAException("Can't read configuration file '%s'", _conf_filename.c_str());
	}

	_interaction->init();

	// check number of particles
	int N = _interaction->get_N_from_topology();
	_particles.resize(N);
	_interaction->read_topology(&_N_strands, _particles);

	_rcut = _interaction->get_rcut();
	_sqr_rcut = SQR(_rcut);

	// check that the interaction has filled the array of "affected" particles
	for(int i = 0; i < N; i++) {
		BaseParticle *p = _particles[i];
		if(p->n3 != P_VIRTUAL || p->n5 != P_VIRTUAL) {
			if(p->affected.size() < 1)
				throw oxDNAException("Found an interaction with bonded interactions that did not set the affected attribute for particle %d. Aborting\n", p->index);
		}
	}

	_conf_input.seekg(0, ios::beg);

	// we need to skip a certain number of lines, depending on how many
	// particles we have and how many configurations we want to skip
	if(_confs_to_skip > 0) {
		OX_LOG(Logger::LOG_INFO, "Skipping %d configuration(s)", _confs_to_skip);
		int i;
		for(i = 0; i < _confs_to_skip && _conf_input.good(); i++) {
			if(!read_next_configuration(_initial_conf_is_binary)) {
				throw oxDNAException("Skipping %d configuration(s) is not possible, as the initial trajectory file only contains %d configurations", _confs_to_skip, i);
			}
		}
	}
	else if(_bytes_to_skip > 0) {
		_conf_input.seekg(_bytes_to_skip, std::ios_base::beg);
	}

	bool check = read_next_configuration(_initial_conf_is_binary);
	if(!check) {
		throw oxDNAException("Could not read the initial configuration, aborting");
	}

	// initialise the molecules
	Molecule::reset_id();
	for(int i = 0; i < _N_strands; i++) {
		_molecules.push_back(std::make_shared<Molecule>());
	}

	for(auto p : _particles) {
		int mol_id = p->strand_id;
		_molecules[mol_id]->add_particle(p);
	}

	if((int) _molecules.size() != _N_strands) {
		throw oxDNAException("There is a mismatch between the number of strands (molecules) reported in the topology file (%d) and the number that results from parsing the particle information (%u)", _N_strands, _molecules.size());
	}

	start_step_from_file = (_restart_step_counter) ? 0 : _read_conf_step;
	_config_info->curr_step = start_step_from_file;

	// initialise external forces
	ForceFactory::instance()->make_forces(_particles, _box.get());

	// now that molecules and external forces have been initialised we can check what strands can be shifted
	// and undo the shifting done to particles that belong to those strands that shouldn't be shifted
	std::vector<bool> printed(_molecules.size(), false); // used to print the message once per strand
	for(auto p : _particles) {
		if(!_molecules[p->strand_id]->shiftable()) {
			p->pos = _box->get_abs_pos(p);
			p->set_pos_shift(0, 0, 0);
			if(!printed[p->strand_id] && _enable_fix_diffusion) {
				OX_LOG(Logger::LOG_INFO, "Strand %d is not shiftable and therefore 'fix_diffusion' will not act on it", p->strand_id);
				printed[p->strand_id] = true;
			}
		}
	}

	_interaction->set_box(_box.get());

	_lists->init(_rcut);
	CONFIG_INFO->subscribe(_box->INIT_EVENT, [this]() { this->_lists->change_box(); });

	_config_info->set(_interaction.get(), &_backend_info, _lists.get(), _box.get());

	// initializes the observable output machinery. This part has to follow read_topology() since _particles has to be initialized
	for(auto const &element : _obs_outputs) {
		element->init();
	}

	OX_LOG(Logger::LOG_INFO, "N: %d, N molecules: %d", N, _molecules.size());
}

LR_vector SimBackend::_read_next_binary_vector() {
	LR_vector res;
	double tmpf;
	_conf_input.read((char*) &tmpf, sizeof(double));
	res.x = tmpf;
	_conf_input.read((char*) &tmpf, sizeof(double));
	res.y = tmpf;
	_conf_input.read((char*) &tmpf, sizeof(double));
	res.z = tmpf;

	return res;
}

// here we cannot use _molecules because it has not been initialised yet
bool SimBackend::read_next_configuration(bool binary) {
	double Lx, Ly, Lz;
	// parse headers. Binary and ascii configurations have different headers, and hence
	// we have to separate the two procedures
	if(binary) {

		// first bytes: step
		// if there's nothing to read, _read_conf_step is unchanged and equal to -1
		_read_conf_step = -1;
		_conf_input.read((char*) &_read_conf_step, sizeof(llint));
		if(_read_conf_step == -1) {
			OX_DEBUG("End of binary configuration file reached.");
			return false;
		}

		unsigned short rndseed[3];
		_conf_input.read((char*) rndseed, 3 * sizeof(unsigned short));
		// we only use the seed if:
		// a) seed was not specified;
		// b) restart_step_counter  == 0;
		if(_reseed) {
			OX_LOG(Logger::LOG_INFO,"Overriding seed: restoring seed: %hu %hu %hu from binary conf", rndseed[0], rndseed[1], rndseed[2]);
			seed48(rndseed);
		}

		double tmpf;
		_conf_input.read((char*) &Lx, sizeof(double));
		_conf_input.read((char*) &Ly, sizeof(double));
		_conf_input.read((char*) &Lz, sizeof(double));

		// energies and such; discarded
		_conf_input.read((char*) &tmpf, sizeof(double));
		_conf_input.read((char*) &tmpf, sizeof(double));
		_conf_input.read((char*) &tmpf, sizeof(double));
	}
	else {
		bool malformed_headers = false;
		std::string line;
		std::getline(_conf_input, line);
		if(_conf_input.eof()) {
			return false;
		}
		int res = sscanf(line.c_str(), "t = %lld", &_read_conf_step);
		std::stringstream error_message;
		// handle the case when t can't be read
		if(res != 1) {
			error_message << "Malformed headers found in an input configuration.\"t = <int>\" was expected, but \"" << line << "\" was found instead.";
			malformed_headers = true;

		}
		if(!malformed_headers) {
			std::getline(_conf_input, line);
			// handle the case when the box_size can't be read
			res += sscanf(line.c_str(), "b = %lf %lf %lf", &Lx, &Ly, &Lz);
			if(res != 4) {
				error_message << "Malformed headers found in an input configuration.\"b = <float> <float> <float>\" was expected, but \"" << line << "\" was found instead.";
				malformed_headers = true;
			}
		}
		if(malformed_headers) {
			double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
			res = sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, &t10, &t11, &t12, &t13, &t14, &t15);
			if(res == 15) {
				error_message << "\nSince the line contains 15 floats, likely the configuration contains more particles than specified in the topology file, or the header has been trimmed away.";
			}
			throw oxDNAException(error_message.str().c_str());
		}

		std::getline(_conf_input, line);
	}

	_box->init(Lx, Ly, Lz);

	// the following part is always carried out in double precision since we want to be able to restart from confs that have
	// large numbers in the conf file and use float precision later
	int k, i;
	std::vector<int> nins(_N_strands, 0);
	std::vector<LR_vector> scdm(_N_strands, LR_vector((double) 0., (double) 0., (double) 0.));

	i = 0;
	std::string line;
	while(!_conf_input.eof() && i < N()) {
		BaseParticle *p = _particles[i];

		if(!binary) {
			std::getline(_conf_input, line);
			auto spl_line = Utils::split_to_numbers(line, " ");

			p->pos = LR_vector(spl_line[0], spl_line[1], spl_line[2]);
			p->orientation.v1 = LR_vector(spl_line[3], spl_line[4], spl_line[5]);
			p->orientation.v3 = LR_vector(spl_line[6], spl_line[7], spl_line[8]);

			// get v2 from v1 and v3 and orthonormalise
			p->orientation.v1.normalize();
			p->orientation.v3.normalize();
			p->orientation.v1 -= p->orientation.v3 * (p->orientation.v1 * p->orientation.v3);
			p->orientation.v1.normalize();
			p->orientation.v2 = p->orientation.v3.cross(p->orientation.v1);
			p->orientation.v2.normalize();

			if(spl_line.size() == 15) {
				// read the momenta
				p->vel = LR_vector(spl_line[9], spl_line[10], spl_line[11]);
				p->L = LR_vector(spl_line[12], spl_line[13], spl_line[14]);
			}
		}
		else {
			p->pos = _read_next_binary_vector();

			int x, y, z;
			_conf_input.read((char*) &x, sizeof(int));
			_conf_input.read((char*) &y, sizeof(int));
			_conf_input.read((char*) &z, sizeof(int));
			p->set_pos_shift(x, y, z);

			p->orientation.v1 = _read_next_binary_vector();
			p->orientation.v2 = _read_next_binary_vector();
			p->orientation.v3 = _read_next_binary_vector();

			p->vel = _read_next_binary_vector();
			p->L = _read_next_binary_vector();
		}

		// v1, v2 and v3 should have length 1. If they don't it means that they are null vectors
		if(p->orientation.v1.module() < 0.9 || p->orientation.v2.module() < 0.9 || p->orientation.v3.module() < 0.9) {
			throw oxDNAException("Invalid orientation for particle %d: at least one of the vectors is a null vector", p->index);
		}
		p->orientation.transpone();

		k = p->strand_id;
		scdm[k] += p->pos;
		nins[k]++;

		p->init();
		p->orientationT = p->orientation.get_transpose();
		p->set_positions();

		i++;
	}

	// discarding the final '\n' in the binary file...
	if(binary && !_conf_input.eof()) {
		char tmpc;
		_conf_input.read((char*) &tmpc, sizeof(char));
	}

	if(i != N()) {
		if(_confs_to_skip > 0) {
			throw oxDNAException("Wrong number of particles (%d) found in configuration. Maybe you skipped too many configurations?", i);
		}
		else {
			throw oxDNAException("The number of lines found in configuration file (%d) doesn't match the parsed number of particles (%d)", i, N());
		}
	}

	for(k = 0; k < _N_strands; k++) {
		scdm[k] /= (double) nins[k];
	}

	for(i = 0; i < N(); i++) {
		BaseParticle *p = _particles[i];
		if(_enable_fix_diffusion && !binary) {
			// we need to manually set the particle shift so that the particle absolute position is the right one
			_box->shift_particle(p, scdm[p->strand_id]);
		}
	}

	_interaction->check_input_sanity(_particles);

	return true;
}

void SimBackend::_on_T_update() {
	_T = _config_info->temperature();
}

void SimBackend::apply_simulation_data_changes() {

}

void SimBackend::apply_changes_to_simulation_data() {

}

void SimBackend::add_output(ObservableOutputPtr new_output) {
	_obs_outputs.push_back(new_output);
}

void SimBackend::remove_output(std::string output_file) {
	auto search = std::find_if(_obs_outputs.begin(), _obs_outputs.end(), [output_file](ObservableOutputPtr p) { return p->get_output_name() == output_file; });
	if(search == _obs_outputs.end()) {
		throw oxDNAException("The output '%s' does not exist, can't remove it", output_file.c_str());
	}

	_obs_outputs.erase(search);
}

void SimBackend::print_observables() {
	_mytimer->resume();
	_obs_timer->resume();

	bool someone_ready = false;
	for(auto const &element : _obs_outputs) {
		if(element->is_ready(current_step()))
			someone_ready = true;
	}

	if(someone_ready) {
		apply_simulation_data_changes();

		llint total_bytes = 0;
		for(auto const &element : _obs_outputs) {
			if(element->is_ready(current_step())) {
				element->print_output(current_step());
			}
			total_bytes += element->get_bytes_written();
		}

		// here we control the timings; we leave the code a 30-second grace time to print the initial configuration
		double time_passed = (double) _mytimer->get_time() / (double) CLOCKS_PER_SEC;
		if(time_passed > 30) {
			double MBps = (total_bytes / (1024. * 1024.)) / time_passed;
			OX_DEBUG("Current data production rate: %g MB/s (total bytes: %lld, time_passed: %g)" , MBps, total_bytes, time_passed);
			if(MBps > _max_io) {
				throw oxDNAException(
						"Aborting because the program is generating too much data (%g MB/s).\n\t\t\tThe current limit is set to %g MB/s;\n\t\t\tyou can change it by setting max_io=<float> in the input file.",
						MBps, _max_io);
			}
		}

		apply_changes_to_simulation_data();
	}

	_backend_info = std::string("");

	_obs_timer->pause();
	_mytimer->pause();
}

void SimBackend::update_observables_data() {
	_obs_timer->resume();

	bool updated = false;
	for(auto const &obs : _config_info->observables) {
		if(obs->need_updating(current_step())) {
			if(!updated && obs->require_data_on_CPU()) {
				apply_simulation_data_changes();
				updated = true;
			}

			obs->update_data(current_step());
		}
	}

	_obs_timer->pause();
}

void SimBackend::fix_diffusion() {
	if(!_enable_fix_diffusion) {
		return;
	}

	apply_simulation_data_changes();

	for(auto mol : _molecules) {
		mol->update_com();
	}

	static std::vector<LR_vector> stored_pos(N());
	static std::vector<LR_matrix> stored_or(N());
	static std::vector<LR_matrix> stored_orT(N());

	// we can't exclude the possibility that N() might change during the course of the simulation
	stored_pos.resize(N());
	stored_or.resize(N());
	stored_orT.resize(N());

	number E_before = _interaction->get_system_energy(_particles, _lists.get());
	for(int i = 0; i < N(); i++) {
		BaseParticle *p = _particles[i];
		stored_pos[i] = p->pos;
		stored_or[i] = p->orientation;
		stored_orT[i] = p->orientationT;
	}

	// change particle position and fix orientation matrix;
	for(int i = 0; i < N(); i++) {
		BaseParticle *p = _particles[i];
		if(_molecules[p->strand_id]->shiftable()) {
			_box->shift_particle(p, _molecules[p->strand_id]->com);
			p->orientation.orthonormalize();
			p->orientationT = p->orientation.get_transpose();
			p->set_positions();
		}
	}

	number E_after = _interaction->get_system_energy(_particles, _lists.get());
	if(_interaction->get_is_infinite() == true || fabs(E_before - E_after) > 1.e-6 * fabs(E_before)) {
		OX_LOG(Logger::LOG_INFO, "Fix diffusion went too far... %g, %g, %d, (%g>%g)", E_before, E_after, _interaction->get_is_infinite(), fabs(E_before - E_after), 1.e-6 * fabs(E_before));
		_interaction->set_is_infinite(false);
		for(int i = 0; i < N(); i++) {
			BaseParticle *p = _particles[i];
			LR_vector dscdm = -_molecules[p->strand_id]->com;
			_box->shift_particle(p, dscdm);
			p->orientation = stored_or[i];
			p->orientationT = stored_orT[i];
			p->set_positions();
		}

		// we try to normalize a few of them...
		std::vector<int> apply;
		for(auto molecule : _molecules) {
			if(drand48() < 0.005) {
				apply.push_back(1);
			}
			else
				apply.push_back(0);
		}

		for(int i = 0; i < N(); i++) {
			BaseParticle *p = _particles[i];
			if(apply[p->strand_id]) {
				p->orientation.orthonormalize();
				p->orientationT = p->orientation.get_transpose();
				p->set_positions();
			}
		}

		number E_after2 = _interaction->get_system_energy(_particles, _lists.get());
		if(_interaction->get_is_infinite() == true || fabs(E_before - E_after2) > 1.e-6 * fabs(E_before)) {
			OX_LOG(Logger::LOG_INFO, " *** Fix diffusion hopeless... %g, %g, %d, (%g>%g)", E_before, E_after2, _interaction->get_is_infinite(), fabs(E_before - E_after2), 1.e-6 * fabs(E_before));
			_interaction->set_is_infinite(false);
			for(int i = 0; i < N(); i++) {
				BaseParticle *p = _particles[i];
				if(apply[p->strand_id]) {
					p->pos = stored_pos[i];
					p->orientation = stored_or[i];
					p->orientationT = stored_orT[i];
					p->set_positions();
				}
			}
		}
		else {
			OX_LOG(Logger::LOG_INFO, "diffusion PARTIALLY fixed");
		}
	}
	else {
		OX_LOG(Logger::LOG_INFO, "diffusion fixed.. %g %g", E_before, E_after);
	}

	apply_changes_to_simulation_data();
}

void SimBackend::print_conf(bool reduced, bool only_last) {
	apply_simulation_data_changes();

	if(reduced) {
		std::stringstream conf_name;
		conf_name << _reduced_conf_output_dir << "/reduced_conf" << current_step() << ".dat";
		_obs_output_reduced_conf->change_output_file(conf_name.str().c_str());
		_obs_output_reduced_conf->print_output(current_step());
	}
	else {
		if(!only_last) {
			_obs_output_trajectory->print_output(current_step());
		}
		_obs_output_last_conf->print_output(current_step());
		if(_obs_output_last_conf_bin != nullptr) {
			_obs_output_last_conf_bin->print_output(current_step());
		}

		// serialise the observables
		for(auto const &obs : _config_info->observables) {
			obs->serialise();
		}
	}
}

void SimBackend::print_equilibration_info() {
	// whoever overloads this will print something;
	return;
}
