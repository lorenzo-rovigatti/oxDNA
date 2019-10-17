/*
 * MC_CPUBackend2.cpp
 *
 *  Created on: 02/apr/2014
 *      Author: Flavio
 */

#include "MC_CPUBackend2.h"
#include <sstream>

MC_CPUBackend2::MC_CPUBackend2() :
				MCBackend() {
	_is_CUDA_sim = false;
	_info_str = std::string("dummy");
	_N_moves = -1;
	_MC_Info = NULL;
	_accumulated_prob = 0.; // total weight

	_MC_Info = ConfigInfo::instance();
}

MC_CPUBackend2::~MC_CPUBackend2() {
	for(typename std::vector<BaseMove *>::iterator it = _moves.begin(); it != _moves.end(); it++)
		delete *it;
}

void MC_CPUBackend2::add_move(std::string move_string, input_file &sim_inp) {
	input_file * move_inp = Utils::get_input_file_from_string(move_string);

	BaseMove * new_move = MoveFactory::make_move(*move_inp, sim_inp);

	_moves.push_back(new_move);

	cleanInputFile(move_inp);
	delete move_inp;
}

void MC_CPUBackend2::get_settings(input_file &inp) {
	MCBackend::get_settings(inp);

	//_MC_Info = ConfigInfo(_particles, &(_box_side), _interaction, &(_N), &_info_str, _lists);
	_MC_Info->set(_particles, _interaction, &(_N), &_info_str, _lists.get(), _box.get());

	//_MC_Info.lists = _lists;

	std::vector<std::string> move_strings;
	_N_moves = getInputKeys(&inp, std::string("move_"), &move_strings, 0);
	if(_N_moves < 1) throw oxDNAException("(MC_CPUBackend2) No moves found in the input file");
	for(int i = 0; i < _N_moves; i++) {
		std::string tmps;
		getInputString(&inp, move_strings[i].c_str(), tmps, 1);
		add_move(tmps, inp);
	}

	typename vector<BaseMove *>::iterator it;
	//for(it = _moves.begin(); it != _moves.end(); it ++) (*it)->get_settings(inp);
	//for(it = _moves.begin(); it != _moves.end(); it ++) (*it)->get_settings(inp);

	for(it = _moves.begin(); it != _moves.end(); it++)
		_accumulated_prob += (*it)->prob;

	OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend2.cpp) accumulated prob: %g", _accumulated_prob);

	//printf ("aborting here %s %d...\n", __FILE__, __LINE__);
	//abort ();
}

void MC_CPUBackend2::init() {
	//OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend2) Initializing backend...");
	MCBackend::init();

	for(int i = 0; i < _N; i++) {
		// this is needed for the first _compute_energy()
		_particles[i]->set_positions();
		_particles[i]->orientationT = _particles[i]->orientation.get_transpose();
	}

	// needed to fill un the pointers....
	_MC_Info->particles = _particles;
	_MC_Info->N = &_N;
	_MC_Info->interaction = _interaction;
	_MC_Info->box = _box.get();

	_lists->global_update();

	_U = _interaction->get_system_energy(_particles, _N, _lists.get());
	if(_interaction->get_is_infinite() == true) {
		_interaction->set_is_infinite(false);
		for(int i = 0; i < _N; i++) {
			for(int j = 0; j < i; j++) {
				BaseParticle * p = _particles[i];
				BaseParticle * q = _particles[j];
				_interaction->pair_interaction(p, q);
				if(_interaction->get_is_infinite() == true) {
					OX_LOG(Logger::LOG_INFO, "   overlap %d %d", i, j);
					//printf ("%g %g %g\n", p->pos.x, p->pos.y, p->pos.z);
					//printf ("%g %g %g %g %g %g\n", q->pos.x, q->pos.y, q->pos.z, q->orientation.v3.x, q->orientation.v3.y, q->orientation.v3.z);
					_interaction->set_is_infinite(false);
				}
			}
		}
		throw oxDNAException("(MC_CPUBackend2) There is an overlap in the initial configuration. Aborting");
	}

	// we initialize the moves
	typename vector<BaseMove *>::iterator it;
	for(it = _moves.begin(); it != _moves.end(); it++) {
		//OX_LOG(Logger::LOG_DEBUG, "(MC_CPUBackend2) Initializing move...");
		(*it)->init();
	}
}

void MC_CPUBackend2::sim_step(llint curr_step) {
	for(int i = 0; i < _N; i++) {
		// pick a move with a given probability
		number choice = drand48() * _accumulated_prob;
		int j = 0;
		number tmp = _moves[0]->prob;
		while(choice > tmp) {
			j++;
			tmp += _moves[j]->prob;
		}

		// now j is the chosen move
		_moves[j]->apply(curr_step);

	}
}

void MC_CPUBackend2::print_observables(llint curr_step) {
	std::string tmpstr("");
	typename std::vector<BaseMove *>::iterator it;
	for(it = _moves.begin(); it != _moves.end(); it++) {
		number ratio = (*it)->get_acceptance();
		tmpstr += Utils::sformat(" %5.3f", ratio);
	}
	_backend_info.insert(0, tmpstr + "  ");

	SimBackend::print_observables(curr_step);
}

void MC_CPUBackend2::print_equilibration_info() {
	typename std::vector<BaseMove *>::iterator it;
	for(it = _moves.begin(); it != _moves.end(); it++) {
		(*it)->log_parameters();
	}
}

void MC_CPUBackend2::_compute_energy() {
	printf("somebody called me...\n");
	abort();
}
