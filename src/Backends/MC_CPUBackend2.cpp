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
	_info_str = std::string("dummy");
	_N_moves = -1;
	_accumulated_prob = 0.; // total weight
}

MC_CPUBackend2::~MC_CPUBackend2() {

}

void MC_CPUBackend2::add_move(std::string move_string, input_file &sim_inp) {
	input_file *move_inp = Utils::get_input_file_from_string(move_string);

	MovePtr new_move = MoveFactory::make_move(*move_inp, sim_inp);
	_moves.push_back(new_move);
	delete move_inp;
}

void MC_CPUBackend2::get_settings(input_file &inp) {
	MCBackend::get_settings(inp);

	std::vector<std::string> move_strings;
	_N_moves = getInputKeys(&inp, std::string("move_"), &move_strings, 0);
	if(_N_moves < 1) throw oxDNAException("(MC_CPUBackend2) No moves found in the input file");
	for(int i = 0; i < _N_moves; i++) {
		std::string tmps;
		getInputString(&inp, move_strings[i].c_str(), tmps, 1);
		add_move(tmps, inp);
	}

	for(auto move : _moves) {
		_accumulated_prob += move->prob;
	}

	OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend2.cpp) accumulated prob: %g", _accumulated_prob);
}

void MC_CPUBackend2::init() {
	MCBackend::init();

	for(auto p: _particles) {
		// this is needed for the first _compute_energy()
		p->set_positions();
		p->orientationT = p->orientation.get_transpose();
	}

	_lists->global_update();

	_U = _interaction->get_system_energy(_particles, _lists.get());
	if(_interaction->get_is_infinite() == true) {
		_interaction->set_is_infinite(false);
		for(int i = 0; i < N(); i++) {
			for(int j = 0; j < i; j++) {
				BaseParticle *p = _particles[i];
				BaseParticle *q = _particles[j];
				_interaction->pair_interaction(p, q);
				if(_interaction->get_is_infinite() == true) {
					OX_LOG(Logger::LOG_INFO, "   overlap %d %d", i, j);
					_interaction->set_is_infinite(false);
				}
			}
		}
		throw oxDNAException("(MC_CPUBackend2) There is an overlap in the initial configuration. Aborting");
	}

	// we initialize the moves
	for(auto move : _moves) {
		move->init();
	}
}

void MC_CPUBackend2::sim_step() {
	_mytimer->resume();

	for(int i = 0; i < N(); i++) {
		// pick a move with a given probability
		number choice = drand48() * _accumulated_prob;
		int j = 0;
		number tmp = _moves[0]->prob;
		while(choice > tmp) {
			j++;
			tmp += _moves[j]->prob;
		}

		// now j is the chosen move
		_moves[j]->apply(current_step());
	}

	_mytimer->pause();
}

void MC_CPUBackend2::print_observables() {
	std::string tmpstr("");
	for(auto move : _moves) {
		number ratio = move->get_acceptance();
		tmpstr += Utils::sformat(" %5.3f", ratio);
	}
	_backend_info.insert(0, tmpstr + "  ");

	SimBackend::print_observables();
}

void MC_CPUBackend2::print_equilibration_info() {
	for(auto move : _moves) {
		move->log_parameters();
	}
}

void MC_CPUBackend2::_compute_energy() {
	throw oxDNAException("Somebody called MC_CPUBackend2's _compute_energy method...\n");
}
