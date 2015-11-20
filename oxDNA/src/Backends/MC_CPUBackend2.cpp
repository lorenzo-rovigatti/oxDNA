/*
 * MC_CPUBackend2.cpp
 *
 *  Created on: 02/apr/2014
 *      Author: Flavio 
 */

#include "MC_CPUBackend2.h"
#include <sstream>

template<typename number>
MC_CPUBackend2<number>::MC_CPUBackend2() : MCBackend<number>() {
	this->_is_CUDA_sim = false;
	_info_str = std::string("dummy");
	_N_moves = -1;
	_MC_Info = ConfigInfo<number>();
	_accumulated_prob = 0.; // total weight
	
}

template<typename number>
MC_CPUBackend2<number>::~MC_CPUBackend2() {

}

template<typename number>
void MC_CPUBackend2<number>::add_move (std::string move_string, input_file &sim_inp) {
	input_file * move_inp = Utils::get_input_file_from_string(move_string);

	BaseMove<number> * new_move = MoveFactory::make_move<number> (*move_inp, sim_inp, &_MC_Info);

	_moves.push_back (new_move);

	cleanInputFile(move_inp);
	delete move_inp;
}

template<typename number>
void MC_CPUBackend2<number>::get_settings(input_file &inp) {
	MCBackend<number>::get_settings(inp);

	//_MC_Info = ConfigInfo<number>(this->_particles, &(this->_box_side), this->_interaction, &(this->_N), &_info_str, this->_lists);
	_MC_Info = ConfigInfo<number>(this->_particles, &(this->_box_side), this->_interaction, &(this->_N), &_info_str, this->_lists, this->_box);

	//_MC_Info.lists = this->_lists;

	std::vector<std::string> move_strings;
	_N_moves = getInputKeys (&inp, std::string("move_"), &move_strings, 0); 
	if (_N_moves < 1) throw oxDNAException ("(MC_CPUBackend2) No moves found in the input file");
	for (int i = 0; i < _N_moves; i ++) {
		std::string tmps;
		getInputString (&inp, move_strings[i].c_str(), tmps, 1);
		add_move (tmps, inp);
	}

	typename vector<BaseMove<number> *>::iterator it;
	//for(it = _moves.begin(); it != _moves.end(); it ++) (*it)->get_settings(inp);
	//for(it = _moves.begin(); it != _moves.end(); it ++) (*it)->get_settings(inp);

	for(it = _moves.begin(); it != _moves.end(); it ++) _accumulated_prob += (*it)->prob;

	OX_LOG (Logger::LOG_INFO, "(MC_CPUBackend2.cpp) accumulated prob: %g", _accumulated_prob);
	
	//printf ("aborting here %s %d...\n", __FILE__, __LINE__);
	//abort ();
}

template<typename number>
void MC_CPUBackend2<number>::init() {
	OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend2) Initializing bakend...");
	MCBackend<number>::init();

	for(int i = 0; i < this->_N; i++) {
		// this is needed for the first _compute_energy()
		this->_particles[i]->set_positions();
		this->_particles[i]->orientationT = this->_particles[i]->orientation.get_transpose();
	}

	this->_U = this->_interaction->get_system_energy (this->_particles, this->_N, this->_lists);
	if(this->_interaction->get_is_infinite() == true) throw oxDNAException("There is an overlap in the initial configuration. Aborting");

	// needed to fill un the pointers....
	_MC_Info.particles = this->_particles;
	_MC_Info.N = &this->_N;
	_MC_Info.interaction = this->_interaction;
	_MC_Info.box_side = &(this->_box_side);
	_MC_Info.box = this->_box;

	//this->update_lists;
	this->_lists->global_update();

	// we initialize the moves
	typename vector<BaseMove<number> *>::iterator it;
	for(it = _moves.begin(); it != _moves.end(); it ++) {
		OX_LOG(Logger::LOG_INFO, "(MC_CPUBackend2) Initializing move...");
		(*it)->init();
	}
}

template<typename number>
void MC_CPUBackend2<number>::sim_step(llint curr_step) {
	for(int i = 0; i < this->_N; i++) {
		// pick a move with a given probability
		number choiche = drand48() * _accumulated_prob;
		int j = 0;
		number tmp = _moves[0]->prob;
		while (choiche > tmp) {
			j ++;
			tmp += _moves[j]->prob;
		}

		// now j is the chosen move
		_moves[j]->apply (curr_step);

		//for (int j=0; j < _N_moves; j ++) {
		//	_moves[j]->apply(curr_step);
		//}
	}
}

template<typename number>
void MC_CPUBackend2<number>::print_observables(llint curr_step) {
	std::string tmpstr("");
	typename std::vector<BaseMove<number> *>::iterator it; 
	for(it = _moves.begin(); it != _moves.end(); it ++) {
		number ratio = (*it)->get_acceptance();
		tmpstr += Utils::sformat (" %5.3f", ratio);
	}
	this->_backend_info.insert(0, tmpstr + "  ");

	SimBackend<number>::print_observables(curr_step);
}


template<typename number>
void MC_CPUBackend2<number>::_compute_energy() {
	printf ("somebody called me...\n");
	abort ();
}

template class MC_CPUBackend2<float>;
template class MC_CPUBackend2<double>;

