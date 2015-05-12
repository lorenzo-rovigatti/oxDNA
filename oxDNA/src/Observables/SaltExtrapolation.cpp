/*
 * SaltExtrapolation.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: petr
 */

#include "SaltExtrapolation.h"
//#include "../Interactions/RNAInteraction_V1.h"
#include <sstream>
#include <map>

#include "../Interactions/InteractionFactory.h"
#include "../Interactions/DNA2Interaction.h"

template<typename number>
SaltExtrapolation<number>::SaltExtrapolation() {
	_skip_zeros = false;
}

template<typename number>
SaltExtrapolation<number>::~SaltExtrapolation() {
	unsigned int i, j;
	for (i = 0; i < _temps.size(); i ++) {
		for (j = 0; j < _salts.size(); j ++) {
			delete _interactions[i][j];
		}
	}
}

template<typename number>
void SaltExtrapolation<number>::init (ConfigInfo<number> &config_info) {
   BaseObservable<number>::init(config_info);

   ifstream tmpf(_op_file.c_str());
   if(!tmpf.good ()) throw oxDNAException ("(SaltExtrapolation.cpp) Can't read file '%s'", _op_file.c_str());

   _op.set_log_level (Logger::LOG_DEBUG);
   _op.init_from_file (_op_file.c_str(), this->_config_info.particles, *(this->_config_info.N));
	
   _weights.init((const char *) _weights_file.c_str(), &_op, false, 1.);
   
   // get the dimension of the order parameter
   int n_op_dim = _op.get_all_parameters_count();
   int * sizes = new int[n_op_dim];
   memcpy (sizes, _op.get_state_sizes(), n_op_dim * sizeof (int));
   int n_op_states = 1;
   for (int i = 0; i < n_op_dim; i ++) n_op_states *= sizes[i];
   delete [] sizes;
   
   _hists = std::vector< std::vector< std::vector <long double> > > (_temps.size(), std::vector< std::vector<long double> > (_salts.size(), std::vector<long double> (n_op_states, (long double) 0.)));
	
   OX_LOG(Logger::LOG_INFO, "(SaltExtrapolation.cpp) Estimate of memory occupied by histogram: ~ %g KB", _temps.size() * _salts.size() * n_op_states * sizeof(long double) / 1024.);

}


template<typename number>
void SaltExtrapolation<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	std::string tmps;

	getInputString (&my_inp, "salts", tmps, 1);
	std::vector<string> salt_strings = Utils::split(tmps, ',');
	for (unsigned int i = 0; i < salt_strings.size(); i ++) _salts.push_back((number) atof(salt_strings[i].c_str()));

	getInputString (&my_inp, "temps", tmps, 1);
	std::vector<string> temp_strings = Utils::split(tmps, ',');
	for (unsigned int i = 0; i < temp_strings.size(); i ++) _temps.push_back(Utils::get_temperature<number>((char*)temp_strings[i].c_str()));

	OX_LOG(Logger::LOG_INFO, "(SaltExtrapolation.cpp) Extrapolating to %d temperatures and %d salts...", _temps.size(), _salts.size());
	
	// set the simulation temperature
	getInputString (&sim_inp, "T", tmps, 1);
	_sim_temp = Utils::get_temperature<number> ((char *)tmps.c_str());

	// check that the interaction is DNA2
	getInputString (&sim_inp, "interaction_type", tmps, 1);
	if (tmps.compare ("DNA2") != 0) throw oxDNAException ("(SaltExtrapolation.cpp) Salt extrapolation only works with interaction_type = DNA2. Found interaction_type = %s. Aborting", tmps.c_str());

	// read the order parameter file with respect to which use the histograms...
	if (getInputString (&my_inp, "op_file", _op_file, 0) == KEY_NOT_FOUND) getInputString(&sim_inp, "op_file", _op_file, 1);
	if (getInputString (&my_inp, "weights_file", _weights_file, 0) == KEY_NOT_FOUND) getInputString(&sim_inp, "weights_file", _weights_file, 1);

	OX_LOG (Logger::LOG_INFO, "(SaltExtrapolation.cpp) Getting order parameter from file %s and weights from file %s...", _op_file.c_str(), _weights_file.c_str());

	// find the interaction with the largest cutoff
	_interactions = std::vector<std::vector <IBaseInteraction<number> *> > (_temps.size(), std::vector<IBaseInteraction<number> *> (_salts.size()));

	// we read the topology to set up the interactions...
	std::string topology;
	getInputString (&sim_inp, "topology", topology, 1);
	unsigned int i, j;
	number smaller_salt = _salts[0];
	number larger_temp = _temps[0];
	for (i = 0; i < _temps.size(); i ++) {
		for (j = 0; j < _salts.size(); j ++) {
			std::stringstream fake_input_stream;
			fake_input_stream << "interaction_type = DNA2" << endl;
			fake_input_stream << "topology = " << topology << endl;
			fake_input_stream << "T = " << _temps[i] << endl;
			fake_input_stream << "salt_concentration = " << _salts[j] << endl;

			input_file * fake_input = Utils::get_input_file_from_string(fake_input_stream.str());

			_interactions[i][j] = InteractionFactory::make_interaction<number>(*fake_input);
			_interactions[i][j]->get_settings (*fake_input);
			_interactions[i][j]->init();

			cleanInputFile (fake_input);
			delete fake_input;

			if (_salts[j] < smaller_salt) smaller_salt = _salts[j];
		}
		if (_temps[i] > larger_temp) larger_temp = _temps[i];
	}
	// we store the interaction that has the largest cutoff
	// [! maybe we should check that the real interaction does not have anything to do with it
	_ref_interaction = _interactions[larger_temp][smaller_salt];

}

template<typename number>
std::string SaltExtrapolation<number>::get_output_string(llint curr_step) {
	unsigned int i, j, k;
	std::string output_str;

	// we need to set the box side for all the interactions
	for (i = 0; i < _temps.size(); i++)
		for (j = 0; j < _salts.size(); j ++)
			_interactions[i][j]->set_box_side(*this->_config_info.box_side);
	
	// we get the potential interactions once from the interaction
	// with the largest cutoff
	std::vector<ParticlePair<number> > neighbour_pairs = this->_config_info.lists->get_potential_interactions();

	number e_sim = 0.;
	number e0 = 0.;
	std::vector<number> es = std::vector<number> (_temps.size(), 0.);
	std::vector< std::vector <number> > edhs = std::vector<std::vector <number> > (_temps.size(), std::vector<number> (_salts.size(), (number) 0.));
	_op.reset();
	_op.fill_distance_parameters(this->_config_info.particles,*this->_config_info.box_side);
	for (k = 0; k < neighbour_pairs.size(); k ++) {
		BaseParticle<number> * p = neighbour_pairs[k].first;
		BaseParticle<number> * q = neighbour_pairs[k].second;
	
		// interaction terms that do not depend on temperature... they are the
		// same for all the interactions
		e0 += this->_config_info.interaction->pair_interaction_term (DNA2Interaction<number>::BACKBONE, p, q, NULL, false);
		e0 += this->_config_info.interaction->pair_interaction_term (DNA2Interaction<number>::BONDED_EXCLUDED_VOLUME, p, q, NULL, false);
		e0 += this->_config_info.interaction->pair_interaction_term (DNA2Interaction<number>::NONBONDED_EXCLUDED_VOLUME, p, q, NULL, false);
		e0 += this->_config_info.interaction->pair_interaction_term (DNA2Interaction<number>::CROSS_STACKING, p, q, NULL, false);
		e0 += this->_config_info.interaction->pair_interaction_term (DNA2Interaction<number>::COAXIAL_STACKING, p, q, NULL, false);
		// the following line also updates the state of the order parameter
		number hb_energy = this->_config_info.interaction->pair_interaction_term(DNA2Interaction<number>::HYDROGEN_BONDING, p, q, NULL, false);
		if(hb_energy < HB_CUTOFF) _op.add_hb(p->index,q->index);
		e0 += hb_energy;
	
		// now we update the terms that depend on temperature and salt separately
		for (i = 0; i < _temps.size(); i ++ ) {
			// the following line uses _interactions[i][0] since all the interactions with the 
			// same i have the same temperature
			es[i] += _interactions[i][0]->pair_interaction_term (DNA2Interaction<number>::STACKING, p, q, NULL, false);
			for (j = 0; j < _salts.size(); j ++) edhs[i][j] += _interactions[i][j]->pair_interaction_term (DNA2Interaction<number>::DEBYE_HUCKEL, p, q, NULL, false);
		}

		// we update stacking and electrostatic for the simulation energy
		e_sim += this->_config_info.interaction->pair_interaction_term (DNA2Interaction<number>::STACKING, p, q, NULL, false);
		e_sim += this->_config_info.interaction->pair_interaction_term (DNA2Interaction<number>::DEBYE_HUCKEL, p, q, NULL, false);
	}
	
	// here we take care of the external energy
	number eext = 0.;
	for (k = 0; k < (unsigned int)*this->_config_info.N; k ++) {
		BaseParticle<number> *p = this->_config_info.particles[k];
		p->set_ext_potential(curr_step, *this->_config_info.box_side);
		eext += p->ext_potential;
	}

	// we get index and weight from the order parameter
	int idx;
	double wgt = _weights.get_weight (_op.get_all_states(), &idx);
	for (i = 0; i < _temps.size(); i ++)
		for (j = 0; j < _salts.size(); j ++)
			_hists[i][j][idx] += expl (-(e0 + es[i] + edhs[i][j] + eext) / _temps[i] + (e0 + e_sim + eext) / _sim_temp) / wgt;

	// for each salt, we print a histogram
	int n_op_dim = _op.get_all_parameters_count();
	int * sizes = new int[n_op_dim];
	memcpy (sizes, _op.get_state_sizes(), n_op_dim * sizeof (int));
	for (j = 0; j < _salts.size(); j ++) {
		// header per each salt
		output_str += Utils::sformat ("t: %llu, salt: %g, Ts: ", curr_step, _salts[j]);
		for (i = 0; i < _temps.size(); i ++) output_str += Utils::sformat ("%g ", _temps[i]);
		output_str += Utils::sformat ("\n");

		// now we cycle through the histogram values
		for (k = 0; k < _hists[0][j].size(); k ++) { // k runs over the op states
			for (int jj = 0; jj < n_op_dim; jj ++) { // this cycle gets the OP VALUE (such as 1, 0, 3 for a three-dim OP)
				int pindex = 1;
				for (int kk = 0; kk < jj; kk ++) pindex *= sizes[kk];
				output_str += Utils::sformat ("%d ", (k / pindex) % sizes[jj]);
			}
			for (i = 0; i < _temps.size(); i ++) output_str += Utils::sformat ("%Le ", _hists[i][j][k]); // j runs over the temperatures
			output_str += Utils::sformat("\n");
		}
	}
	output_str += Utils::sformat ("\n");
	output_str += Utils::sformat ("\n");
	
	delete [] sizes;

    return output_str;
}

template class SaltExtrapolation<float>;
template class SaltExtrapolation<double>;

