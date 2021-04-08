/*
 * MCBackend.cpp
 *
 *  Created on: 25/nov/2010
 *      Author: lorenzo
 */

#include "MCBackend.h"
#include "IOManager.h"

template<typename number>
MCBackend<number>::MCBackend(IOManager *IO) : SimBackend<number>(IO), _overlap(false), _check_energy_counter(0) {
	this->_sim_type = SIM_MC;
	for(int i = 0; i < MC_MOVES; i++) {
		_tries[i] = _accepted[i] = 0;
	}
}

template<typename number>
MCBackend<number>::~MCBackend() {

}

template<>
void MCBackend<float>::_get_number_settings(input_file &inp) {
	if(getInputFloat(&inp, "check_energy_threshold", &_check_energy_threshold, 0) == KEY_NOT_FOUND)
		_check_energy_threshold = 1e-4f;
	getInputFloat(&inp, "delta_translation", &_delta[MC_MOVE_TRANSLATION], 1);
	if(_ensemble == MC_ENSEMBLE_NPT) getInputFloat(&inp, "delta_volume", &_delta[MC_MOVE_VOLUME], 1);
}

template<>
void MCBackend<double>::_get_number_settings(input_file &inp) {
	if(getInputDouble(&inp, "check_energy_threshold", &_check_energy_threshold, 0) == KEY_NOT_FOUND)
		_check_energy_threshold = 1e-4;
	getInputDouble(&inp, "delta_translation", &_delta[MC_MOVE_TRANSLATION], 1);
	getInputDouble(&inp, "delta_rotation", &_delta[MC_MOVE_ROTATION], 1);
	if(_ensemble == MC_ENSEMBLE_NPT) getInputDouble(&inp, "delta_volume", &_delta[MC_MOVE_VOLUME], 1);
}

template<typename number>
void MCBackend<number>::get_settings(input_file &inp) {
	char tmp[256];

	SimBackend<number>::get_settings(inp);

	if(getInputInt(&inp, "check_energy_every", &_check_energy_every, 0) == KEY_NOT_FOUND)
		_check_energy_every = 10;
	getInputString(&inp, "ensemble", tmp, 1);
	if(strncasecmp(tmp, "npt", 256) == 0) {
		_ensemble = MC_ENSEMBLE_NPT;
		this->_IO->die("Ensemble NPT will be implemented soon");
	}
	else if(strncasecmp(tmp, "nvt", 256) == 0) _ensemble = MC_ENSEMBLE_NVT;
	else this->_IO->die("Ensemble '%s' not supported\n", tmp);

	_get_number_settings(inp);
}

template<typename number>
//void MCBackend<number>::init(ifstream &conf_input) {
void MCBackend<number>::init(char conf_filename[256]) {
	//SimBackend<number>::init(conf_input);
	SimBackend<number>::init(conf_filename);
}

template<typename number>
void MCBackend<number>::print_energy(llint curr_step) {
	if(_check_energy_counter == _check_energy_every) {
		_check_energy_counter = 1;
		number old_energy = this->_U;
		_compute_energy();
		this->_IO->debug("Checking energy: %lf vs %lf", old_energy, this->_U);
		if(fabs((this->_U - old_energy)/old_energy) > _check_energy_threshold) this->_IO->die("Check of the energy failed: %f != %f (and the relative difference is %f, more than the threshold %f)", old_energy, this->_U, fabs((this->_U - old_energy)/old_energy), this->_check_energy_threshold);
		if(_overlap == true) this->_IO->die("Check of the energy failed: there is an overlap in the configuration");
	}
	else _check_energy_counter++;

	this->_IO->print_energy(*this, curr_step);
}

template<typename number>
void MCBackend<number>::print_conf(llint curr_step, bool reduced, bool only_last) {
	if(reduced == false) this->_IO->print_conf(*this, curr_step, only_last);
	else this->_IO->print_reduced_conf(*this, curr_step);
}

template class MCBackend<float>;
template class MCBackend<double>;
