/*
 * Step.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "Step.h"
#include "../Utilities/Utils.h"

template<typename number>
Step<number>::Step() {
    _dt = (number) 1.f;
    _units = STEP_UNIT_ONE_PER_STEP;
}

template<typename number>
Step<number>::~Step() {

}

template<typename number>
std::string Step<number>::get_output_string(llint curr_step) {
	if (_units == STEP_UNIT_ONE_PER_STEP) return Utils::sformat("%12lld", curr_step);
	else return Utils::sformat("%14.4lf", curr_step * _dt);
}

template<typename number>
void Step<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
    // if no units are specified, we assume STEP_UNIT_ONE_PER_STEP
    char tmpstr[256];
    if (getInputString(&my_inp, "units", tmpstr, 0) == KEY_FOUND) {
    	if(!strncmp (tmpstr, "MD", 256)) _units = STEP_UNIT_HONOUR_DT;
    	else if (!strncmp (tmpstr, "steps", 256)) _units = STEP_UNIT_ONE_PER_STEP;
    	else throw oxDNAException("Could not understand syntax when parsing the configuration of the step observable. Aborting");
    }
    else _units = STEP_UNIT_ONE_PER_STEP;

    // we need to read dt; if it's not there, we assume dt=1, which is the
    // sensible choice since we don't want this to break for MC, for
    // example;
    if (_units != STEP_UNIT_ONE_PER_STEP) {
		double tmp_dt;
		getInputDouble (&sim_inp, "dt", &tmp_dt, 1);
		_dt = (number)tmp_dt;
    }
}

template class Step<float>;
template class Step<double>;
