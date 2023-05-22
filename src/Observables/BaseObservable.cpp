/*
 * BaseObservable.cpp
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#include "BaseObservable.h"

BaseObservable::BaseObservable() :
				_config_info(ConfigInfo::instance().get()) {
}

BaseObservable::~BaseObservable() {

}

bool BaseObservable::is_update_every_set() {
	return _update_every > 0;
}

bool BaseObservable::need_updating(llint curr_step) {
	return (_update_every > 0 && (curr_step % _update_every) == 0);
}

void BaseObservable::update_data(llint curr_step) {

}

void BaseObservable::get_settings(input_file &my_inp, input_file &sim_inp) {
	getInputString(&my_inp, "id", _id, 0);
	getInputLLInt(&my_inp, "update_every", &_update_every, 0);
}

void BaseObservable::init() {

}
