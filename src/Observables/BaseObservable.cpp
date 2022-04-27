/*
 * BaseObservable.cpp
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#include "BaseObservable.h"

BaseObservable::BaseObservable() :
				_config_info(ConfigInfo::instance()) {
}

BaseObservable::~BaseObservable() {

}

void BaseObservable::get_settings(input_file &my_inp, input_file &sim_inp) {
	getInputString(&my_inp, "id", _id, 0);
}

void BaseObservable::init() {

}
