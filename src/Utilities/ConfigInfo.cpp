/*
 * ConfigInfo.cpp
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#include "ConfigInfo.h"

#include "oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Interactions/BaseInteraction.h"

std::shared_ptr<ConfigInfo> ConfigInfo::_config_info = nullptr;

ConfigInfo::ConfigInfo(std::vector<BaseParticle *> *ps, std::vector<std::shared_ptr<Molecule>> *mols) :
				particles_pointer(ps),
				molecules_pointer(mols) {

}

ConfigInfo::~ConfigInfo() {

}

void ConfigInfo::update_temperature(number new_T) {
	_temperature = new_T;
	notify("T_updated");
}

void ConfigInfo::subscribe(std::string event, std::function<void()> callback) {
	_event_callbacks[event].emplace_back(callback);
}

void ConfigInfo::notify(std::string event) {
	for(auto callback : _event_callbacks[event]) {
		callback();
	}
}

void ConfigInfo::set(BaseInteraction *i, std::string *info, BaseList *l, BaseBox *abox) {
	interaction = i;
	backend_info = info;
	lists = l;
	box = abox;
}

void ConfigInfo::init(std::vector<BaseParticle *> *ps, std::vector<std::shared_ptr<Molecule>> *mols) {
	if(_config_info != nullptr) {
		throw oxDNAException("The ConfigInfo object have been already initialised");
	}

	_config_info = std::shared_ptr<ConfigInfo>(new ConfigInfo(ps, mols));
}

void ConfigInfo::clear() {
	_config_info.reset();
	_config_info = nullptr;
}

const FlattenedConfigInfo &ConfigInfo::flattened_conf() {
	_flattened_conf.update(curr_step, particles());
	return _flattened_conf;
}
