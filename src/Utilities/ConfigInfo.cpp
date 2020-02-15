/*
 * ConfigInfo.cpp
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#include "ConfigInfo.h"

#include "oxDNAException.h"

std::shared_ptr<ConfigInfo> ConfigInfo::_config_info = nullptr;

ConfigInfo::ConfigInfo(std::vector<BaseParticle *> &ps) :
				particles(ps) {

}

ConfigInfo::~ConfigInfo() {

}

void ConfigInfo::set(IBaseInteraction *i, std::string *info, BaseList *l, BaseBox *abox) {
	interaction = i;
	backend_info = info;
	lists = l;
	box = abox;
}

void ConfigInfo::init(std::vector<BaseParticle *> &ps) {
	if(_config_info != nullptr) {
		throw oxDNAException("The ConfigInfo object have been already initialised");
	}

	_config_info = std::shared_ptr<ConfigInfo>(new ConfigInfo(ps));
}
