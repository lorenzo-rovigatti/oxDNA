/*
 * ConfigInfo.cpp
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#include "ConfigInfo.h"

#include "oxDNAException.h"

ConfigInfo *ConfigInfo::_config_info = NULL;

ConfigInfo::ConfigInfo() :
				particles(NULL),
				interaction(NULL),
				N(NULL),
				backend_info(NULL),
				lists(NULL),
				box(NULL),
				curr_step(0) {

}

ConfigInfo::~ConfigInfo() {

}

void ConfigInfo::set(BaseParticle **p, IBaseInteraction *i, int *Nn, std::string *info, BaseList *l, BaseBox *abox) {
	particles = p;
	interaction = i;
	N = Nn;
	backend_info = info;
	lists = l;
	box = abox;
}

void ConfigInfo::init() {
	if(_config_info != NULL) throw oxDNAException("The ConfigInfo object have been already initialised");

	_config_info = new ConfigInfo();
}

void ConfigInfo::clear() {
	if(_config_info != NULL) delete _config_info;
}
