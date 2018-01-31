/*
 * ConfigInfo.cpp
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#include "ConfigInfo.h"

#include "oxDNAException.h"

template<typename number>
ConfigInfo<number> *ConfigInfo<number>::_config_info = NULL;

template<typename number>
ConfigInfo<number>::ConfigInfo() :
				particles(NULL),
				interaction(NULL),
				N(NULL),
				backend_info(NULL),
				lists(NULL),
				box(NULL),
				curr_step(0) {

}

template<typename number>
ConfigInfo<number>::~ConfigInfo() {

}

template<typename number>
void ConfigInfo<number>::set(BaseParticle<number> **p, IBaseInteraction<number> *i, int *Nn, std::string *info, BaseList<number> *l, BaseBox<number> *abox) {
	particles = p;
	interaction = i;
	N = Nn;
	backend_info = info;
	lists = l;
	box = abox;
}

template<typename number>
void ConfigInfo<number>::init() {
	if(_config_info != NULL) throw oxDNAException("The ConfigInfo object have been already initialised");

	_config_info = new ConfigInfo();
}

template<typename number>
void ConfigInfo<number>::clear() {
	if(_config_info != NULL) delete _config_info;
}

template class ConfigInfo<float> ;
template class ConfigInfo<double> ;
