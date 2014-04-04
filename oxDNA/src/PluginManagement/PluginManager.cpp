/*
 * PluginManager.cpp
 *
 *  Created on: 09/dic/2013
 *      Author: lorenzo
 */

#include <dlfcn.h>

#include "PluginManager.h"

#include "../Utilities/Logger.h"
#include "../Utilities/Utils.h"
#include "../Utilities/oxDNAException.h"

using std::string;

PluginManager *PluginManager::_manager = NULL;

typedef BaseObservable<float>* make_float_obs();
typedef BaseObservable<double>* make_double_obs();

typedef IBaseInteraction<float>* make_float_inter();
typedef IBaseInteraction<double>* make_double_inter();

PluginManager::PluginManager() {
	_path.push_back(string("."));
	_initialised = false;
}

PluginManager::~PluginManager() {

}

void PluginManager::init(input_file &sim_inp) {
	_initialised = true;

	// load the plugin path from the input file
	char ppath[10000];
	if(getInputString(&sim_inp, "plugin_search_path", ppath, 0) == KEY_FOUND) {
		std::vector<string> paths = Utils::split(string(ppath), ':');
		for(std::vector<string>::iterator it = paths.begin(); it < paths.end(); it++) add_to_path(*it);
	}

}

void PluginManager::add_to_path(std::string s) {
	_path.push_back(s);
}

void *PluginManager::_get_handle(std::string &name) {
	if(!_initialised) throw oxDNAException("PluginManager not initialised, aborting");

	void *handle = NULL;
	for(std::vector<string>::iterator it = _path.begin(); it != _path.end() && handle == NULL; it++) {
		string path = *it + "/" + name + ".so";
		OX_DEBUG("Looking for plugin '%s' in '%s'", name.c_str(), it->c_str());
		handle = dlopen(path.c_str(), RTLD_LAZY);
	}
	if(!handle) throw oxDNAException("Shared library '%s.so' not found", name.c_str());

	_handles.push(handle);

	return handle;
}

template<typename number>
BaseObservable<number> *PluginManager::get_observable(string name) {
	void *handle = _get_handle(name);

	// we do this c-like because dynamic linking can be done only in c and thus
	// we have no way of using templates
	void *temp_obs;
	// choose between float and double
	const char *dlsym_error;
	if(sizeof(number) == 4) {
		make_float_obs *make_obs = (make_float_obs *) dlsym(handle, "make_float");
		dlsym_error = dlerror();
		if(!dlsym_error) temp_obs = (void *)make_obs();
	}
	else {
		make_double_obs *make_obs = (make_double_obs *) dlsym(handle, "make_double");
		dlsym_error = dlerror();
		if(!dlsym_error) temp_obs = (void *)make_obs();
	}

	if(dlsym_error) {
		OX_LOG(Logger::LOG_WARNING, "Cannot load symbol from plugin observable library '%s'", name.c_str());
		return NULL;
	}

	// now we cast it back to the type required by the code
	return static_cast<BaseObservable<number> *>(temp_obs);
}

template<typename number>
IBaseInteraction<number> *PluginManager::get_interaction(std::string name) {
	void *handle = _get_handle(name);

	// we do this c-like because dynamic linking can be done only in c and thus
	// we have no way of using templates
	void *temp_inter;
	// choose between float and double
	const char *dlsym_error;
	if(sizeof(number) == 4) {
		make_float_inter *make_inter = (make_float_inter *) dlsym(handle, "make_float");
		dlsym_error = dlerror();
		if(!dlsym_error) temp_inter = (void *)make_inter();
	}
	else {
		make_double_inter *make_inter = (make_double_inter *) dlsym(handle, "make_double");
		dlsym_error = dlerror();
		if(!dlsym_error) temp_inter = (void *)make_inter();
	}

	if(dlsym_error) {
		OX_LOG(Logger::LOG_WARNING, "Cannot load symbol from plugin interaction library '%s'", name.c_str());
		return NULL;
	}

	// now we cast it back to the type required by the code
	return static_cast<IBaseInteraction<number> *>(temp_inter);
}

void PluginManager::clear() {
	if(_manager != NULL) {
		while(!_manager->_handles.empty()) {
			dlclose(_manager->_handles.top());
			_manager->_handles.pop();
		}
		delete _manager;
	}
}

PluginManager *PluginManager::instance() {
	if(_manager == NULL) _manager = new PluginManager();

	return _manager;
}

template BaseObservable<float> *PluginManager::get_observable(string name);
template BaseObservable<double> *PluginManager::get_observable(string name);

template IBaseInteraction<float> *PluginManager::get_interaction(string name);
template IBaseInteraction<double> *PluginManager::get_interaction(string name);
