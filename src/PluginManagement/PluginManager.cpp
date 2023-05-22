/*
 * PluginManager.cpp
 *
 *  Created on: 09/dic/2013
 *      Author: lorenzo
 */

#include "PluginManager.h"

#include <dlfcn.h>
#include <sys/stat.h>

#include "../Utilities/Logger.h"
#include "../Utilities/Utils.h"
#include "../Utilities/oxDNAException.h"

using std::string;
using std::vector;

std::shared_ptr<PluginManager> PluginManager::_manager = nullptr;

typedef BaseObservable* make_obs();
typedef BaseInteraction* make_inter();
typedef BaseMove* make_move();

PluginManager::PluginManager() :
				_initialised(false),
				_do_cleanup(true) {
	_path.push_back(string("."));
}

PluginManager::~PluginManager() {
	if(_manager != nullptr && _manager->_do_cleanup) {
		while(!_manager->_handles.empty()) {
			dlclose(_manager->_handles.top());
			_manager->_handles.pop();
		}
	}
}

void PluginManager::init(input_file &sim_inp) {
	_initialised = true;

	getInputBool(&sim_inp, "plugin_do_cleanup", &_do_cleanup, 0);

	// load the plugin path from the input file
	string ppath;
	if(getInputString(&sim_inp, "plugin_search_path", ppath, 0) == KEY_FOUND) {
		vector<string> paths = Utils::split(ppath, ':');
		for(vector<string>::iterator it = paths.begin(); it < paths.end(); it++)
			add_to_path(*it);
	}

	string entries;
	// for observables
	if(getInputString(&sim_inp, "plugin_observable_entry_points", entries, 0) == KEY_FOUND) {
		vector<string> v_entries = Utils::split(entries, ':');
		for(vector<string>::iterator it = v_entries.begin(); it < v_entries.end(); it++)
			_obs_entry_points.push_back(*it);
	}
	// for interactions
	if(getInputString(&sim_inp, "plugin_interaction_entry_points", entries, 0) == KEY_FOUND) {
		vector<string> v_entries = Utils::split(entries, ':');
		for(vector<string>::iterator it = v_entries.begin(); it < v_entries.end(); it++)
			_inter_entry_points.push_back(*it);
	}

	if(getInputString(&sim_inp, "plugin_move_entry_points", entries, 0) == KEY_FOUND) {
		vector<string> v_entries = Utils::split(entries, ':');
		for(vector<string>::iterator it = v_entries.begin(); it < v_entries.end(); it++)
			_move_entry_points.push_back(*it);
	}

	// these are the default entry point names
	_obs_entry_points.push_back(string("make"));
	_obs_entry_points.push_back(string("make_observable"));

	_inter_entry_points.push_back(string("make"));
	_inter_entry_points.push_back(string("make_interaction"));

	_move_entry_points.push_back(string("make"));
	_move_entry_points.push_back(string("make_move"));
}

void PluginManager::add_to_path(string s) {
	_path.push_back(s);
}

void *PluginManager::_get_handle(string &name) {
	if(!_initialised) throw oxDNAException("PluginManager not initialised, aborting");

	void *handle = NULL;
	vector<string> possible_paths;
	for(vector<string>::iterator it = _path.begin(); it != _path.end() && handle == NULL; it++) {
		string path = *it + "/" + name + ".so";
		OX_DEBUG("Looking for plugin '%s' in '%s'", name.c_str(), it->c_str());
		struct stat buffer;
		if(stat(path.c_str(), &buffer) == 0) {
			possible_paths.push_back(path);
		}
	}

	if(possible_paths.size() == 0) {
		OX_LOG(Logger::LOG_WARNING, "Plugin shared library '%s.so' not found", name.c_str());
	}
	else {
		if(possible_paths.size() > 1) {
			OX_LOG(Logger::LOG_WARNING, "Multiple (%d) plugin shared libraries named '%s.so' were found. Using %s", possible_paths.size(), name.c_str(), possible_paths[0].c_str());
		}

		handle = dlopen(possible_paths[0].c_str(), RTLD_LAZY);
		const char *dl_error = dlerror();
		if(dl_error != NULL) {
			throw oxDNAException("Caught an error while opening shared library '%s.so': %s", name.c_str(), dl_error);
		}

		_handles.push(handle);
	}

	return handle;
}

void *PluginManager::_get_entry_point(void *handle, string name, vector<string> entry_points) {
	void *res = NULL;

	// we add the make_NAME entry to the list of possible entry point names
	string named_entry = string("make_") + name;
	entry_points.insert(entry_points.begin(), named_entry);

	for(auto entry_point : entry_points) {
		res = dlsym(handle, entry_point.c_str());
		if(!dlerror()) {
			return res;
		}
	}

	return res;
}

ObservablePtr PluginManager::get_observable(string name) {
	void *handle = _get_handle(name);
	if(handle == NULL) return NULL;

	void *temp_obs;
	bool found = false;
	make_obs *make_new_obs = (make_obs *) _get_entry_point(handle, name, _obs_entry_points);
	if(make_new_obs != NULL) {
		temp_obs = (void *) make_new_obs();
		found = true;
	}

	if(!found) {
		OX_LOG(Logger::LOG_WARNING, "Cannot load symbol from plugin observable library '%s'", name.c_str());
		return NULL;
	}

	// now we cast it back to the type required by the code
	return ObservablePtr((BaseObservable *) temp_obs);
}

InteractionPtr PluginManager::get_interaction(string name) {
	void *handle = _get_handle(name);
	if(handle == NULL) return NULL;

	void *temp_inter;
	bool found = false;
	make_inter *make_new_inter = (make_inter *) _get_entry_point(handle, name, _inter_entry_points);
	if(make_new_inter != NULL) {
		temp_inter = (void *) make_new_inter();
		found = true;
	}

	if(!found) {
		OX_LOG(Logger::LOG_WARNING, "Cannot load symbol from plugin interaction library '%s'", name.c_str());
		return NULL;
	}

	// now we cast it back to the type required by the code
	return InteractionPtr((BaseInteraction *)temp_inter);
}

MovePtr PluginManager::get_move(std::string name) {
	void *handle = _get_handle(name);
	if(handle == NULL) return NULL;

	void *temp_move;
	bool found = false;
	make_move *make_new_move = (make_move *) _get_entry_point(handle, name, _move_entry_points);
	if(make_new_move != NULL) {
		temp_move = (void *) make_new_move();
		found = true;
	}

	if(!found) {
		OX_LOG(Logger::LOG_WARNING, "Cannot load symbol from plugin interaction library '%s'", name.c_str());
		return NULL;
	}

	// now we cast it back to the type required by the code
	return MovePtr((BaseMove *)temp_move);

}

std::shared_ptr<PluginManager> PluginManager::instance() {
	if(_manager == nullptr) {
		// we can't use std::make_shared because PluginManager's constructor is private
		_manager = std::shared_ptr<PluginManager>(new PluginManager());
	}

	return _manager;
}
