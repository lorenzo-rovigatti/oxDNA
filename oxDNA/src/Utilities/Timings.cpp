/*
 * timings.cpp
 *
 *  Created on: 21/mar/2014
 *      Author: flavio 
 */

#include "Timings.h"
#include "oxDNAException.h"
#include "Logger.h"

#include <algorithm>

Timer::Timer() {
	_time = (clock_t) 0;
	_last = (clock_t) 0;
	//_last = clock();
	_active = false;
	_desc = std::string("Uninitialized timer");
}

Timer::Timer(std::string arg) {
	_desc = std::string(arg);
	_time = (clock_t) 0;
	_last = (clock_t) 0;
	//_last = clock();
	_active = false;
}

Timer::~Timer() {

}

void Timer::resume() {
	if(_active) throw oxDNAException("resuming already active timer %s", _desc.c_str());
	_last = clock();
	_active = true;
}

void Timer::pause() {
	if(!_active) throw oxDNAException("pausing resuming already inactive timer %s", _desc.c_str());
	_time += (clock() - _last);
	_active = false;
}

// this should work regardless of the timers being active
long long int Timer::get_time() {
	if(_active)
		return (long long int) (_time + (clock() - _last));
	else
		return (long long int) _time;
}

/***************** END OF TIMER CLASS *********************/

// singleton
TimingManager * TimingManager::_timingManager = NULL;

// time manager class
TimingManager::TimingManager() {

}

TimingManager::~TimingManager() {
	for(unsigned int i = 0; i < _timers.size(); i++) delete _timers[i];
}

void TimingManager::init() {
	if(_timingManager != NULL) throw oxDNAException("initializing an already initialized TimingManager");
	_timingManager = new TimingManager();
}

void TimingManager::clear() {
	if(_timingManager != NULL) delete _timingManager;
}

TimingManager * TimingManager::instance() {
	if(_timingManager == NULL) throw oxDNAException("accessing uninitialized TimingManager");
	return _timingManager;
}

void TimingManager::add_timer(Timer * arg) {
	_timers.push_back(arg);
	_parents.insert(std::make_pair(arg, (Timer *) NULL));
	_desc_map.insert(std::make_pair(arg->get_desc(), arg));
}

Timer * TimingManager::new_timer(std::string desc) {
	if(_desc_map.count(desc) != 0) throw oxDNAException("timer %s already used! Aborting", desc.c_str());

	Timer * timer = new Timer(desc);

	_timers.push_back(timer);
	//_parents.insert (std::make_pair (timer, (Timer *)NULL));
	_parents[timer] = (Timer *) NULL;
	//_desc_map.insert (std::make_pair (desc, timer));
	_desc_map[desc] = timer;

	OX_DEBUG("Adding new timer with description %s and no parent", desc.c_str());

	return timer;
}

Timer * TimingManager::new_timer(std::string desc, std::string parent_desc) {
	if(_desc_map.count(desc) != 0) throw oxDNAException("timer %s already used! Aborting", desc.c_str());
	if(_desc_map.count(parent_desc) == 0) throw oxDNAException("Cannot add timer %s because parent timer %s does not exist", desc.c_str(), parent_desc.c_str());

	Timer * timer = new Timer(desc);
	_timers.push_back(timer);
	//_parents.insert (std::make_pair (timer, get_timer_by_desc(parent_desc)));
	_parents[timer] = get_timer_by_desc(parent_desc);
	//_desc_map.insert (std::make_pair (desc, timer));
	_desc_map[desc] = timer;

	OX_DEBUG("Adding new timer with description %s and parent %s", desc.c_str(), parent_desc.c_str());

	return timer;
}

void TimingManager::add_timer(Timer * arg, std::string parent_desc) {
	std::string my_parent_desc;
	Timer * my_parent_ptr;
	if(_desc_map.count(parent_desc) > 0) {
		my_parent_desc = std::string(parent_desc);
		my_parent_ptr = _desc_map[parent_desc];
	}
	else {
		OX_LOG(Logger::LOG_WARNING, "Trying to add timer \"%s\" with an unknown parent \"%s\". Setting parent to \"None\"", arg->get_desc().c_str(), parent_desc.c_str());
		my_parent_desc = std::string("None");
		my_parent_ptr = NULL;
	}

	_timers.push_back(arg);
	_parents.insert(std::make_pair(arg, my_parent_ptr));
	_desc_map.insert(std::make_pair(arg->get_desc(), arg));
}

void TimingManager::print() {
	// times (including children) 
	std::map<Timer *, long long int> totaltimes;
	for(unsigned int i = 0; i < _timers.size(); i++) totaltimes[_timers[i]] = _timers[i]->get_time();

	// times in children 
	std::map<Timer *, long long int> sum_of_children;
	for(unsigned int i = 0; i < _timers.size(); i++) sum_of_children[_timers[i]] = 0;
	for(unsigned int i = 0; i < _timers.size(); i++) {
		Timer * t = _timers[i];
		Timer * p = _parents[t];
		while(p != NULL) {
			sum_of_children[p] += totaltimes[t];
			p = _parents[p];
		}
	}

	// own time (not in children)
	std::map<Timer *, long long int> own_time;
	for(unsigned int i = 0; i < _timers.size(); i++) {
		Timer * t = _timers[i];
		own_time[t] = totaltimes[t] - sum_of_children[t];
	}

	// mylist will be ordered as a tree
	std::vector<std::string> mylist;
	while(mylist.size() < _timers.size()) {
		for(unsigned int i = 0; i < _timers.size(); i++) {
			Timer * t = _timers[i];
			Timer * p = _parents[t];

			if(p == NULL) {
				mylist.push_back(t->get_desc());
			}
			else {
				// troviamo il nome del parente
				std::vector<std::string>::iterator it = std::find(mylist.begin(), mylist.end(), p->get_desc());
				if(it != mylist.end()) {
					it++;
					mylist.insert(it, t->get_desc());
				}
			}
		}
	}

	// now the list is ordered in the order we want to print it
	OX_LOG(Logger::LOG_NOTHING, "");
	OX_LOG(Logger::LOG_INFO, "Timings, in seconds, by Timer (total, own, spent in children)");
	for(unsigned int i = 0; i < mylist.size(); i++) {
		char mystr[512] = "";
		Timer * t = get_timer_by_desc(mylist[i]);
		Timer * p = _parents[t];
		int generations = 0;
		while(p != NULL) {
			generations++;
			p = _parents[p];
		}
		for(int j = 0; j < generations; j++) {
			strcat(mystr, "***");
		}
		strcat(mystr, "> ");
		strcat(mystr, t->get_desc().c_str());
		//printf ("%s %lld %lld %lld %s\n", t->get_desc().c_str(), totaltimes[t], own_time[t], sum_of_children[t]);
		OX_LOG(Logger::LOG_NOTHING, "%-30s %12.3lf %12.3lf %12.3lf", (char *) mystr, totaltimes[t] / CPSF, own_time[t] / CPSF, sum_of_children[t] / CPSF);
	}
	OX_LOG(Logger::LOG_NOTHING, "");

	return;
}

