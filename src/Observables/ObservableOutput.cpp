/*
 * ObservableOutput.cpp
 *
 *  Created on: Feb 11, 2013
 *      Author: rovigatti
 */

#include <sstream>
#include <iostream>

#include "ObservableOutput.h"
#include "ObservableFactory.h"
#include "../Utilities/Utils.h"

using namespace std;

template<typename number>
ObservableOutput<number>::ObservableOutput(std::string &stream_string, input_file &sim_inp) : _prefix("") {
	_sim_inp = sim_inp;
	_start_from = 0;
	_stop_at = -1;
	input_file *obs_input = Utils::get_input_file_from_string(stream_string);

	string out_name;
	getInputString(obs_input, "name", out_name, 1);
	// the prefix should be used only when the output is nor stdout nor stderr
	if(out_name != "stdout" && out_name != "stderr") getInputString(&sim_inp, "output_prefix", _prefix, 0);
	sprintf(_output_name, "%s%s", _prefix.c_str(), out_name.c_str());
	getInputLLInt(obs_input, "print_every", &_print_every, 1);

	getInputLLInt(obs_input, "start_from", &_start_from, 0);
	getInputLLInt(obs_input, "stop_at", &_stop_at, 0);

	int restart_step_counter = 0;
	getInputBoolAsInt(&sim_inp, "restart_step_counter", &restart_step_counter, 0);
	_append = !restart_step_counter;

	_only_last = 0;
	getInputBoolAsInt(obs_input, "only_last", &_only_last, 0);

	_is_binary = false;
	int tmpi;
	if(getInputBoolAsInt(obs_input, "binary", &tmpi, 0) == KEY_FOUND) _is_binary = (tmpi > 0);

	int i = 1;
	bool found = true;
	while(found) {
		stringstream ss;
		ss << "col_" << i;
		string obs_string;
		if(getInputString(obs_input, ss.str().c_str(), obs_string, 0) == KEY_FOUND) {
			add_observable(obs_string);
		}
		else found = false;
		i++;
	}

	cleanInputFile(obs_input);
	delete obs_input;
}

template<typename number>
ObservableOutput<number>::~ObservableOutput() {
	if(_output_stream.is_open()) _output_stream.close();

	typename vector<BaseObservable<number> *>::iterator it;
	for(it = _obss.begin(); it != _obss.end(); it++) delete *it;
}

template<typename number>
void ObservableOutput<number>::_open_output() {
	if(!strncmp (_output_name, "stderr", 512)) _output = &std::cerr;
	else if (!strncmp (_output_name, "stdout", 512)) _output = &std::cout;
	else {
		if(!_only_last) {
			if(_append)	{
				if (_is_binary) _output_stream.open(_output_name, ios::binary | ios_base::app);
				else _output_stream.open(_output_name, ios::binary | ios_base::app);
			}
			else {
				if (_is_binary) _output_stream.open(_output_name, ios::binary);
				else _output_stream.open(_output_name);
			}
		}

		_output = &_output_stream;
	}
}

template<typename number>
void ObservableOutput<number>::init(ConfigInfo<number> &config_info) {
	_open_output();

	typename vector<BaseObservable<number> *>::iterator it;
	for(it = _obss.begin(); it != _obss.end(); it++) (*it)->init(config_info);
}

template<typename number>
void ObservableOutput<number>::add_observable(std::string obs_string) {
	input_file *obs_inp = Utils::get_input_file_from_string(obs_string);

	BaseObservable<number> *new_obs = ObservableFactory::make_observable<number>(*obs_inp, _sim_inp);
	_obss.push_back(new_obs);

	cleanInputFile(obs_inp);
	delete obs_inp;
}

template<typename number>
void ObservableOutput<number>::change_output_file(string new_filename) {
	sprintf(_output_name, "%s%s", _prefix.c_str(), new_filename.c_str());
	if(_output_stream.is_open()) _output_stream.close();
	_open_output();
}

template<typename number>
bool ObservableOutput<number>::is_ready(llint step) {
	if(_print_every < 1) return false;
	if(_stop_at > -1 && step > _stop_at) return false;
	return ((step >= _start_from) && (step % _print_every == 0));
}

template<typename number>
void ObservableOutput<number>::print_output(llint step) {
	stringstream ss;
	typename vector<BaseObservable<number> *>::iterator it;
	for(it = _obss.begin(); it != _obss.end(); it++) {
		if(it != _obss.begin()) ss << " ";
		ss << (*it)->get_output_string(step);
	}

	if(_only_last) _output_stream.open(_output_name);
	*_output << ss.str() << endl;
	if(_only_last) _output_stream.close();
}

template class ObservableOutput<float>;
template class ObservableOutput<double>;
