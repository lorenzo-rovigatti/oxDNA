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
	_bytes_written = 0;
	_output_name = std::string("");
	input_file *obs_input = Utils::get_input_file_from_string(stream_string);

	string out_name;
	getInputString(obs_input, "name", out_name, 1);
	// the prefix should be used only when the output is nor stdout nor stderr
	if(out_name != "stdout" && out_name != "stderr") getInputString(&sim_inp, "output_prefix", _prefix, 0);
	//sprintf(_output_name, "%s%s", _prefix.c_str(), out_name.c_str());
	_output_name = _prefix + out_name;

	getInputLLInt(obs_input, "start_from", &_start_from, 0);
	getInputLLInt(obs_input, "stop_at", &_stop_at, 0);

	int restart_step_counter = 0;
	getInputBoolAsInt(&sim_inp, "restart_step_counter", &restart_step_counter, 0);
	_append = !restart_step_counter;

	_only_last = 0;
	getInputBoolAsInt(obs_input, "only_last", &_only_last, 0);

	_update_name_with_time = false;
	getInputBool(obs_input, "update_name_with_time", &_update_name_with_time, 0);
	if(_update_name_with_time) _base_name = _output_name;

	if(_only_last && _update_name_with_time) throw oxDNAException("only_last and update_name_with_time are mutually exclusive");

	_is_binary = false;
	getInputBool(obs_input, "binary", &_is_binary, 0);

	_linear = true;
	getInputBool(obs_input, "linear", &_linear, 0);
	if(!_linear) {
		getInputInt(obs_input, "log_ppc", &_log_ppc, 1);
		getInputInt(obs_input, "log_n0", &_log_n0, 1);
		getInputNumber(obs_input, "log_fact", &_log_fact, 1);
		_log_next = _log_n0;
		_log_pos_in_cycle = 1;
		_log_n_cycle = 0;
		_log_tot_cycle = (llint) round(_log_n0 * pow(_log_fact, (number)(_log_ppc - 1.)));
	}
	else getInputLLInt(obs_input, "print_every", &_print_every, 1);

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
	if(!strncmp (_output_name.c_str(), "stderr", 512)) _output = &std::cerr;
	else if (!strncmp (_output_name.c_str(), "stdout", 512)) _output = &std::cout;
	else {
		if(!_only_last) {
			if(_append && !_update_name_with_time)	{
				if(_is_binary) _output_stream.open(_output_name.c_str(), ios::binary | ios_base::app);
				else _output_stream.open(_output_name.c_str(), ios::binary | ios_base::app);
			}
			else {
				if (_is_binary) _output_stream.open(_output_name.c_str(), ios::binary);
				else _output_stream.open(_output_name.c_str());
			}
		}

		_output = &_output_stream;
	}

	if(_output->bad() || !_output->good()) throw oxDNAException("Stream %s not writable", _output_name.c_str());
}

template<typename number>
void ObservableOutput<number>::init(ConfigInfo<number> &config_info) {
	if(!_update_name_with_time) _open_output();

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
	_output_name = _prefix + new_filename;
	if(_output_stream.is_open()) _output_stream.close();
	_open_output();
}

template<typename number>
void ObservableOutput<number>::_set_next_log_step() {
	_log_next = (llint)round((_log_n0*pow(_log_fact, _log_pos_in_cycle))) + _log_tot_cycle*_log_n_cycle;
	_log_pos_in_cycle++;
	if(_log_pos_in_cycle == _log_ppc) {
		_log_n_cycle++;
		_log_pos_in_cycle = 0;
	}
}

template<typename number>
bool ObservableOutput<number>::is_ready(llint step) {
	if(_stop_at > -1 && step > _stop_at) return false;

	if(_linear) {
		if(_print_every < 1) return false;
		return ((step >= _start_from) && (step % _print_every == 0));
	}
	else {
		while(_log_next < step) _set_next_log_step();
		return (step == _log_next);
	}
}

template<typename number>
void ObservableOutput<number>::print_output(llint step) {
	stringstream ss;
	typename vector<BaseObservable<number> *>::iterator it;
	for(it = _obss.begin(); it != _obss.end(); it++) {
		if(it != _obss.begin()) ss << " ";
		ss << (*it)->get_output_string(step);
	}

	if(_update_name_with_time) {
		string new_name = Utils::sformat("%s%lld", _base_name.c_str(), step);
		change_output_file(new_name);
	}
	else if(_only_last) _output_stream.open(_output_name.c_str());
	
	ss << endl;
	std::string towrite = ss.str();
	_bytes_written += (llint) towrite.length();
	*_output << towrite;
	_output->flush();
	
	if(_only_last) _output_stream.close();
	if(!_linear) _set_next_log_step();
}

template class ObservableOutput<float>;
template class ObservableOutput<double>;
