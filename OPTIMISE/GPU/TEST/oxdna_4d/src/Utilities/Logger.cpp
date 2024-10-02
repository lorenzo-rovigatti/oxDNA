/*
 * Logger.cpp
 *
 *  Created on: 04/ott/2013
 *      Author: lorenzo
 */

#include "oxDNAException.h"
#include "Logger.h"

std::shared_ptr<Logger> Logger::_logger = nullptr;

Logger::Logger() {
	_log_open = false;
	_log_stream = stderr;
	_allow_log = true;
	_debug = false;

	sprintf(_log_level_strings[LOG_INFO], "%s", "INFO");
	sprintf(_log_level_strings[LOG_WARNING], "%s", "WARNING");
	sprintf(_log_level_strings[LOG_DEBUG], "%s", "DEBUG");
	sprintf(_log_level_strings[LOG_ERROR], "%s", "ERROR");
}

Logger::~Logger() {
	if(_log_open) {
		fclose(_log_stream);
	}
}

void Logger::_log(int log_level, const char *format, va_list &ap) {
	if(!_allow_log)	{
		return;
	}

	if(log_level != LOG_NOTHING) {
		fprintf(_log_stream, "%s: ", _log_level_strings[log_level]);
	}
	vfprintf(_log_stream, format, ap);
	va_end(ap);

	fprintf(_log_stream, "\n");
	fflush(_log_stream);
}

void Logger::_set_stream(const char *filename) {
	FILE *buff = fopen(filename, "w");
	if(buff == nullptr) {
		throw oxDNAException("Log file '%s' is not writable", filename);
	}

	_log_stream = buff;
	_log_open = true;
}

void Logger::log(int log_level, const char *format, ...) {
	va_list ap;
	va_start(ap, format);
	_log(log_level, format, ap);
}

void Logger::debug(const char *format, ...) {
	if(_debug == true) {
		va_list ap;
		va_start(ap, format);
		_log(LOG_DEBUG, format, ap);
	}
}

void Logger::get_settings(input_file &inp) {
	char filename[256];

	if(getInputString(&inp, "log_file", filename, 0) == KEY_FOUND) _set_stream(filename);
	getInputBool(&inp, "debug", &_debug, 0);
}

void Logger::init() {
	if(_logger != nullptr) {
		throw oxDNAException("The logger has been already initialised");
	}

	_logger = std::shared_ptr<Logger>(new Logger());
}

std::shared_ptr<Logger> Logger::instance() {
	if(_logger == nullptr) {
		throw oxDNAException("Trying to access an uninitialised logger");
	}

	return _logger;
}
