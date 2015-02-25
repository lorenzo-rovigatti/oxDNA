/*
 * Logger.h
 *
 *  Created on: 04/ott/2013
 *      Author: lorenzo
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#define OX_LOG Logger::instance()->log
#define OX_DEBUG Logger::instance()->debug

#include <cstdarg>

#include "parse_input/parse_input.h"

/**
 * @brief Simple logger
 *
 * This logger is a singleton (http://www.yolinux.com/TUTORIALS/C++Singleton.html) used to log
 * infos, warnings and errors either on a log file or to stderr.
 */
class Logger {
private:
	/// Static pointer to the only allowed Logger instance
	static Logger *_logger;
	/// If true, the logger prints also debug informations
	bool _debug;
	/// If false, the logger won't print anything
	bool _allow_log;

	/// Output stream
	FILE *_log_stream;
	/// Output string for each log level
	char _log_level_strings[4][256];
	/// If true, the log stream can be used
	bool _log_open;

	void _log(int log_level, const char *format, va_list &ap);
	void _set_stream(const char *filename);

	/**
	 * @brief Default constructor. It is kept private to enforce the singleton pattern.
	 *
	 * @return
	 */
	Logger();

	/**
	 * @brief Copy constructor. It is kept private to enforce the singleton pattern.
	 *
	 * @param
	 * @return
	 */
	Logger(Logger const&) {};

public:
	enum {
		LOG_INFO = 0,
		LOG_WARNING = 1,
		LOG_DEBUG = 2,
		LOG_ERROR = 3,
		LOG_NOTHING = 4
	} log_levels;

	/**
	 * @brief Read options from the input file
	 *
	 * @param inp reference to an input_file struct containing a dictionary of key/values
	 */
	void get_settings(input_file &inp);

	/**
	 * @brief Initialize the logger. It is static to avoid more than one instantiation
	 */
	static void init();

	/**
	 * @brief Destroy the logger
	 */
	static void clear();

	/**
	 * @brief Variadic method. Does the actual logging
	 *
	 * @param log_level it should be one of the log_levels
	 * @param format printf-like string
	 */
	void log(int log_level, const char *format, ...);

	/**
	 * @brief Variadic method. Used to log debug messages.
	 *
	 * @param format
	 */
	void debug(const char *format, ...);

	/**
	 * @brief Disables logging.
	 */
	void disable_log() { _allow_log = false; }

	/**
	 * @brief Enabled logging.
	 */
	void enable_log() { _allow_log = true; }

	/**
	 * @brief Returns the pointer to the log stream used by this Logger
	 *
	 * @return pointer to the log stream (it can be stdout, stderr or an open file)
	 */
	FILE *get_log_stream() { return _log_stream; }

	/**
	 * @brief Returns the actual logger. Static method to enforce the singleton pattern.
	 *
	 * @return Pointer to an already initialized logger
	 */
	static Logger *instance();

	virtual ~Logger();
};

#endif /* LOGGER_H_ */
