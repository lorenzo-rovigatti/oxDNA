/*
 * oxDNAException.h
 *
 *  Created on: 04/mar/2013
 *      Author: lorenzo
 */

#ifndef OXDNAEXCEPTION_H_
#define OXDNAEXCEPTION_H_

#include <string>
#include <cstdarg>

/**
 * @brief Generic oxDNA exception. It should be handled by either SimManager or AnalysisManager.
 */
class oxDNAException : public std::exception {
private:
	std::string _error;
public:
	oxDNAException(const std::string &ss, ...);
	// I'm not sure why, but the throw() bit is needed to avoid a 'looser throw specifier' error
	virtual ~oxDNAException() throw() {};

	/**
	 * @brief Returns the actual error message.
	 *
	 * @return error message associated to the exception
	 */
	virtual const char* error() const throw();
};

#endif /* OXDNAEXCEPTION_H_ */
