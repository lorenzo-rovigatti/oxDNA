/*
 * oxDNAException.cpp
 *
 *  Created on: 04/mar/2013
 *      Author: lorenzo
 */

#include "oxDNAException.h"
#include "Utils.h"

oxDNAException::oxDNAException(std::string ss, ...) {
	va_list ap;
	va_start(ap, ss);
	_error = Utils::sformat_ap(ss, ap);
	va_end(ap);
}

const char* oxDNAException::what() const throw() {
	return _error.c_str();
}
