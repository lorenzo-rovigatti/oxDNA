/*
 * Utils.cpp
 *
 *  Created on: 04/set/2010
 *      Author: lorenzo
 */

#include "Utils.h"
#include "oxDNAException.h"

#include <sstream>

#include <errno.h>
extern int errno;

using std::string;

Utils::Utils() {

}

Utils::~Utils() {

}

int Utils::decode_base(char c) {
	c = toupper(c);
	switch (c) {
	case 'A':
		return N_A;
	case 'C':
		return N_C;
	case 'G':
		return N_G;
	case 'T':
		return N_T;
	case 'U':
		return N_T;
	case 'D':
		return N_DUMMY;
	default:
		return P_INVALID;
	}
}

char Utils::encode_base(int b) {
	switch (b) {
	case N_A:
		return 'A';
	case N_C:
		return 'C';
	case N_G:
		return 'G';
	case N_T:
		return 'T';
	case N_DUMMY:
		return 'D';
	default:
		return 'X';
	}
}

/*
std::vector<int> Utils::get_list_from_string (std::string arg) {
	vector<int> ret;

	vector<string> sections = Utils::split (arg.c_str(), ',');

	if (sections.size() == 0) throw oxDNAException ("Invalid string %s");

	for (int i = 0; i < sections.size(); i ++) {
		vector<string> subsections = Utils::split(sections[i].c_str(), '-');
		for (int j = 0; j < subsections.size(); j ++) {
			ret.push_back (atoi (subsections[j].c_str));
		}
	}

	return ret;
}
*/

std::vector<std::string> Utils::split(const string &s, char delim) {
	string s_copy(s);
	if (delim == ' ')
		Utils::trim(s_copy);
	std::vector<string> elems;
	std::stringstream ss(s_copy);
	string item;

	while (getline(ss, item, delim)) {
		if (delim == ' ') {
			Utils::trim(item);
			if (item.length() > 0)
				elems.push_back(item);
		} else
			elems.push_back(item);
	}

	return elems;
}

std::string Utils::sformat(const std::string &fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	std::string str = Utils::sformat_ap(fmt, ap);
	va_end(ap);
	return str;
}

// c++ wrapper around sprintf (WTF?)
std::string Utils::sformat_ap(const std::string &fmt, va_list &ap) {
	int size = 500;
	std::string str;
	while (1) {
		str.resize(size);
		int n = vsnprintf((char *) str.c_str(), size, fmt.c_str(), ap);
		if (n > -1 && n < size) {
			str.resize(n);
			return str;
		}
		if (n > -1)
			size = n + 1;
		else
			size *= 2;
	}
	return str;
}

input_file *Utils::get_input_file_from_string(const std::string &inp) {
	std::string real_inp(inp);

	if (inp[0] == '{') {
		int sum = 0;
		for (unsigned int i = 0; i < inp.size(); i++) {
			if (inp[i] == '{')
				sum += 1;
			else if (inp[i] == '}')
				sum -= 1;
			if (sum == 0) {
				real_inp = inp.substr(1, i - 1);
				break;
			}
		}
	}

	errno = 0;

	FILE *temp = NULL;
	const int max_tries = 100;
	for(int i = 0; i < max_tries && temp == NULL; i++) temp = tmpfile();
	if(temp == NULL) throw oxDNAException("Failed to create a temporary file, exiting");
	int check = fprintf(temp, "%s", real_inp.c_str());
	if (check != (int)real_inp.size()) throw oxDNAException ("Failed to write to temporary file...; maybe /tmp has no space left? Aborting");

	rewind(temp);
	if (errno == ENOSPC) throw oxDNAException ("Failed to write to temporary file. No space left on device. maybe /tmp has no space left? Aborting");

	input_file *ret = new input_file;
	loadInput(ret, temp);

	fclose(temp);

	return ret;
}

template<typename number>
number Utils::get_temperature(char * raw_T) {
	char deg;
	double tmp_T;
	number T;
	int res = sscanf(raw_T, "%lf %c", &tmp_T, &deg);
	if (res == 2) {
		deg = tolower(deg);
		switch (deg) {
		case 'c':
			T = (number) ((tmp_T + 273.15) * 0.1 / 300.); // convert to kelvin and then to simulation units
			break;
		case 'k':
			T = (number) (tmp_T * 0.1 / 300.); // convert to simulation units
			break;
		default:
			throw oxDNAException("Unrecognizable temperature '%s'", raw_T);
			/* no break */
		}
	} else
		T = (number) tmp_T;

	return T;
}

std::string Utils::bytes_to_human (llint bytes) {
	llint base = 1024;
	int ctr = 0;
	while (bytes / base > 0 && ctr < 4) {
		base *= 1024;
		ctr ++;
	}
	base /= 1024;
	std::string ret = Utils::sformat ("%7.3lf ", bytes / (double) base);
	switch (ctr) {
		case 0:
			ret += std::string(" B");
			break;
		case 1:
			ret += std::string("KB");
			break;
		case 2:
			ret += std::string("MB");
			break;
		case 3:
			ret += std::string("GB");
			break;
		default:
			throw oxDNAException ("Should never get here... (ctr = %d) in %s:%d\n", ctr, __FILE__, __LINE__);
	}
	return ret;
}

/**
 * @brief fills the memory pointed to by seedptr with the current state of
 * the random number generator. Does not handle the memory: it assumes that
 * it can overwrite the first 48 bit of the memory.
 *
 * @param seedptr the memory address to store the 48 bits of the seed into.
 */
void Utils::get_seed(unsigned short * seedptr) {
	unsigned short seme[3] = { 0, 0, 0 };
	unsigned short * tmpptr;
	tmpptr = seed48(seme);
	memcpy(seedptr, tmpptr, 3 * sizeof(unsigned short));
	seed48(seedptr);
	seed48(seedptr);
}

template float Utils::gaussian<float>();
template double Utils::gaussian<double>();

template float Utils::get_temperature<float>(char *);
template double Utils::get_temperature<double>(char *);
