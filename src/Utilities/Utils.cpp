/*
 * Utils.cpp
 *
 *  Created on: 04/set/2010
 *      Author: lorenzo
 */

#include "Utils.h"

#include "../Particles/TEPParticle.h"
#include "../Particles/DNANucleotide.h"
#include "../Particles/RNANucleotide.h"

#include <sstream>

using std::string;

namespace Utils {

int decode_base(char c) {
	c = toupper(c);
	switch(c) {
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

char encode_base(int b) {
	switch(b) {
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

std::vector<int> btypes_from_sequence(const std::string &sequence) {
	std::vector<int> result;
	bool open_parenthesis = false;
	std::string parenthesis_token;

	for(char c : sequence) {
		if(c == '(') {
			open_parenthesis = true;
			parenthesis_token.clear();
			continue;
		}
		else if(c == ')') {
			if(!open_parenthesis) {
				throw oxDNAException("unbalanced parenthesis");
			}
			open_parenthesis = false;
			int btype = std::atoi(parenthesis_token.c_str());
			result.push_back(btype);
		}
		else {
			if(open_parenthesis) {
				parenthesis_token.push_back(c);
			}
			else {
				int btype = Utils::decode_base(c);
				if(btype == P_INVALID) {
					throw oxDNAException("invalid base %c", c);
				}
				result.push_back(btype);
			}
		}
	}

	if(open_parenthesis) {
		throw oxDNAException("missing closing parenthesis");
	}

	return result;
}

std::vector<std::string> split(const string &s, char delim) {
	string s_copy(s);
	if(delim == ' ') {
		trim(s_copy);
	}
	std::vector<string> elems;
	std::stringstream ss(s_copy);
	string item;

	while(getline(ss, item, delim)) {
		if(delim == ' ') {
			trim(item);
			if(item.length() > 0) {
				elems.push_back(item);
			}
		}
		else
			elems.push_back(item);
	}

	return elems;
}

std::vector<number> split_to_numbers(const std::string &str, const std::string &delims) {
	std::vector<number> output;
	output.reserve(15);

	const char *ptr = str.c_str();
	while(ptr) {
		auto base = ptr;
		ptr = std::strpbrk(ptr, delims.c_str());
		if(ptr) {
			// this check makes sure that no empty strings are added to the output
			if(ptr - base) {
				output.emplace_back(lexical_cast(std::string(base, ptr - base)));
			}
			ptr++;
		}
		else {
			std::string remainder(base);
			if(remainder.size() > 0) {
				output.emplace_back(lexical_cast(remainder));
			}
		}
	}

	return output;
}

std::string sformat(std::string fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	std::string str = sformat_ap(fmt, ap);
	va_end(ap);
	return str;
}

// c++ wrapper around sprintf (WTF?)
std::string sformat_ap(const std::string &fmt, va_list &ap) {
	int size = 500;
	std::string str;
	while(1) {
		str.resize(size);
		va_list ap_copy;
		va_copy(ap_copy, ap);
		int n = vsnprintf((char*) str.c_str(), size, fmt.c_str(), ap_copy);
		va_end(ap_copy);
		if(n > -1 && n < size) {
			str.resize(n);
			return str;
		}
		if(n > -1)
			size = n + 1;
		else
			size *= 2;
	}
	return str;
}

void orthonormalize_matrix(LR_matrix &m) {
	number v1_norm2 = m.v1 * m.v1;
	number v2_v1 = m.v2 * m.v1;

	m.v2 -= (v2_v1 / v1_norm2) * m.v1;

	number v3_v1 = m.v3 * m.v1;
	number v3_v2 = m.v3 * m.v2;
	number v2_norm2 = m.v2 * m.v2;

	m.v3 -= (v3_v1 / v1_norm2) * m.v1 + (v3_v2 / v2_norm2) * m.v2;

	m.v1.normalize();
	m.v2.normalize();
	m.v3.normalize();
}

input_file* get_input_file_from_string(const std::string &inp) {
	std::string real_inp(inp);

	if(inp[0] == '{') {
		int sum = 0;
		for(unsigned int i = 0; i < inp.size(); i++) {
			if(inp[i] == '{')
				sum += 1;
			else if(inp[i] == '}')
				sum -= 1;
			if(sum == 0) {
				real_inp = inp.substr(1, i - 1);
				break;
			}
		}
	}

	input_file *ret = new input_file();
	ret->init_from_string(real_inp);

	return ret;
}

number get_temperature(std::string raw_T) {
	static std::set<std::string> converted_temperatures;

	bool print_output = false;
	if(converted_temperatures.find(raw_T) == converted_temperatures.end()) {
		converted_temperatures.insert(raw_T);
		print_output = true;
	}

	char deg;
	double tmp_T;
	number T;
	int res = sscanf(raw_T.c_str(), "%lf %c", &tmp_T, &deg);
	if(res == 2) {
		deg = tolower(deg);
		switch(deg) {
		case 'c':
			T = (number) ((tmp_T + 273.15) * 0.1 / 300.); // convert to kelvin and then to simulation units
			if(print_output) {
				OX_LOG(Logger::LOG_INFO, "Converting temperature from Celsius (%lf C°) to simulation units (%lf)", tmp_T, T);
			}
			break;
		case 'k':
			T = (number) (tmp_T * 0.1 / 300.); // convert to simulation units
			if(print_output) {
				OX_LOG(Logger::LOG_INFO, "Converting temperature from Kelvin (%lf K) to simulation units (%lf)", tmp_T, T);
			}
			break;
		default:
			throw oxDNAException("Unrecognizable temperature '%s'", raw_T.c_str());
			/* no break */
		}
	}
	else {
		T = (number) tmp_T;
	}

	return T;
}

std::string bytes_to_human(llint bytes) {
	llint base = 1024;
	int ctr = 0;
	while(bytes / base > 0 && ctr < 4) {
		base *= 1024;
		ctr++;
	}
	base /= 1024;
	std::string ret = sformat("%7.3lf ", bytes / (double) base);
	switch(ctr) {
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
		throw oxDNAException("Should never get here... (ctr = %d) in %s:%d\n", ctr, __FILE__, __LINE__);
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
void get_seed(unsigned short *seedptr) {
	unsigned short seme[3] = { 0, 0, 0 };
	unsigned short *tmpptr;
	tmpptr = seed48(seme);
	memcpy(seedptr, tmpptr, 3 * sizeof(unsigned short));
	seed48(seedptr);
	seed48(seedptr);
}

void assert_is_valid_particle(int index, int N, std::string identifier) {
	if(index >= N || index < -1) {
		throw oxDNAException("Trying to add a %s on non-existent particle %d. Aborting", identifier.c_str(), index);
	}
}

bool is_integer(std::string s) {
	return s.find_first_not_of("0123456789") == string::npos;

}

std::vector<int> get_particles_from_string(std::vector<BaseParticle*> &particles, std::string particle_string, std::string identifier) {
	// first remove all the spaces from the string, so that the parsing goes well.
	particle_string.erase(remove_if(particle_string.begin(), particle_string.end(), static_cast<int (*)(int)>(isspace)), particle_string.end());

	std::vector<std::string> temp = split(particle_string.c_str(), ',');
	std::vector<int> particles_index;

	// try to understand whether we are dealing with a strand-based system or not
	bool has_strands = false;
	if(dynamic_cast<DNANucleotide*>(particles[0]) != NULL) {
		has_strands = true;
	}
	else if(dynamic_cast<RNANucleotide*>(particles[0]) != NULL) {
		has_strands = true;
	}
	else if(dynamic_cast<TEPParticle*>(particles[0]) != NULL) {
		has_strands = true;
	}

	for(std::vector<std::string>::size_type i = 0; i < temp.size(); i++) {
		bool found_dash = temp[i].find('-') != std::string::npos;
		// if the string contains a dash, then it has to be interpreted as a list of particles
		// unless it's a negative number

		if(found_dash && '-' != temp[i].c_str()[0]) {
			// get the two indices p0 and p1 and check they make sense
			std::vector<std::string> p0_p1_index = split(temp[i].c_str(), '-');

			int p[2] = { 0 };
			// check whether the p0 and p1 keys can be understood, and set them
			for(int ii = 0; ii < 2; ii++) {
				if(is_integer(p0_p1_index[ii])) {
					p[ii] = atoi(p0_p1_index[ii].c_str());
					assert_is_valid_particle(p[ii], particles.size(), identifier);
				}
				else {
					if(p0_p1_index[ii] == "last") {
						p[ii] = particles.size() - 1;
					}
					else {
						throw oxDNAException("In %s I couldn't interpret particle identifier \"%s\" used as a boundary particle.", identifier.c_str(), p0_p1_index[ii].c_str());
					}
				}
			}

			// the behaviour changes whether the particles are arranged on strands (DNA, RNA, TEP) or not (everything else)

			// add all the particles between p0 and p1 (extremes included)
			if(has_strands) {
				int j = p[0];
				bool found_p1 = false;
				do {
					particles_index.push_back(j);
					if(j == p[1]) {
						found_p1 = true;
					}
					if(particles[j]->n5 == P_VIRTUAL)
						break;
					j = particles[j]->n5->index;

				} while(j != p[0] && !found_p1);
				// check that it hasn't got to either the end of the strand or back to p1
				if(!found_p1) {
					throw oxDNAException("In %s I couldn't get from particle %d to particle %d.", identifier.c_str(), p[0], p[1]);
				}
			}
			else {
				if(p[0] >= p[1])
					throw oxDNAException("%s: the two indexes in a particle range (here %d and %d) should be sorted (the first one should be smaller than the second one).", identifier.c_str(), p[0], p[1]);
				for(int p_idx = p[0]; p_idx <= p[1]; p_idx++) {
					particles_index.push_back(p_idx);
				}
			}

		}
		else if(temp[i] == "last") {
			particles_index.push_back(particles.size() - 1);
		}
		else if(temp[i] == "all") {
			particles_index.push_back(-1);
		}
		// add it to the vector, and make sure that the identifier is not an unidentified string
		else {
			if(temp[i] != "-1" && !is_integer(temp[i])) {
				throw oxDNAException("In %s I couldn't interpret particle identifier \"%s\".", identifier.c_str(), temp[i].c_str());

			}
			int j = atoi(temp[i].c_str());

			assert_is_valid_particle(j, particles.size(), identifier);
			particles_index.push_back(j);
		}

	}
	// check that if -1 is present then that's the only key - something must be wrong if you
	// specified -1 (all particles) and then some more particles.
	if(std::find(particles_index.begin(), particles_index.end(), -1) != particles_index.end() && particles_index.size() > 1) {
		throw oxDNAException("In %s there is more than one particle identifier, including -1 or \"all\". If either -1 or \"all\" are used as particle identifiers then they have to be the only one, as both translate to \"all the particles\". Dying badly.", identifier.c_str());
	}
	// check that no particle appears twice
	for(std::vector<int>::size_type i = 0; i < particles_index.size(); i++) {
		for(std::vector<int>::size_type j = i + 1; j < particles_index.size(); j++) {
			if(particles_index[i] == particles_index[j]) {
				throw oxDNAException("In %s particle index %d appears twice (both at position %d and at position %d), but each index can only appear once. Dying badly.", identifier.c_str(), particles_index[i], i + 1, j + 1);
			}
		}
	}
	// finally return the vector.
	return particles_index;
}

}

