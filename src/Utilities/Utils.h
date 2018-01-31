/**
 * @file    Utils.h
 * @date    04/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>
#include <functional>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <string>
#include <cstdarg>
#include <cstdio>
#include <vector>

#include "../defs.h"

template<typename number> class BaseParticle;

/**
 * @brief Utility class. It mostly contains static methods.
 */
class Utils {
public:
	Utils();
	virtual ~Utils();

	static int decode_base(char c);
	static char encode_base(int b);

	template<typename number> static number gaussian();
	template<typename number> static number gamma(number alpha, number beta);
	template<typename number> static number sum(number *v, int N) {
		number res = (number) 0.;
		for(int i = 0; i < N; i++)
			res += v[i];
		return res;
	}

	/**
	 * @brief split a string into tokens, according to the given delimiter
	 *
	 * @param s string to be splitted
	 * @param delim delimiter, defaults to a space
	 * @return a vector of strings containing all the tokens
	 */
	static std::vector<std::string> split(const std::string &s, char delim = ' ');

	// trim from both ends, it works like Python's own trim
	static inline std::string &trim(std::string &s) {
		return ltrim(rtrim(s));
	}

	// trim from start
	static inline std::string &ltrim(std::string &s) {
		s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
		return s;
	}

	// trim from end
	static inline std::string &rtrim(std::string &s) {
		s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
		return s;
	}

	/**
	 * @brief sprintf c++ wrapper (I love c++...).
	 *
	 * @param fmt
	 * @return
	 */
	static std::string sformat(const std::string &fmt, ...);
	/**
	 * @brief vsprintf c++ wrapper.
	 *
	 * This method can be called by other variadic. It is used, for example, by oxDNAException
	 * @param fmt format string
	 * @param ap variadic parameter list initialized by the caller
	 * @return
	 */
	static std::string sformat_ap(const std::string &fmt, va_list &ap);

	/**
	 * @brief Generates a random vector having module 1.
	 *
	 * @return
	 */
	template<typename number> static LR_vector<number> get_random_vector();

	/**
	 * @brief Generates a random vector inside a sphere of given radius.
	 *
	 * @param r sphere radius
	 * @return
	 */
	template<typename number> static LR_vector<number> get_random_vector_in_sphere(number r);

	/**
	 * @brief Applies the Gram-Schmidt orthonormalization to the given matrix.
	 *
	 * @param M the matrix to be orthonormalized
	 */
	template<typename number> static void orthonormalize_matrix(LR_matrix<number> &M);

	/**
	 * @brief Returns a matrix which generates a rotation around a random axis of a random angle, extracted between 0 and max_angle.
	 *
	 * @param max_angle
	 * @return
	 */
	template<typename number> static LR_matrix<number> get_random_rotation_matrix(number max_angle = 2 * M_PI);

	/**
	 * @brief Returns a matrix which generates a rotation around a random axis of the given angle.
	 *
	 * @param angle
	 * @return
	 */
	template<typename number> static LR_matrix<number> get_random_rotation_matrix_from_angle(number angle);

	/**
	 * @brief Creates a temporary file and loads it in an input_file.
	 *
	 * If the string parameter starts with a curly opening bracket, this method will print in the temporary file
	 * only the part of the string which is enclosed by the outer bracket pair.
	 * The calling method is responsible for the calling of cleanInputFile and the deletion of the returned pointer
	 * @param inp string to be written in the temporary file and then loaded in the input_file
	 * @return pointer to the newly loaded input_file
	 */
	static input_file *get_input_file_from_string(const std::string &inp);

	/**
	 * @brief Parses the string passed as the only argument and try to interpret it as a temperature.
	 *
	 * This static method tries to convert raw_T into a temperature. This is mostly used in DNA and
	 * RNA simulations, as the temperature for these can be specified in Celsius or Kelvin degrees.
	 *
	 * For DNA, for example, all the following are equivalent:
	 *
	 * Value 	Simulation Units
	 * 0.1		0.1
	 * 300 K	0.1
	 * 300k		0.1
	 * 26.85c	0.1
	 * 26.85 C	0.1
	 *
	 * @param raw_T c-string containing the text to be parsed
	 * @return
	 */
	template<typename number> static number get_temperature(char *raw_T);

	/**
	 * @brief fills the memory pointed to by seedptr with the current
	 * state of the random number generator.
	 *
	 * This method does not handle the memory: it assumes that it can overwrite the first 48 bit of the
	 * memory.
	 *
	 * @param seedptr the memory address to store the 48 bits of the
	 * seed into.
	 */
	static void get_seed(unsigned short * seedptr);

	/**
	 * @brief Utility function that returns a string from a number converting to megabytes,
	 * kilobytes, etc. Examples: 1010813664 --> "963.987 MB", 783989 --> "765.614 kB" 
	 *
	 * @param arg number to be coverted, assumed in bytes
	 */
	static std::string bytes_to_human(llint arg);

	/**
	 * @brief Utility function that sets the centre of mass velocity of the system to 0.
	 *
	 * @param particles pointer to array of particle pointers
	 * @param N number of particles
	 */
	template<typename number>
	static void stop_com(BaseParticle<number> **particles, int N);

	/**
	 * @brief Utility function that reads a string like "10-16,18" and returns a vector of integers.
	 * @param particles pointer to array of particle pointers
	 * @param N number of particles
	 * @param particles_string string to process
	 * @param identifier the identifier of the calling item (to display to the user in case problems arise).
	 */
	template<typename number>
	static std::vector<int> getParticlesFromString(BaseParticle<number> **particles, int N, std::string particle_string, char const *identifier);

	/**
	 * @brief Utility function that checks if an integer is a valid particle index, or -1.
	 * @param n integer to check.
	 * @param N number of particles
	 * @param identifier the identifier of the calling item (to display to the user in case problems arise).
	 */
	static void assert_is_valid_particle(int n, int N, char const * identifier);
	/**
	 *	@brief Utility function that returns true if a string only contains digits, and false otherwise.
	 *	@param s string to be checked.
	 */
	static bool is_integer(std::string s);

};

template<typename number>
inline LR_vector<number> Utils::get_random_vector() {
	number ransq = 1.;
	number ran1, ran2;

	while(ransq >= 1) {
		ran1 = 1. - 2. * drand48();
		ran2 = 1. - 2. * drand48();
		ransq = ran1 * ran1 + ran2 * ran2;
	}

	number ranh = 2. * sqrt(1. - ransq);
	return LR_vector<number>(ran1 * ranh, ran2 * ranh, 1. - 2. * ransq);
}

template<typename number>
inline LR_vector<number> Utils::get_random_vector_in_sphere(number r) {
	number r2 = SQR(r);
	LR_vector<number> res = LR_vector<number>(r, r, r);

	while(res.norm() > r2) {
		res = LR_vector<number>(2. * r * (drand48() - 0.5), 2. * r * (drand48() - 0.5), 2. * r * (drand48() - 0.5));
	}

	return res;
}

template<typename number>
void Utils::orthonormalize_matrix(LR_matrix<number> &m) {
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

template<typename number>
inline LR_matrix<number> Utils::get_random_rotation_matrix_from_angle(number angle) {
	LR_vector<number> axis = Utils::get_random_vector<number>();

	number t = angle;
	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = 1. - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix<number> R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin, xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin, xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

	return R;
}

template<typename number>
inline LR_matrix<number> Utils::get_random_rotation_matrix(number max_angle) {
	number t = max_angle * (drand48() - 0.5);
	return get_random_rotation_matrix_from_angle(t);
}

template<typename number>
inline number Utils::gaussian() {
	static unsigned int isNextG = 0;
	static number nextG;
	number toRet;
	number u, v, w;

	if(isNextG) {
		isNextG = 0;
		return nextG;
	}

	w = 2.;
	while(w >= 1.0) {
		u = 2. * drand48() - 1.0;
		v = 2. * drand48() - 1.0;
		w = u * u + v * v;
	}

	w = sqrt((-2. * log(w)) / w);
	toRet = u * w;
	nextG = v * w;
	isNextG = 1;

	return toRet;
}

#endif /* UTILS_H_ */
