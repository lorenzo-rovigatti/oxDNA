/**
 * @file    GenericConstantForce.h
 * @date    22/jun/2017
 * @author  Lorenzo
 *
 */

#ifndef GENERICCENTRALFORCE_H_
#define GENERICCENTRALFORCE_H_

#include "BaseForce.h"

template<typename number>
class LookupTable;

/**
 * @brief A customizable radial force.
 *
 @verbatim
 particle = <int> (comma-separated list of indices of particles to apply the force to. -1 applies it to all particles. Entries separated by a dash "-" get expanded in a list of all the particles on a same strand comprised between the two indices. E.g., particle= 1,2,5-7 applies the force to 1,2,5,6,7 if 5 and 7 are on the same strand.)
 center = <float>,<float>,<float> (the centre from which the force originates from)
 @endverbatim
 */
template<typename number>
class GenericCentralForce: public BaseForce<number> {
private:
	std::string _particles_string;

	enum {
		GRAVITY, INTERPOLATED
	};

	std::map<std::string, int> _supported_types;
	int _type;
	number _E_shift;

	std::string _table_filename;
	int _table_N;
	LookupTable<number> _table;

public:
	GenericCentralForce();
	virtual ~GenericCentralForce();

	void get_settings(input_file &);
	void init(BaseParticle<number> **, int, BaseBox<number> *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);

	LR_vector<number> center;
	number inner_cut_off, inner_cut_off_sqr;
	number outer_cut_off, outer_cut_off_sqr;
};

template<typename number>
class LookupTable {
public:
	LookupTable() {
		_N = -1;
	}

	void init(std::string filename, int N);
	number query_function(number x);
	number query_derivative(number x);
private:
	number _linear_interpolation(number x, std::vector<number> &x_data, std::vector<number> &fx_data);

	int _N;
	number _delta, _inv_sqr_delta, _xlow, _xupp;
	std::vector<number> _A, _B, _C, _D;
};

#endif /* GENERICCENTRALFORCE_H_ */
