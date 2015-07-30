/*
 * FSConf.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef FSCONF_H_
#define FSCONF_H_

#include "Configurations/Configuration.h"

/**
 * @brief Prints a configuration in the FS format.
 *
 * @verbatim
[in_box = <bool> (if true all the positions are brought back between -L/2 and L/2. Defaults to false)]
@endverbatim
 */
template<typename number>
class FSConf: public Configuration<number>  {
protected:
	int _N, _N_A, _N_B;
	bool _in_box;
	bool _also_patch;
	bool _print_bonds;
	number _bond_threshold;

	std::vector<std::map<int, int> > _bonds;

	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle<number> *p);

public:
	FSConf();
	virtual ~FSConf();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);

	std::string _configuration(llint step);
};

extern "C" BaseObservable<float> *make_float() { return new FSConf<float>(); }
extern "C" BaseObservable<double> *make_double() { return new FSConf<double>(); }

#endif /* FSCONF_H_ */
