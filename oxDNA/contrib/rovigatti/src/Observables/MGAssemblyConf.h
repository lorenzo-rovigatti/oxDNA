/*
 * MGAssemblyConf.h
 *
 *  Created on: 17/may/2019
 *      Author: lorenzo
 */

#ifndef MGASSEMBLYCONF_H_
#define MGASSEMBLYCONF_H_

#include "Observables/Configurations/Configuration.h"

/**
 * @brief Prints a configuration in the FS format.
 *
 * @verbatim
[in_box = <bool> (if true all the positions are brought back between -L/2 and L/2. Defaults to false)]
@endverbatim
 */
template<typename number>
class MGAssemblyConf: public Configuration<number>  {
protected:
	int _N, _N_A, _N_in_polymers;
	bool _in_box;
	number _bond_threshold;
	bool _with_polymers;

	std::vector<std::map<int, int> > _bonds;

	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle<number> *p);

public:
	MGAssemblyConf();
	virtual ~MGAssemblyConf();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);

	std::string _configuration(llint step);
};

extern "C" BaseObservable<float> *make_float() { return new MGAssemblyConf<float>(); }
extern "C" BaseObservable<double> *make_double() { return new MGAssemblyConf<double>(); }

#endif /* MGASSEMBLYCONF_H_ */
