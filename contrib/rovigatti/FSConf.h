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
 * @brief Prints an mgl configuration for patchy systems.
 */
template<typename number>
class FSConf: public Configuration<number>  {
protected:
	int _N, _N_A, _N_B;

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
