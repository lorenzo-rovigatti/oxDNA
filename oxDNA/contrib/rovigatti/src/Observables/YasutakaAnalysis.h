/*
 * YasutakaAnalysis.h
 *
 *  Created on: 06/oct/2016
 *      Author: lorenzo
 */

#ifndef YASUTAKAANALYSIS_H_
#define YASUTAKAANALYSIS_H_

#include "Observables/Configurations/Configuration.h"

#include <complex>

/**
 * @brief Analyse Yasutaka's systems.
 */
template<typename number>
class YasutakaAnalysis: public Configuration<number>  {
protected:
	enum {
		BONDS
	};

	int _mode;
	number _threshold;

	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle<number> *p);
	virtual std::string _configuration(llint step);

public:
	YasutakaAnalysis();
	virtual ~YasutakaAnalysis();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new YasutakaAnalysis<float>(); }
extern "C" BaseObservable<double> *make_double() { return new YasutakaAnalysis<double>(); }

#endif /* YASUTAKAANALYSIS_H_ */
