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

class YasutakaAnalysis: public Configuration {
protected:
	enum {
		BONDS
	};

	int _mode;
	number _threshold;

	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle *p);
	virtual std::string _configuration(llint step);

public:
	YasutakaAnalysis();
	virtual ~YasutakaAnalysis();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init();
};

extern "C" BaseObservable *make_YasutakaAnalysis() {
	return new YasutakaAnalysis();
}

#endif /* YASUTAKAANALYSIS_H_ */
