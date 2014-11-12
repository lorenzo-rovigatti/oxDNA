/*
 * TSPAnalysis.h
 *
 *  Created on: 6 Nov 2014
 *      Author: lorenzo
 */

#ifndef TSPANALYSIS_H_
#define TSPANALYSIS_H_

#include "BaseObservable.h"

#include <vector>

template<typename number>
class TSPAnalysis : public BaseObservable<number>{
protected:
	int _N_stars;
	std::string _topology_filename;

	std::vector<std::vector<BaseParticle<number> *> > _attractive_monomers;
public:
	TSPAnalysis();
	virtual ~TSPAnalysis();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);

	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new TSPAnalysis<float>(); }
extern "C" BaseObservable<double> *make_double() { return new TSPAnalysis<double>(); }

#endif /* TSPANALYSIS_H_ */
