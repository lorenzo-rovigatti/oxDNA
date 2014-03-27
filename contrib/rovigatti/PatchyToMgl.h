/*
 * PatchyToMgl.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef PATCHYTOMGL_H_
#define PATCHYTOMGL_H_

#include "Configurations/Configuration.h"

/**
 * @brief Prints an mgl configuration for patchy systems.
 */
template<typename number>
class PatchyToMgl: public Configuration<number>  {
protected:
	bool _first_neighbours;
	bool _second_neighbours;
	number _patch_size;
	number _threshold;

	std::set<BaseParticle<number> *> _printed;

	std::string _mgl_patchy_line(BaseParticle<number> *p, const char *color, bool print_p, const char *p_color="");
	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle<number> *p);
	virtual std::string _configuration(llint step);

public:
	PatchyToMgl();
	virtual ~PatchyToMgl();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new PatchyToMgl<float>(); }
extern "C" BaseObservable<double> *make_double() { return new PatchyToMgl<double>(); }

#endif /* PATCHYTOMGL_H_ */
