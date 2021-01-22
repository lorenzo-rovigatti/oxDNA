/*
 * PatchyToMgl.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef PATCHYTOMGL_H_
#define PATCHYTOMGL_H_

#include "Observables/Configurations/Configuration.h"

/**
 * @brief Prints an mgl configuration for patchy systems.
 */

class PatchyToMgl: public Configuration {
protected:
	bool _first_neighbours;
	bool _second_neighbours;
	number _patch_size;
	number _threshold;

	std::set<BaseParticle *> _printed;

	std::string _mgl_patchy_line(BaseParticle *p, const char *color, bool print_p, const char *p_color = "");
	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle *p);
	virtual std::string _configuration(llint step);

public:
	PatchyToMgl();
	virtual ~PatchyToMgl();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init();
};

extern "C" BaseObservable *make_PatchyToMgl() {
	return new PatchyToMgl();
}

#endif /* PATCHYTOMGL_H_ */
