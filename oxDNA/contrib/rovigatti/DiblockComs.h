/*
 * Gyradius.h
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#ifndef DIBLOCKCOMS_H_
#define DIBLOCKCOMS_H_

#include "BaseObservable.h"

/**
 * @brief Computes the A-A, A-B, B-B and intramolecular g(r) of a two-chain diblock copolymer system.
 *
 */
template<typename number>
class DiblockComs : public BaseObservable<number> {
protected:
	bool _only_intra;

public:
	DiblockComs();
	virtual ~DiblockComs();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
	virtual std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new DiblockComs<float>(); }
extern "C" BaseObservable<double> *make_double() { return new DiblockComs<double>(); }

#endif /* DIBLOCKCOMS_H_ */
