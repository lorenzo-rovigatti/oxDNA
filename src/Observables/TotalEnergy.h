/*
 * TotalEnergy.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef TOTALENERGY_H_
#define TOTALENERGY_H_

#include "BaseObservable.h"
#include "PotentialEnergy.h"
#include "KineticEnergy.h"

/**
 * @brief Outputs the total (kinetic + potential) energy of the system.
 */

class TotalEnergy: public BaseObservable {
private:
	PotentialEnergy _pot_energy;
	KineticEnergy _kin_energy;
public:
	TotalEnergy();
	virtual ~TotalEnergy();

	void get_settings(input_file &my_inp, input_file &sim_inp);

	virtual void init();
	virtual std::string get_output_string(llint curr_step);

	number get_U(llint curr_step);
	number get_K(llint curr_step);
};

#endif /* TOTALENERGY_H_ */
