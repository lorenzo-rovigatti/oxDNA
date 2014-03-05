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
template<typename number>
class TotalEnergy: public BaseObservable<number> {
private:
	static number _U;
	static number _K;
	static llint _saved_on;

	PotentialEnergy<number> _pot_energy;
	KineticEnergy<number> _kin_energy;
public:
	TotalEnergy();
	virtual ~TotalEnergy();

	virtual void init(ConfigInfo<number> &config_info);
	virtual std::string get_output_string(llint curr_step);

	number get_U(llint curr_step);
	number get_K(llint curr_step);
};

#endif /* TOTALENERGY_H_ */
