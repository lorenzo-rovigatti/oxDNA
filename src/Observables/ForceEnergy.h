/*
 * ForceEnergy.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef FORCEENERGY_H_
#define FORCEENERGY_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the energy of the system due to external forces acting upon particles.
 */
template<typename number>
class ForceEnergy : public BaseObservable<number> {
protected:
	std::string _group_name;
public:
	ForceEnergy();
	virtual ~ForceEnergy();

	virtual void get_settings (input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
};

#endif /* FORCEENERGY_H_ */
