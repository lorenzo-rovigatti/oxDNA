/*
 * Checkpoint.h
 *
 *  Created on: May 16, 2015
 *      Author: flavio
 */

#ifndef CHECKPOINT_H_
#define CHECKPOINT_H_

#include "BaseObservable.h"
#include "Configurations/BinaryConfiguration.h"

/**
 * @brief Outputs the total (kinetic + potential) energy of the system.
 */

class Checkpoint: public BaseObservable {
private:
	BinaryConfiguration _conf;

public:
	Checkpoint();
	virtual ~Checkpoint();

	virtual std::string get_output_string(llint curr_step);
};

#endif /* CHECKPOINT_H_ */
