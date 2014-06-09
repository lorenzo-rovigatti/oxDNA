/**
 * @file    MDBackend.h
 * @date    23/nov/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MDBACKEND_H_
#define MDBACKEND_H_

#include "SimBackend.h"

using namespace std;

/**
 * @brief Abstract class for MD simulations.
 *
 * This class sets up a basic MD simulation. It does not do any sensible computation but
 * takes care of the most basic input/output operations associated with MD simulations.
 */
template<typename number>
class MDBackend: public SimBackend<number> {
protected:
	number _dt;
	bool _refresh_velocities;

	void _generate_vel();

public:
	MDBackend();
	virtual ~MDBackend();

	void get_settings(input_file &inp);
	void init(char conf_filename[256]);
};

#endif /* MDBACKEND_H_ */
