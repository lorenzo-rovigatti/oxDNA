/**
 * @file    RepulsionPlaneMoving.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef REPULSIONPLANEMOVING_H_
#define REPULSIONPLANEMOVING_H_

#include "BaseForce.h"

/// plane moving with particle
template<typename number>
class RepulsionPlaneMoving : public BaseForce<number> {
private:
	number * box_side_ptr;
	std::string _particles_string;
	int _ref_id;

public:
	BaseParticle<number> * _p_ptr;

	RepulsionPlaneMoving ();
	virtual ~RepulsionPlaneMoving() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);


	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif // REPULSIONPLANEMOVING_H_
