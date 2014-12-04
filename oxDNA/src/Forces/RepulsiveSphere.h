/**
 * @file    RepulsiveSphere.h
 * @date    28/nov/2014
 * @author  Flavio 
 *
 */

#ifndef REPULSIVESPHERE_H_
#define REPULSIVESPHERE_H_

#include "BaseForce.h"

/// pos0 * (x,y,z) + position = 0 is the definition of the plane.
/// The pos0 vector is pointing to the halfplane where the
/// repulsion is not acting!
template<typename number>
class RepulsiveSphere : public BaseForce<number> {
private:
	int _particle;
	LR_vector<number> _center;
	number _r0, _rate;
	number * _box_side_ptr;

public:
	RepulsiveSphere();
	virtual ~RepulsiveSphere() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif // REPULSIVESPHERE
