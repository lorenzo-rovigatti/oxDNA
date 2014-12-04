/**
 * @file    RepulsionPlane.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef REPULSIONPLANE_H_
#define REPULSIONPLANE_H_

#include "BaseForce.h"

/// pos0 * (x,y,z) + position = 0 is the definition of the plane.
/// The pos0 vector is pointing to the halfplane where the
/// repulsion is not acting!
template<typename number>
class RepulsionPlane : public BaseForce<number> {
private:
	int _particle;

public:
	number _position;

	RepulsionPlane ();
	virtual ~RepulsionPlane() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif // MOVINGTRAP_H
