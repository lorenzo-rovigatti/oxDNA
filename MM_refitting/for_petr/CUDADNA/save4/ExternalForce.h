/*
 * ExternalForce.h
 *
 *  Created on: 18/lug/2011
 *      Author: lorenzo
 */

#ifndef EXTERNALFORCE_H_
#define EXTERNALFORCE_H_

#include "defs.h"

template<typename number>
class ExternalForce {
public:
	// we need these members to be public because 
	// we need access in order to copy these numbers
	// to the GPU memory
	number _rate;
	number _F0;
	LR_vector<number> _direction;
	LR_vector<number> _pos0;
	number _stiff;
	
	ExternalForce(number base, number rate, LR_vector<number> dir);
	virtual ~ExternalForce();
	
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos) = 0;
	virtual number potential (llint step, LR_vector<number> &pos) = 0;
};

template<typename number>
class ConstantRateForce : public ExternalForce<number> {
public:
	ConstantRateForce(number base, number rate, LR_vector<number> dir);
	virtual ~ConstantRateForce();
	
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential (llint step, LR_vector<number> &pos);
};

template<typename number> class Particle;

// mutual trap
template<typename number>
class MutualTrap : public ExternalForce<number> {
public:
	Particle<number> * _p_ptr;
	number _r0;

	MutualTrap (number stiff, number r0, Particle<number> *p);
	virtual ~MutualTrap() {}
	
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

template<typename number>
class MovingTrap : public ExternalForce<number> {
public:
	MovingTrap (number stiff, LR_vector<number> pos0, number rate, LR_vector<number> dir);
	virtual ~MovingTrap() {}
	
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

//pos0 * (x,y,z) + position = 0 is the definition of the plane.
//The pos0 vector is pointing to the halfplane where the
// repulsion is not acting!
template<typename number>
class RepulsionPlane : public ExternalForce<number> {
public:
	number _position;

	RepulsionPlane (number stiff, number position, LR_vector<number> dir);
	virtual ~RepulsionPlane() {}

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};


// plane moving with particle
template<typename number>
class RepulsionPlaneMoving : public ExternalForce<number> {
public:
	Particle<number> * _p_ptr;

	RepulsionPlaneMoving (number stiff, Particle<number> * p, LR_vector<number> dir);
	virtual ~RepulsionPlaneMoving() {}

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

// Constant Torque implemented as moving traps
template<typename number>
class ConstantRateTorque : public ExternalForce<number> {
public:
	LR_vector<number> _center, _pos0, _axis, _mask;
	number _base, _rate;
	ConstantRateTorque (number, number, number, LR_vector<number>, LR_vector<number>, LR_vector<number>, LR_vector<number>);
	virtual ~ConstantRateTorque() {}
	
	// it returns a force, so it is easier to set in Particle.cpp
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif /* EXTERNALFORCE_H_ */
