/**
 * @file    ExternalForce.h
 * @date    18/lug/2011
 * @author  lorenzo
 *
 *
 */

#ifndef EXTERNALFORCE_H_
#define EXTERNALFORCE_H_

#include <string>

#include "../defs.h"

/**
 * @brief Interface for external forces.
 *
 * This class contains public members with names starting with underscores.
 * We have to change scope policy due to the fact that GPU classes
 * require access to these members.
 */
template<typename number>
class ExternalForce {
private:
	/**
	 * @brief A name of the group this force belongs to.
	 *
	 * Different forces can be grouped under the same name. This can be used to
	 * act separately on different forces. For example, it can be used together
	 * with the observable ForceEnergy to print only the energy due to specific
	 * groups of forces.
	 */
	std::string _group_name;

public:
	/// we need these members to be public because
	/// we need access in order to copy these numbers
	/// to the GPU memory
	number _rate;
	number _F0;
	LR_vector<number> _direction;
	LR_vector<number> _pos0;
	number _stiff;

	ExternalForce(number base, number rate, LR_vector<number> dir);
	virtual ~ExternalForce();

	virtual void set_group_name(std::string &name) { _group_name = name; }
	virtual std::string get_group_name() { return _group_name; }

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

template<typename number> class BaseParticle;

// mutual trap
template<typename number>
class MutualTrap : public ExternalForce<number> {
public:
	BaseParticle<number> * _p_ptr;
	number _r0;
	bool PBC;
	number * box_side_ptr;

	//MutualTrap (number stiff, number r0, BaseParticle<number> *p);
	//MutualTrap (number stiff, number r0, BaseParticle<number> *p, number * box_side_ptr);
	MutualTrap (number stiff, number r0, BaseParticle<number> *p, number * box_side_ptr, bool use_PBC);
	virtual ~MutualTrap() {}

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
protected:
	LR_vector<number> _distance(LR_vector<number> u, LR_vector<number> v);
};

template<typename number>
class MovingTrap : public ExternalForce<number> {
public:
	MovingTrap (number stiff, LR_vector<number> pos0, number rate, LR_vector<number> dir);
	virtual ~MovingTrap() {}

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

template<typename number>
class LowdimMovingTrap : public ExternalForce<number> {
public:
	LowdimMovingTrap (number stiff, LR_vector<number> pos0, number rate, LR_vector<number> dir, bool vis_X, bool vis_Y, bool vis_Z);
	virtual ~LowdimMovingTrap() {}

	bool _visX;
	bool _visY;
	bool _visZ;

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};


/// pos0 * (x,y,z) + position = 0 is the definition of the plane.
/// The pos0 vector is pointing to the halfplane where the
/// repulsion is not acting!
template<typename number>
class RepulsionPlane : public ExternalForce<number> {
public:
	number _position;

	RepulsionPlane (number stiff, number position, LR_vector<number> dir);
	virtual ~RepulsionPlane() {}

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};


/// plane moving with particle
template<typename number>
class RepulsionPlaneMoving : public ExternalForce<number> {
private:
	number * box_side_ptr;
public:
	BaseParticle<number> * _p_ptr;

	RepulsionPlaneMoving (number stiff, BaseParticle<number> * p, LR_vector<number> dir, number * ptr);
	virtual ~RepulsionPlaneMoving() {}

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

/// Constant Torque implemented as moving traps
template<typename number>
class ConstantRateTorque : public ExternalForce<number> {
public:
	LR_vector<number> _center, _pos0, _axis, _mask;
	number _base, _rate;
	ConstantRateTorque (number, number, number, LR_vector<number>, LR_vector<number>, LR_vector<number>, LR_vector<number>);
	virtual ~ConstantRateTorque() {}

	// it returns a force, so it is easier to set in BaseParticle.cpp
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif /* EXTERNALFORCE_H_ */
