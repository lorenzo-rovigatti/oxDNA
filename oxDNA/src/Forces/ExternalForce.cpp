/*
 * ExternalForce.cpp
 *
 *  Created on: 18/lug/2011
 *      Author: lorenzo
 */

#include "ExternalForce.h"
#include "../Particles/BaseParticle.h"
#include <cstdio>
#include <cstdlib>

template<typename number>
ExternalForce<number>::ExternalForce(number base, number rate, LR_vector<number> dir) {
	_F0 = base;
	_rate = rate;
	_direction = dir;
}

template<typename number>
ExternalForce<number>::~ExternalForce() {

}

template class ExternalForce<double>;
template class ExternalForce<float>;

template<typename number>
ConstantRateForce<number>::ConstantRateForce(number base, number rate, LR_vector<number> dir) : ExternalForce<number>(base, rate, dir) {
}

template<typename number>
ConstantRateForce<number>::~ConstantRateForce() {

}

template<typename number>
LR_vector<number> ConstantRateForce<number>::value(llint step, LR_vector<number> &pos) {
	number x = (this->_F0 + this->_rate * step) * this->_direction.x;
	number y = (this->_F0 + this->_rate * step) * this->_direction.y;
	number z = (this->_F0 + this->_rate * step) * this->_direction.z;
	return LR_vector<number>(x, y, z);
}

template<typename number>
number ConstantRateForce<number>::potential(llint step, LR_vector<number> &pos) {
	return (number) -(this->_F0 + this->_rate * step) * (pos * this->_direction);
}

template class ConstantRateForce<double>;
template class ConstantRateForce<float>;

// constant force between two particles, to make them come together
template<typename number>
ConstantTrap<number>::ConstantTrap(number stiff, number r0, BaseParticle<number> *p, number *box_ptr, bool use_PBC) : ExternalForce<number>(0, 0, LR_vector<number> (0, 0, 0)) {
	this->_stiff = stiff; // not really a stiffnets, units are energy / distance, not energy / (distance^2)
	this->_r0 = r0;
	this->_p_ptr = p;
	this->box_side_ptr = box_ptr;
	this->PBC = use_PBC;
}

template<typename number>
LR_vector<number> ConstantTrap<number>::_distance(LR_vector<number> u, LR_vector<number> v) {
	if (this->PBC) return v.minimum_image(u, *(this->box_side_ptr));
	else return v - u;
}

template<typename number>
LR_vector<number> ConstantTrap<number>::value (llint step, LR_vector<number> &pos) {
	LR_vector<number> dr = this->_distance(pos, _p_ptr->get_abs_pos(*(this->box_side_ptr))); // other - self
	number sign = copysign (1., (double)(dr.module() - _r0));
	return (this->_stiff * sign) * (dr / dr.module());
}

template <typename number>
number ConstantTrap<number>::potential (llint step, LR_vector<number> &pos) {
	LR_vector<number> dr = this->_distance(pos, _p_ptr->get_abs_pos(*(this->box_side_ptr))); // other - self
	return this->_stiff * (dr.module () - _r0);
}

template class ConstantTrap<double>;
template class ConstantTrap<float>;

// Mutual trap: to make things come together...
template<typename number>
MutualTrap<number>::MutualTrap(number stiff, number r0, BaseParticle<number> * p, number * box_ptr, bool use_PBC) : ExternalForce<number>(0, 0, LR_vector<number> (0, 0, 0)) {
	this->_stiff = stiff;
	this->_r0 = r0;
	this->_p_ptr = p;
	//this->box_side = box_side;
	this->box_side_ptr = box_ptr;
	this->PBC = use_PBC;
}

template<typename number>
LR_vector<number> MutualTrap<number>::_distance(LR_vector<number> u, LR_vector<number> v) {
	if (this->PBC) return v.minimum_image(u, *(this->box_side_ptr));
	else return v - u;
}

template<typename number>
LR_vector<number> MutualTrap<number>::value (llint step, LR_vector<number> &pos) {
	//LR_vector<number> dr = this->_p_ptr->get_abs_pos() - pos; // other - self
	LR_vector<number> dr = this->_distance(pos, _p_ptr->get_abs_pos(*(this->box_side_ptr)));
	return (dr / dr.module()) * (dr.module() - _r0) * this->_stiff;
}

template<typename number>
number MutualTrap<number>::potential (llint step, LR_vector<number> &pos) {
	//LR_vector<number> dr = this->_p_ptr->get_abs_pos() - pos;
	LR_vector<number> dr = this->_distance(pos, _p_ptr->get_abs_pos(*(this->box_side_ptr)));
//printf ("@@@ %lf %lf %lf %lf %i\n", dr.module(), dr.x, dr.y, dr.z, this->_p_ptr->index);
	return pow (dr.module() - _r0, 2) * ((number) 0.5) * this->_stiff;
}
template class MutualTrap<double>;
template class MutualTrap<float>;

// added by Flavio a pene di segugio per Debayan Chakraborty
template<typename number>
MovingTrap<number>::MovingTrap(number stiff, LR_vector<number> pos0, number rate, LR_vector<number> dir) : ExternalForce<number>(0,rate,dir) {
    this->_stiff = stiff;
    this->_pos0 = pos0;
}


template<typename number>
LR_vector<number> MovingTrap<number>::value(llint step, LR_vector<number> &pos) {
    LR_vector<number> postrap;
    number x, y, z;

    postrap.x = this->_pos0.x + (this->_rate * step) * this-> _direction.x;
    postrap.y = this->_pos0.y + (this->_rate * step) * this-> _direction.y;
    postrap.z = this->_pos0.z + (this->_rate * step) * this-> _direction.z;

    x = - this->_stiff * (pos.x - postrap.x);
    y = - this->_stiff * (pos.y - postrap.y);
    z = - this->_stiff * (pos.z - postrap.z);

    return LR_vector<number>(x, y, z);
}

template<typename number>
number MovingTrap<number>::potential (llint step, LR_vector<number> &pos) {
    LR_vector<number> postrap;

    postrap.x = this->_pos0.x + (this->_rate * step) * this-> _direction.x;
    postrap.y = this->_pos0.y + (this->_rate * step) * this-> _direction.y;
    postrap.z = this->_pos0.z + (this->_rate * step) * this-> _direction.z;

    return (number) (0.5 * this->_stiff * (pos - postrap).norm());
}

template class MovingTrap<double>;
template class MovingTrap<float>;

template<typename number>
LowdimMovingTrap<number>::LowdimMovingTrap(number stiff, LR_vector<number> pos0, number rate, LR_vector<number> dir, bool vis_X, bool vis_Y, bool vis_Z) : ExternalForce<number>(0,rate,dir) {
  this->_stiff = stiff;
  this->_pos0 = pos0;
  this->_visX = vis_X;
  this->_visY = vis_Y;
  this->_visZ = vis_Z;
}

template<typename number>
LR_vector<number> LowdimMovingTrap<number>::value(llint step, LR_vector<number> &pos) {
  LR_vector<number> postrap;
  number x, y, z;
  postrap.x = (this->_visX) ? this->_pos0.x + (this->_rate * step) * this-> _direction.x : pos.x;
  postrap.y = (this->_visY) ? this->_pos0.y + (this->_rate * step) * this-> _direction.y : pos.y;
  postrap.z = (this->_visZ) ? this->_pos0.z + (this->_rate * step) * this-> _direction.z : pos.z;

  x = - this->_stiff * (pos.x - postrap.x);
  y = - this->_stiff * (pos.y - postrap.y);
  z = - this->_stiff * (pos.z - postrap.z);

  return LR_vector<number>(x, y, z);
}

template<typename number>
number LowdimMovingTrap<number>::potential (llint step, LR_vector<number> &pos) {
  LR_vector<number> postrap;

  postrap.x = (this->_visX) ? this->_pos0.x + (this->_rate * step) * this-> _direction.x : pos.x;
  postrap.y = (this->_visY) ? this->_pos0.y + (this->_rate * step) * this-> _direction.y : pos.y;
  postrap.z = (this->_visZ) ? this->_pos0.z + (this->_rate * step) * this-> _direction.z : pos.z;

  return (number) (0.5 * this->_stiff * (pos - postrap).norm());
}

template class LowdimMovingTrap<double>;
template class LowdimMovingTrap<float>;

// Repulsion plane moving
template<typename number>
RepulsionPlaneMoving<number>::RepulsionPlaneMoving(number stiff, BaseParticle<number> *p, LR_vector<number> dir, number * ptr) : ExternalForce<number>(0,0,dir) {
	this->_stiff = stiff;
	this->_direction = dir / dir.module(); //this is the normal to the plane
	this->_p_ptr = p;
	this->box_side_ptr = ptr;
}

template<typename number>
LR_vector<number> RepulsionPlaneMoving<number>::value(llint step, LR_vector<number> &pos) {
	number distance = (pos - _p_ptr->get_abs_pos(*(this->box_side_ptr))) * this->_direction;
	if (distance >=  0)
		return LR_vector<number>(0,0,0);
	else
		return -distance * this->_stiff * this->_direction;
}

template<typename number>
number RepulsionPlaneMoving<number>::potential (llint step, LR_vector<number> &pos) {
	number distance = (pos - _p_ptr->get_abs_pos(*(this->box_side_ptr))) * this->_direction;
//	number distance = _p_ptr->pos * this->_direction; //distance from the plane
	if (distance >= 0)
		return 0;
	else
		return (number) (0.5 * this->_stiff * distance * distance);
}

template class RepulsionPlaneMoving<double>;
template class RepulsionPlaneMoving<float>;

// Repulsion plane addendum
template<typename number>
RepulsionPlane<number>::RepulsionPlane(number stiff, number position, LR_vector<number> dir) : ExternalForce<number>(0,0,dir) {
    this->_stiff = stiff;
    this->_direction = dir / dir.module(); //this is the normal to the plane
    this->_position = position;
}

template<typename number>
LR_vector<number> RepulsionPlane<number>::value(llint step, LR_vector<number> &pos) {
    number distance = this->_direction * pos + this->_position;
    if(distance >=  0 )
    {
    	return LR_vector<number>(0,0,0);
    }
    else
    {
    	 return -distance * this->_stiff * this->_direction;
    }
}

template<typename number>
number RepulsionPlane<number>::potential (llint step, LR_vector<number> &pos) {
	number distance = this->_direction * pos + this->_position; //distance from the plane
	if(distance >= 0) {
		return 0;
	}
	else {
		return (number) (0.5 * this->_stiff * distance * distance);
	}
}

template class RepulsionPlane<double>;
template class RepulsionPlane<float>;

// ConstantRateTorque
template<typename number>
ConstantRateTorque<number>::ConstantRateTorque(number stiff, number base, number rate, LR_vector<number> center, LR_vector<number> axis, LR_vector<number> pos0, LR_vector<number> mask) : ExternalForce<number> (0, 0, pos0) {
	this->_stiff = stiff;
	this->_base = base;
	this->_rate = rate;
	this->_center = center;
	this->_axis = axis / axis.module();
	this->_pos0 = pos0;
	this->_mask = mask;
}

template<typename number>
number ConstantRateTorque<number>::potential(llint step, LR_vector<number> &pos) {
	number t = (this->_base + this->_rate * (number) step);

	number sintheta = sin (t);
	number costheta = cos (t);
	number olcos = ((number) 1.) - costheta;

	number xyo = this->_axis.x * this->_axis.y * olcos;
	number xzo = this->_axis.x * this->_axis.z * olcos;
	number yzo = this->_axis.y * this->_axis.z * olcos;
	number xsin = this->_axis.x * sintheta;
	number ysin = this->_axis.y * sintheta;
	number zsin = this->_axis.z * sintheta;

	LR_matrix<number> R(this->_axis.x * this->_axis.x * olcos + costheta,
						xyo - zsin,
						xzo + ysin,

						xyo + zsin,
						this->_axis.y * this->_axis.y * olcos + costheta,
						yzo - xsin,

						xzo - ysin,
						yzo + xsin,
						this->_axis.z * this->_axis.z * olcos + costheta);

	LR_vector<number> postrap = R * (_pos0 - _center) + _center;

	// we "mask" the resulting vector;
	LR_vector<number> dr = pos - postrap;
	dr.x *= _mask.x;
	dr.y *= _mask.y;
	dr.z *= _mask.z;

	return (number) (0.5 * this->_stiff * (dr * dr));
}


template<typename number>
LR_vector<number> ConstantRateTorque<number>::value(llint step, LR_vector<number> &pos) {
	//
	// dobbiamo ruotare pos0 di (base + rate * t) radianti
	// intorno all'asse passante per il centro.
	// e quello ci da la pos della trappola;
	//
	number t = (this->_base + this->_rate * (number) step);

	number sintheta = sin (t);
	number costheta = cos (t);
	number olcos = ((number) 1.) - costheta;

	number xyo = this->_axis.x * this->_axis.y * olcos;
	number xzo = this->_axis.x * this->_axis.z * olcos;
	number yzo = this->_axis.y * this->_axis.z * olcos;
	number xsin = this->_axis.x * sintheta;
	number ysin = this->_axis.y * sintheta;
	number zsin = this->_axis.z * sintheta;

	LR_matrix<number> R(this->_axis.x * this->_axis.x * olcos + costheta,
						xyo - zsin,
						xzo + ysin,

						xyo + zsin,
						this->_axis.y * this->_axis.y * olcos + costheta,
						yzo - xsin,

						xzo - ysin,
						yzo + xsin,
						this->_axis.z * this->_axis.z * olcos + costheta);

	LR_vector<number> postrap = R * (_pos0 - _center) + _center;

	// we "mask" the resulting vector;
	number x = - this->_stiff * (pos.x - postrap.x) * _mask.x;
	number y = - this->_stiff * (pos.y - postrap.y) * _mask.y;
	number z = - this->_stiff * (pos.z - postrap.z) * _mask.z;

	return LR_vector<number>(x, y, z);
}

template class ConstantRateTorque<double>;
template class ConstantRateTorque<float>;

