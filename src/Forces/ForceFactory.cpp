/*
 * ForceFactory.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "ForceFactory.h"

// forces to include
#include "COMForce.h"
#include "ConstantRateForce.h"
#include "ConstantRateTorque.h"
#include "ConstantTrap.h"
#include "LowdimMovingTrap.h"
#include "MovingTrap.h"
#include "MutualTrap.h"
#include "RepulsionPlane.h"
#include "RepulsionPlaneMoving.h"
#include "RepulsiveSphere.h"

using namespace std;

template <typename number>
ForceFactory<number> * ForceFactory<number>::_ForceFactoryPtr = NULL;

template <typename number>
ForceFactory<number>::ForceFactory() {

}

template <typename number>
ForceFactory<number>::~ForceFactory() {
	for (unsigned int i =0; i < _forces.size(); i ++) delete _forces[i];
}

template <typename number>
ForceFactory<number> * ForceFactory<number>::instance() {
	if (_ForceFactoryPtr == NULL) _ForceFactoryPtr = new ForceFactory();

	return _ForceFactoryPtr;
}

template<typename number>
void ForceFactory<number>::add_force(input_file &inp, BaseParticle<number> **particles, int N, bool is_CUDA, number * box_side_ptr) {

	string type_str;
	getInputString (&inp, "type", type_str, 1);

	BaseForce<number> * extF;
	
	if (type_str.compare("string") == 0) extF = new ConstantRateForce<number> ();
	else if (type_str.compare("twist") == 0) extF = new ConstantRateTorque<number> ();
	else if (type_str.compare("trap") == 0) extF = new MovingTrap<number> ();
	else if (type_str.compare("repulsion_plane") == 0) extF = new RepulsionPlane<number> ();
	else if (type_str.compare("repulsion_plane_moving") == 0) extF = new RepulsionPlaneMoving<number> ();
	else if (type_str.compare("mutual_trap") == 0) extF = new MutualTrap<number> ();
	else if (type_str.compare("lowdim_trap") == 0) extF = new LowdimMovingTrap<number> ();
	else if (type_str.compare("constant_trap") == 0) extF = new ConstantTrap<number> ();
	else if (type_str.compare("sphere") == 0) extF = new RepulsiveSphere<number> ();
	else if (type_str.compare("com") == 0) extF = new COMForce<number> ();
	else throw oxDNAException ("Invalid force type `%s\'", type_str.c_str());

	string group = string("default");
	getInputString(&inp, "group_name", group, 0);
	
	extF->get_settings(inp);
	extF->init(particles, N, box_side_ptr); // here the force is added to the particle
	extF->set_group_name(group);

	_forces.push_back(extF);
}

template<typename number>
void ForceFactory<number>::clear () {
	if (_ForceFactoryPtr != NULL) delete _ForceFactoryPtr; // this if is to make oxDNA not segfault if an exception is thrown in the code

	return;
}

template class ForceFactory<float>;
template class ForceFactory<double>;
