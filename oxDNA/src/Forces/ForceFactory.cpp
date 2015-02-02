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

#include <fstream>

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
void ForceFactory<number>::read_external_forces(std::string external_filename, BaseParticle<number> ** particles, int N, bool is_CUDA, number * box) {
	OX_LOG(Logger::LOG_INFO, "Parsing Force file %s", external_filename.c_str());

	//char line[512], typestr[512];
	int open, justopen, a;
	ifstream external(external_filename.c_str());

	if(!external.good ()) throw oxDNAException ("Can't read external_forces_file '%s'", external_filename.c_str());

	justopen = open = 0;
	a = external.get();
	while(external.good()) {
		justopen = 0;
		if (a == '{') {
			open ++;
			justopen = 1;
		}
		if (a == '}') {
			if (justopen) throw oxDNAException ("Syntax error in '%s': nothing between parentheses", external_filename.c_str());
			open --;
		}
		if (open > 1 || open < 0) throw oxDNAException ("Syntax error in '%s': parentheses do not match", external_filename.c_str());
		a = external.get();
	}
	external.clear();
	external.seekg(0, ios::beg);

	a = external.get();
	while(external.good()) {
		while (a != '{' && external.good()) a = external.get();
		if(!external.good()) break;
		// this function create a temporary file which is destroyed upon calling fclose
		// the temporary file is opened with "wb+" flags
		FILE *temp = tmpfile();
		OX_LOG(Logger::LOG_INFO, "   Using temporary file");
		a = external.get();
		while (a != '}' && external.good()) {
			fprintf (temp, "%c", a);
			a = external.get();
		}
		rewind(temp);
		input_file input;
		loadInput(&input, temp);

		//ForceFactory::add_force<number>(input, _particles, _N, _is_CUDA_sim, &_box_side);
		ForceFactory<number>::instance()->add_force(input, particles, N, is_CUDA, box);

		cleanInputFile (&input);
		fclose (temp);
	}
	OX_LOG(Logger::LOG_INFO, "   Force file parsed", external_filename.c_str());
}

template<typename number>
void ForceFactory<number>::clear () {
	if (_ForceFactoryPtr != NULL) delete _ForceFactoryPtr; // this if is to make oxDNA not segfault if an exception is thrown in the code

	return;
}

template class ForceFactory<float>;
template class ForceFactory<double>;
