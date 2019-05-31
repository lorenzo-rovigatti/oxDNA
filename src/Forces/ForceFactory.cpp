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
#include "SawtoothForce.h"
#include "ConstantRateTorque.h"
#include "ConstantTrap.h"
#include "LowdimMovingTrap.h"
#include "MovingTrap.h"
#include "MutualTrap.h"
#include "RepulsionPlane.h"
#include "RepulsionPlaneMoving.h"
#include "LJWall.h"
#include "HardWall.h"
#include "RepulsiveSphere.h"
#include "RepulsiveSphereSmooth.h"
#include "AlignmentField.h"
#include "GenericCentralForce.h"
#include "LJCone.h"
#include <fstream>
#include <sstream>
#include "LJCone.h"

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
void ForceFactory<number>::add_force(input_file &inp, BaseParticle<number> **particles, int N, bool is_CUDA, BaseBox<number> * box_ptr) {

	string type_str;
	getInputString (&inp, "type", type_str, 1);

	BaseForce<number> * extF;
	
	if (type_str.compare("string") == 0) extF = new ConstantRateForce<number>();
	else if (type_str.compare("sawtooth") == 0) extF = new SawtoothForce<number>();
	else if (type_str.compare("twist") == 0) extF = new ConstantRateTorque<number>();
	else if (type_str.compare("trap") == 0) extF = new MovingTrap<number>();
	else if (type_str.compare("repulsion_plane") == 0) extF = new RepulsionPlane<number>();
	else if (type_str.compare("repulsion_plane_moving") == 0) extF = new RepulsionPlaneMoving<number>();
	else if (type_str.compare("mutual_trap") == 0) extF = new MutualTrap<number>();
	else if (type_str.compare("lowdim_trap") == 0) extF = new LowdimMovingTrap<number>();
	else if (type_str.compare("constant_trap") == 0) extF = new ConstantTrap<number>();
	else if (type_str.compare("sphere") == 0) extF = new RepulsiveSphere<number>();
	else if (type_str.compare("sphere_smooth") == 0) extF = new RepulsiveSphereSmooth<number>();
	else if (type_str.compare("com") == 0) extF = new COMForce<number>();
	else if (type_str.compare("LJ_wall") == 0) extF = new LJWall<number>();
	else if (type_str.compare("hard_wall") == 0) extF = new HardWall<number>();
	else if (type_str.compare("alignment_field") == 0) extF = new AlignmentField<number>();
	else if (type_str.compare("generic_central_force") == 0) extF = new GenericCentralForce<number>();
	else if (type_str.compare("LJ_cone") == 0) extF = new LJCone<number>();
	else throw oxDNAException ("Invalid force type `%s\'", type_str.c_str());

	string group = string("default");
	getInputString(&inp, "group_name", group, 0);
	
	extF->get_settings(inp);
	extF->init(particles, N, box_ptr); // here the force is added to the particle
	extF->set_group_name(group);

	_forces.push_back(extF);
}

template<typename number>
void ForceFactory<number>::read_external_forces(std::string external_filename, BaseParticle<number> ** particles, int N, bool is_CUDA,BaseBox<number> * box) {
	OX_LOG(Logger::LOG_INFO, "Parsing Force file %s", external_filename.c_str());

	//char line[512], typestr[512];
	int open, justopen, a;
	ifstream external(external_filename.c_str());

	if(!external.good ()) throw oxDNAException ("Can't read external_forces_file '%s'", external_filename.c_str());

	justopen = open = 0;
	a = external.get();
	bool is_commented = false;
	stringstream external_string;
	while(external.good()) {
		justopen = 0;
		switch(a) {
		case '#':
			is_commented = true;
			break;
		case '\n':
			is_commented = false;
			break;
		case '{':
			if(!is_commented) {
				open++;
				justopen = 1;
			}
			break;
		case '}':
			if(!is_commented) {
				if(justopen) throw oxDNAException ("Syntax error in '%s': nothing between parentheses", external_filename.c_str());
				open--;
			}
			break;
		default:
			break;
		}

		if(!is_commented) external_string << (char)a;
		if(open > 1 || open < 0) throw oxDNAException ("Syntax error in '%s': parentheses do not match", external_filename.c_str());
		a = external.get();
	}
	external.close();

	external_string.clear();
	external_string.seekg(0, ios::beg);
	a = external_string.get();
	while(external_string.good()) {
		while (a != '{' && external_string.good()) a = external_string.get();
		if(!external_string.good()) break;
		// this function create a temporary file which is destroyed upon calling fclose
		// the temporary file is opened with "wb+" flags
		FILE *temp = tmpfile();
		OX_LOG(Logger::LOG_INFO, "   Using temporary file");
		a = external_string.get();
		while (a != '}' && external_string.good()) {
			fprintf (temp, "%c", a);
			a = external_string.get();
		}
		rewind(temp);
		input_file input;
		loadInput(&input, temp);

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
