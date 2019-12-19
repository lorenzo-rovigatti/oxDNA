/*
 * ForceFactory.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "ForceFactory.h"

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

using namespace std;

std::shared_ptr<ForceFactory> ForceFactory::_ForceFactoryPtr = nullptr;

ForceFactory::ForceFactory() {

}

ForceFactory::~ForceFactory() {

}

std::shared_ptr<ForceFactory> ForceFactory::instance() {
	// we can't use std::make_shared because ForceFactory's constructor is private
	if(_ForceFactoryPtr == nullptr) {
		_ForceFactoryPtr = std::shared_ptr<ForceFactory>(new ForceFactory());
	}

	return _ForceFactoryPtr;
}

void ForceFactory::add_force(input_file &inp, BaseParticle **particles, int N, bool is_CUDA, BaseBox * box_ptr) {

	string type_str;
	getInputString(&inp, "type", type_str, 1);

	ForcePtr extF;

	if(type_str.compare("string") == 0) extF = std::make_shared<ConstantRateForce>();
	else if(type_str.compare("sawtooth") == 0) extF = std::make_shared<SawtoothForce>();
	else if(type_str.compare("twist") == 0) extF = std::make_shared<ConstantRateTorque>();
	else if(type_str.compare("trap") == 0) extF = std::make_shared<MovingTrap>();
	else if(type_str.compare("repulsion_plane") == 0) extF = std::make_shared<RepulsionPlane>();
	else if(type_str.compare("repulsion_plane_moving") == 0) extF = std::make_shared<RepulsionPlaneMoving>();
	else if(type_str.compare("mutual_trap") == 0) extF = std::make_shared<MutualTrap>();
	else if(type_str.compare("lowdim_trap") == 0) extF = std::make_shared<LowdimMovingTrap>();
	else if(type_str.compare("constant_trap") == 0) extF = std::make_shared<ConstantTrap>();
	else if(type_str.compare("sphere") == 0) extF = std::make_shared<RepulsiveSphere>();
	else if(type_str.compare("sphere_smooth") == 0) extF = std::make_shared<RepulsiveSphereSmooth>();
	else if(type_str.compare("com") == 0) extF = std::make_shared<COMForce>();
	else if(type_str.compare("LJ_wall") == 0) extF = std::make_shared<LJWall>();
	else if(type_str.compare("hard_wall") == 0) extF = std::make_shared<HardWall>();
	else if(type_str.compare("alignment_field") == 0) extF = std::make_shared<AlignmentField>();
	else if(type_str.compare("generic_central_force") == 0) extF = std::make_shared<GenericCentralForce>();
	else if(type_str.compare("LJ_cone") == 0) extF = std::make_shared<LJCone>();
	else throw oxDNAException("Invalid force type `%s\'", type_str.c_str());

	string group = string("default");
	getInputString(&inp, "group_name", group, 0);

	extF->get_settings(inp);
	extF->init(particles, N, box_ptr); // here the force is added to the particle
	extF->set_group_name(group);
}

void ForceFactory::read_external_forces(std::string external_filename, BaseParticle ** particles, int N, bool is_CUDA, BaseBox * box) {
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

		ForceFactory::instance()->add_force(input, particles, N, is_CUDA, box);

		cleanInputFile (&input);
		fclose (temp);
	}
	OX_LOG(Logger::LOG_INFO, "   Force file parsed", external_filename.c_str());
}
