/*
 * MoveFactory.cpp
 *
 *  Created on: Jun 25, 2014
 *      Author: flavio
 */

#include "MoveFactory.h"

// we may want to add this later...
//#include "../PluginManagement/PluginManager.h"

#include "MCTras.h"
#include "MCRot.h"

MoveFactory::MoveFactory() {

}

MoveFactory::~MoveFactory() {

}

template<typename number>
BaseMove<number> * MoveFactory::make_move (input_file &inp, input_file &sim_inp, ConfigInfo<number> * Info) {
	BaseMove<number> * ret = NULL;

	std::string move_type;
	getInputString (&inp, "type", move_type, 0);
	
	if (!strncmp (move_type.c_str(), "rotation", move_type.length() + 2)) ret = new MCRot<number>(Info);
	else if (!strncmp (move_type.c_str(), "traslation", move_type.length() + 2)) ret = new MCTras<number>(Info);
	else throw oxDNAException ("(MoveFactory.cpp) Unknown move type %s. Aborting", move_type.c_str());

	/* // plugin management o be added later on... NEEDS to be fixed, copied and pasted from interactionfactory.cpp
	else {
		BaseMove<number> *res = PluginManager::instance()->get_interaction<number>(inter_type);
		if(res == NULL) throw oxDNAException ("Move '%s' not found. Aborting", inter_type);
		return res;
	}*/
	ret->get_settings (inp, sim_inp);
	
	return ret;
}

template BaseMove<float> *MoveFactory::make_move<float>(input_file &inp, input_file &sim_inp, ConfigInfo<float> * Info);
template BaseMove<double> *MoveFactory::make_move<double>(input_file &inp, input_file &sim_inp, ConfigInfo<double> * Info);

