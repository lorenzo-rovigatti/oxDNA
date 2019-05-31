/*
 * MoveFactory.cpp
 *
 *  Created on: Jun 25, 2014
 *      Author: flavio
 */

#include "MoveFactory.h"

// we may want to add this later...
#include "../../PluginManagement/PluginManager.h"

#include "MCTras.h"
#include "MCRot.h"
#include "VolumeMove.h"
#include "VMMC.h"
#include "ShapeMove.h"
#include "RotateSite.h"
//#include "MCMovePatchyShape.h"

MoveFactory::MoveFactory() {

}

MoveFactory::~MoveFactory() {

}

template<typename number>
BaseMove<number> * MoveFactory::make_move (input_file &inp, input_file &sim_inp) {
	BaseMove<number> * ret = NULL;

	std::string move_type;
	getInputString (&inp, "type", move_type, 0);
	
	if (!move_type.compare("rotation")) ret = new MCRot<number>();
	else if (!move_type.compare("translation")) ret = new MCTras<number>();
	else if (!move_type.compare("volume")) ret = new VolumeMove<number>();
	else if (!move_type.compare("VMMC")) ret = new VMMC<number>();
	else if (!move_type.compare("shape")) ret = new ShapeMove<number>();
	else if (!move_type.compare("rotate_site")) ret = new RotateSite<number>();
	//else if (!move_type.compare("PatchyShape")) ret = new MCMovePatchyShape<number>(Info);
	else {
			ret = PluginManager::instance()->get_move<number>(move_type);
			if(ret == NULL) throw oxDNAException ("Move object '%s' not found. Aborting", move_type.c_str());

	}


    //throw oxDNAException ("(MoveFactory.cpp) Unknown move type %s. Aborting", move_type.c_str());

	/* // plugin management o be added later on... NEEDS to be fixed, copied and pasted from interactionfactory.cpp
	else {
		BaseMove<number> *res = PluginManager::instance()->get_interaction<number>(inter_type);
		if(res == NULL) throw oxDNAException ("Move '%s' not found. Aborting", inter_type);
		return res;
	}*/
	ret->get_settings (inp, sim_inp);
	
	return ret;
}

template BaseMove<float> *MoveFactory::make_move<float>(input_file &inp, input_file &sim_inp);
template BaseMove<double> *MoveFactory::make_move<double>(input_file &inp, input_file &sim_inp);

