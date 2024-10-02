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
#include "MoleculeVolumeMove.h"
#include "Pivot.h"
#include "VMMC.h"
#include "ShapeMove.h"
#include "RotateSite.h"
//#include "MCMovePatchyShape.h"

MovePtr MoveFactory::make_move(input_file &inp, input_file &sim_inp) {
	MovePtr ret = NULL;

	std::string move_type;
	getInputString(&inp, "type", move_type, 0);

	if(!move_type.compare("rotation")) ret = std::make_shared<MCRot>();
	else if(!move_type.compare("translation")) ret = std::make_shared<MCTras>();
	else if(!move_type.compare("volume")) ret = std::make_shared<VolumeMove>();
	else if(!move_type.compare("molecule_volume")) ret = std::make_shared<MoleculeVolumeMove>();
	else if(!move_type.compare("VMMC")) ret = std::make_shared<VMMC>();
	else if(!move_type.compare("shape")) ret = std::make_shared<ShapeMove>();
	else if(!move_type.compare("rotate_site")) ret = std::make_shared<RotateSite>();
	else if(!move_type.compare("pivot")) ret = std::make_shared<Pivot>();
	//else if (!move_type.compare("PatchyShape")) ret = std::make_shared<MCMovePatchyShape(Info);
	else {
		ret = PluginManager::instance()->get_move(move_type);
		if(ret == NULL) throw oxDNAException("Move object '%s' not found. Aborting", move_type.c_str());
	}
	ret->get_settings(inp, sim_inp);

	return ret;
}
