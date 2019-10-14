/*
 * ObservableFactory.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "ObservableFactory.h"

#include "../PluginManagement/PluginManager.h"

#include "Step.h"
#include "PotentialEnergy.h"
#include "KineticEnergy.h"
#include "TotalEnergy.h"
#include "Configurations/Configuration.h"
#include "BackendInfo.h"

/*#include "HBEnergy.h"
#include "PairEnergy.h"
#include "PairForce.h"
#include "OrderParameterValues.h"
#include "HBList.h"
#include "StrandwiseBonds.h"
#include "ForceEnergy.h"
#include "Pressure.h"
#include "Density.h"
#include "DensityProfile.h"
#include "ParticlePosition.h"
#include "Rdf.h"
#include "Distance.h"
#include "CoaxVariables.h"
#include "Pitch.h"
#include "SaltExtrapolation.h"
#include "ExternalTorque.h"
#include "MeanVectorCosine.h"
#include "VectorAngle.h"
#include "Checkpoint.h"
#include "Contacts.h"
#include "Writhe.h"
#include "UnstackedList.h"
#include "PlectonemePosition.h"
#include "StretchedBonds.h"
#include "StructureFactor.h"
#include "FormFactor.h"
#include "UnstackedList.h"
#include "PlectonemePosition.h"
#include "TEPPlectonemePosition.h"
#include "AverageEnergy.h"
#include "ContactMap.h"
#include "AllVectors.h"

#include "Configurations/PdbOutput.h"
#include "Configurations/ChimeraOutput.h"
#include "Configurations/BinaryConfiguration.h"
#include "Configurations/TclOutput.h"
#include "Configurations/TEPtclOutput.h"
#include "Configurations/TEPxyzOutput.h"
#include "Configurations/JordanOutput.h"*/

ObservableFactory::ObservableFactory() {

}

ObservableFactory::~ObservableFactory() {

}


BaseObservable *ObservableFactory::make_observable(input_file &obs_inp, input_file &sim_inp) {
	char obs_type[512];
	getInputString(&obs_inp, "type", obs_type, 1);

	BaseObservable *res = NULL;

	if(!strncasecmp(obs_type, "step", 512)) res = new Step();
	else if(!strncasecmp(obs_type, "potential_energy", 512)) res = new PotentialEnergy();
	else if(!strncasecmp(obs_type, "kinetic_energy", 512)) res = new KineticEnergy();
	else if(!strncasecmp(obs_type, "total_energy", 512)) res = new TotalEnergy();
	else if(!strncasecmp(obs_type, "backend_info", 512)) res = new BackendInfo();
	else if(!strncasecmp(obs_type, "configuration", 512)) res = new Configuration();

	/*else if(!strncasecmp(obs_type, "hb_energy", 512)) res = new HBEnergy();
	else if(!strncasecmp(obs_type, "pair_energy", 512)) res = new PairEnergy();
	else if(!strncasecmp(obs_type, "pair_force", 512)) res = new PairForce();
	else if(!strncasecmp(obs_type, "hb_list", 512)) res = new HBList();
	else if(!strncasecmp(obs_type, "order_parameters", 512)) res = new OrderParameterValues();
	else if(!strncasecmp(obs_type, "strandwise_bonds", 512)) res = new StrandwiseBonds();
	else if(!strncasecmp(obs_type, "force_energy", 512)) res = new ForceEnergy();
	else if(!strncasecmp(obs_type, "binary_configuration", 512)) res = new BinaryConfiguration();
	else if(!strncasecmp(obs_type, "tcl_configuration", 512)) res = new TclOutput();
	else if(!strncasecmp(obs_type, "pressure", 512)) res = new Pressure();
	else if(!strncasecmp(obs_type, "density", 512)) res = new Density();
	else if(!strncasecmp(obs_type, "density_profile", 512)) res = new DensityProfile();
	else if(!strncasecmp(obs_type, "rdf", 512)) res = new Rdf();
	else if(!strncasecmp(obs_type, "particle_position", 512)) res = new ParticlePosition();
	else if(!strncasecmp(obs_type, "distance", 512)) res = new Distance();
	else if(!strncasecmp(obs_type, "pdb_configuration", 512)) res = new PdbOutput();
	else if(!strncasecmp(obs_type, "chimera_script", 512)) res = new ChimeraOutput();
	else if(!strncasecmp(obs_type, "coax_variables", 512)) res = new CoaxVariables();
	else if(!strncasecmp(obs_type, "pitch", 512)) res = new Pitch();
	else if(!strncasecmp(obs_type, "salt_extrapolation", 512)) res = new SaltExtrapolation();
	else if(!strncasecmp(obs_type, "external_torque", 512)) res = new ExternalTorque();
	else if(!strncasecmp(obs_type, "mean_vector_cosine", 512)) res = new MeanVectorCosine();
	else if(!strncasecmp(obs_type, "vector_angle", 512)) res = new VectorAngle();
	else if(!strncasecmp(obs_type, "TEPtcl_configuration", 512)) res = new TEPtclOutput();
	else if(!strncasecmp(obs_type, "TEPxyz_configuration", 512)) res = new TEPxyzOutput();
	else if(!strncasecmp(obs_type, "checkpoint", 512)) res = new Checkpoint();
	else if(!strncasecmp(obs_type, "contacts", 512)) res = new Contacts();
	else if(!strncasecmp(obs_type, "writhe", 512)) res = new Writhe();
	else if(!strncasecmp(obs_type, "stretched", 512)) res = new StretchedBonds();
	else if(!strncasecmp(obs_type, "jordan_conf", 512)) res = new JordanOutput();
	else if(!strncasecmp(obs_type, "unstacked_list", 512)) res = new UnstackedList();
	else if(!strncasecmp(obs_type, "plectoneme_position", 512)) res = new PlectonemePosition();
	else if(!strncasecmp(obs_type, "TEP_plectoneme_position", 512)) res = new TEPPlectonemePosition();
	else if(!strncasecmp(obs_type, "Sq", 512)) res = new StructureFactor();
	else if(!strncasecmp(obs_type, "Pq", 512)) res = new FormFactor();
	else if(!strncasecmp(obs_type, "unstacked_list", 512)) res = new UnstackedList();
	else if(!strncasecmp(obs_type, "plectoneme_position", 512)) res = new PlectonemePosition();
	else if(!strncasecmp(obs_type, "average_energy", 512)) res = new AverageEnergy();
	else if(!strncasecmp(obs_type, "contact_map", 512)) res = new ContactMap();
	else if(!strncasecmp(obs_type, "all_vectors", 512)) res = new AllVectors();*/
	else {
		res = PluginManager::instance()->get_observable(obs_type);
		if(res == NULL) throw oxDNAException ("Observable '%s' not found. Aborting", obs_type);
	}
	
	res->get_settings(obs_inp, sim_inp);

	return res;
}
