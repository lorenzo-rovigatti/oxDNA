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
#include "HBEnergy.h"
#include "BackendInfo.h"
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
#include "NematicS.h"
#include "UnstackedList.h"
#include "PlectonemePosition.h"

#include "Configurations/PdbOutput.h"
#include "Configurations/ChimeraOutput.h"
#include "Configurations/Configuration.h"
#include "Configurations/BinaryConfiguration.h"
#include "Configurations/TclOutput.h"
#include "Configurations/TEPtclOutput.h"
#include "Configurations/TEPxyzOutput.h"
#include "Configurations/JordanOutput.h"

ObservableFactory::ObservableFactory() {

}

ObservableFactory::~ObservableFactory() {

}

template<typename number>
BaseObservable<number> *ObservableFactory::make_observable(input_file &obs_inp, input_file &sim_inp) {
	char obs_type[512];
	getInputString(&obs_inp, "type", obs_type, 1);

	BaseObservable<number> *res = NULL;

	if(!strncasecmp(obs_type, "step", 512)) res = new Step<number>();
	else if(!strncasecmp(obs_type, "potential_energy", 512)) res = new PotentialEnergy<number>();
	else if(!strncasecmp(obs_type, "kinetic_energy", 512)) res = new KineticEnergy<number>();
	else if(!strncasecmp(obs_type, "total_energy", 512)) res = new TotalEnergy<number>();
	else if(!strncasecmp(obs_type, "hb_energy", 512)) res = new HBEnergy<number>();
	else if(!strncasecmp(obs_type, "backend_info", 512)) res = new BackendInfo<number>();
	else if(!strncasecmp(obs_type, "pair_energy", 512)) res = new PairEnergy<number>();
	else if(!strncasecmp(obs_type, "pair_force", 512)) res = new PairForce<number>();
	else if(!strncasecmp(obs_type, "hb_list", 512)) res = new HBList<number>();
	else if(!strncasecmp(obs_type, "order_parameters", 512)) res = new OrderParameterValues<number>();
	else if(!strncasecmp(obs_type, "strandwise_bonds", 512)) res = new StrandwiseBonds<number>();
	else if(!strncasecmp(obs_type, "force_energy", 512)) res = new ForceEnergy<number>();
	else if(!strncasecmp(obs_type, "configuration", 512)) res = new Configuration<number>();
	else if(!strncasecmp(obs_type, "binary_configuration", 512)) res = new BinaryConfiguration<number>();
	else if(!strncasecmp(obs_type, "tcl_configuration", 512)) res = new TclOutput<number>();
	else if(!strncasecmp(obs_type, "pressure", 512)) res = new Pressure<number>();
	else if(!strncasecmp(obs_type, "density", 512)) res = new Density<number>();
	else if(!strncasecmp(obs_type, "density_profile", 512)) res = new DensityProfile<number>();
	else if(!strncasecmp(obs_type, "rdf", 512)) res = new Rdf<number>();
	else if(!strncasecmp(obs_type, "particle_position", 512)) res = new ParticlePosition<number>();
	else if(!strncasecmp(obs_type, "distance", 512)) res = new Distance<number>();
	else if(!strncasecmp(obs_type, "pdb_configuration", 512)) res = new PdbOutput<number>();
	else if(!strncasecmp(obs_type, "chimera_script", 512)) res = new ChimeraOutput<number>();
	else if(!strncasecmp(obs_type, "coax_variables", 512)) res = new CoaxVariables<number>();
	else if(!strncasecmp(obs_type, "pitch", 512)) res = new Pitch<number>();
	else if(!strncasecmp(obs_type, "salt_extrapolation", 512)) res = new SaltExtrapolation<number>();
	else if(!strncasecmp(obs_type, "external_torque", 512)) res = new ExternalTorque<number>();
	else if(!strncasecmp(obs_type, "mean_vector_cosine", 512)) res = new MeanVectorCosine<number>();
	else if(!strncasecmp(obs_type, "vector_angle", 512)) res = new VectorAngle<number>();
	else if(!strncasecmp(obs_type, "TEPtcl_configuration", 512)) res = new TEPtclOutput<number>();
	else if(!strncasecmp(obs_type, "TEPxyz_configuration", 512)) res = new TEPxyzOutput<number>();
	else if(!strncasecmp(obs_type, "checkpoint", 512)) res = new Checkpoint<number>();
	else if(!strncasecmp(obs_type, "contacts", 512)) res = new Contacts<number>();
	else if(!strncasecmp(obs_type, "writhe", 512)) res = new Writhe<number>();
	else if(!strncasecmp(obs_type, "nematic_s", 512)) res = new NematicS<number>();
	else if(!strncasecmp(obs_type, "jordan_conf", 512)) res = new JordanOutput<number>();
	else if(!strncasecmp(obs_type, "unstacked_list", 512)) res = new UnstackedList<number>();
	else if(!strncasecmp(obs_type, "plectoneme_position", 512)) res = new PlectonemePosition<number>();
	else {
		res = PluginManager::instance()->get_observable<number>(obs_type);
		if(res == NULL) throw oxDNAException ("Observable '%s' not found. Aborting", obs_type);
	}
	
	res->get_settings(obs_inp, sim_inp);

	return res;
}

template BaseObservable<float> *ObservableFactory::make_observable<float>(input_file &obs_inp, input_file &sim_inp);
template BaseObservable<double> *ObservableFactory::make_observable<double>(input_file &obs_inp, input_file &sim_inp);

