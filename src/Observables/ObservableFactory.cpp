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
#include "HBEnergy.h"
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
#include "StressAutocorrelation.h"
#include "ExternalForce.h"

#include "Configurations/PdbOutput.h"
#include "Configurations/ChimeraOutput.h"
#include "Configurations/BinaryConfiguration.h"
#include "Configurations/TclOutput.h"
#include "Configurations/TEPtclOutput.h"
#include "Configurations/TEPxyzOutput.h"
#include "Configurations/JordanOutput.h"

#ifdef JSON_ENABLED
#include <nlohmann/json.hpp>
#endif

ObservablePtr ObservableFactory::make_observable(input_file &obs_inp) {
	char obs_type[512];
	getInputString(&obs_inp, "type", obs_type, 1);

	ObservablePtr res = nullptr;

	if(!strncasecmp(obs_type, "step", 512)) res = std::make_shared<Step>();
	else if(!strncasecmp(obs_type, "potential_energy", 512)) res = std::make_shared<PotentialEnergy>();
	else if(!strncasecmp(obs_type, "kinetic_energy", 512)) res = std::make_shared<KineticEnergy>();
	else if(!strncasecmp(obs_type, "total_energy", 512)) res = std::make_shared<TotalEnergy>();
	else if(!strncasecmp(obs_type, "backend_info", 512)) res = std::make_shared<BackendInfo>();
	else if(!strncasecmp(obs_type, "configuration", 512)) res = std::make_shared<Configuration>();
	else if(!strncasecmp(obs_type, "hb_energy", 512)) res = std::make_shared<HBEnergy>();
	else if(!strncasecmp(obs_type, "pair_energy", 512)) res = std::make_shared<PairEnergy>();
	else if(!strncasecmp(obs_type, "pair_force", 512)) res = std::make_shared<PairForce>();
	else if(!strncasecmp(obs_type, "hb_list", 512)) res = std::make_shared<HBList>();
	else if(!strncasecmp(obs_type, "order_parameters", 512)) res = std::make_shared<OrderParameterValues>();
	else if(!strncasecmp(obs_type, "strandwise_bonds", 512)) res = std::make_shared<StrandwiseBonds>();
	else if(!strncasecmp(obs_type, "force_energy", 512)) res = std::make_shared<ForceEnergy>();
	else if(!strncasecmp(obs_type, "binary_configuration", 512)) res = std::make_shared<BinaryConfiguration>();
	else if(!strncasecmp(obs_type, "tcl_configuration", 512)) res = std::make_shared<TclOutput>();
	else if(!strncasecmp(obs_type, "pressure", 512)) res = std::make_shared<Pressure>();
	else if(!strncasecmp(obs_type, "density", 512)) res = std::make_shared<Density>();
	else if(!strncasecmp(obs_type, "density_profile", 512)) res = std::make_shared<DensityProfile>();
	else if(!strncasecmp(obs_type, "rdf", 512)) res = std::make_shared<Rdf>();
	else if(!strncasecmp(obs_type, "particle_position", 512)) res = std::make_shared<ParticlePosition>();
	else if(!strncasecmp(obs_type, "distance", 512)) res = std::make_shared<Distance>();
	else if(!strncasecmp(obs_type, "pdb_configuration", 512)) res = std::make_shared<PdbOutput>();
	else if(!strncasecmp(obs_type, "chimera_script", 512)) res = std::make_shared<ChimeraOutput>();
	else if(!strncasecmp(obs_type, "coax_variables", 512)) res = std::make_shared<CoaxVariables>();
	else if(!strncasecmp(obs_type, "pitch", 512)) res = std::make_shared<Pitch>();
	else if(!strncasecmp(obs_type, "salt_extrapolation", 512)) res = std::make_shared<SaltExtrapolation>();
	else if(!strncasecmp(obs_type, "external_torque", 512)) res = std::make_shared<ExternalTorque>();
	else if(!strncasecmp(obs_type, "mean_vector_cosine", 512)) res = std::make_shared<MeanVectorCosine>();
	else if(!strncasecmp(obs_type, "vector_angle", 512)) res = std::make_shared<VectorAngle>();
	else if(!strncasecmp(obs_type, "TEPtcl_configuration", 512)) res = std::make_shared<TEPtclOutput>();
	else if(!strncasecmp(obs_type, "TEPxyz_configuration", 512)) res = std::make_shared<TEPxyzOutput>();
	else if(!strncasecmp(obs_type, "checkpoint", 512)) res = std::make_shared<Checkpoint>();
	else if(!strncasecmp(obs_type, "contacts", 512)) res = std::make_shared<Contacts>();
	else if(!strncasecmp(obs_type, "writhe", 512)) res = std::make_shared<Writhe>();
	else if(!strncasecmp(obs_type, "stretched", 512)) res = std::make_shared<StretchedBonds>();
	else if(!strncasecmp(obs_type, "jordan_conf", 512)) res = std::make_shared<JordanOutput>();
	else if(!strncasecmp(obs_type, "unstacked_list", 512)) res = std::make_shared<UnstackedList>();
	else if(!strncasecmp(obs_type, "plectoneme_position", 512)) res = std::make_shared<PlectonemePosition>();
	else if(!strncasecmp(obs_type, "TEP_plectoneme_position", 512)) res = std::make_shared<TEPPlectonemePosition>();
	else if(!strncasecmp(obs_type, "Sq", 512)) res = std::make_shared<StructureFactor>();
	else if(!strncasecmp(obs_type, "Pq", 512)) res = std::make_shared<FormFactor>();
	else if(!strncasecmp(obs_type, "unstacked_list", 512)) res = std::make_shared<UnstackedList>();
	else if(!strncasecmp(obs_type, "plectoneme_position", 512)) res = std::make_shared<PlectonemePosition>();
	else if(!strncasecmp(obs_type, "average_energy", 512)) res = std::make_shared<AverageEnergy>();
	else if(!strncasecmp(obs_type, "contact_map", 512)) res = std::make_shared<ContactMap>();
	else if(!strncasecmp(obs_type, "all_vectors", 512)) res = std::make_shared<AllVectors>();
	else if(!strncasecmp(obs_type, "stress_autocorrelation", 512)) res = std::make_shared<StressAutocorrelation>();
	else if(!strncasecmp(obs_type, "external_force", 512)) res = std::make_shared<ExternalForce>();
	else {
		res = PluginManager::instance()->get_observable(obs_type);
		if(res == NULL) throw oxDNAException("Observable '%s' not found. Aborting", obs_type);
	}

	res->get_settings(obs_inp, *(CONFIG_INFO->sim_input));

	return res;
}

std::vector<ObservableOutputPtr> ObservableFactory::make_observables(std::string prefix) {
	std::vector<ObservableOutputPtr> result;

	std::string obs_input_base = Utils::sformat("%sdata_output_", prefix.c_str());
	std::string obs_key = Utils::sformat("%sobservables_file", prefix.c_str());

	int i = 1;
	bool found = true;
	while(found) {
		stringstream ss;
		ss << obs_input_base << i;
		string obs_string;
		if(getInputString(CONFIG_INFO->sim_input, ss.str().c_str(), obs_string, 0) == KEY_FOUND) {
			auto new_obs_out = std::make_shared<ObservableOutput>(obs_string);
			result.push_back(new_obs_out);
		}
		else {
			found = false;
		}

		i++;
	}

	std::string obs_filename;
	if(getInputString(CONFIG_INFO->sim_input, obs_key.c_str(), obs_filename, 0) == KEY_FOUND) {
#ifdef JSON_ENABLED
		OX_LOG(Logger::LOG_INFO, "Parsing JSON observable file %s", obs_filename.c_str());

		ifstream external(obs_filename.c_str());
		if(!external.good ()) {
			throw oxDNAException ("Can't read obs_filename '%s'", obs_filename.c_str());
		}

		nlohmann::json my_json = nlohmann::json::parse(external);

		for(auto &obs_json : my_json) {
			input_file obs_input;
			// extract the array containing details about the columns
			nlohmann::json cols_json;
			if(obs_json.contains("cols")) {
				cols_json = obs_json["cols"];
				obs_json.erase("cols");
			}

			// initialise the observable output
			obs_input.init_from_json(obs_json);
			std::string obs_string = obs_input.to_string();
			auto new_obs_out = std::make_shared<ObservableOutput>(obs_string);

			// and add the observables initialised from the cols array
			for(auto &col_json : cols_json) {
				input_file col_input;
				col_input.init_from_json(col_json);
				auto new_obs = make_observable(col_input);
				new_obs_out->add_observable(new_obs);
			}

			result.push_back(new_obs_out);
		}

		external.close();
		OX_LOG(Logger::LOG_INFO, "   Observable file parsed");
#else
		throw oxDNAException("Cannot parse the JSON observable file %s since JSON has been disabled at compile time", obs_filename.c_str());
#endif
	}

	return result;
}
