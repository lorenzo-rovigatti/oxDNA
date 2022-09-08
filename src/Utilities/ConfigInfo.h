/*
 * ConfigInfo.h
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#ifndef SRC_UTILITIES_CONFIGINFO_H_
#define SRC_UTILITIES_CONFIGINFO_H_

#define CONFIG_INFO ConfigInfo::instance()

#include "../defs.h"
#include "oxDNAException.h"
#include "FlattenedConfigInfo.h"

#include <string>
#include <memory>
#include <vector>

#include <functional>
#include <map>

class BaseInteraction;
class BaseParticle;
class BaseList;
class BaseBox;
class BaseForce;
class BaseObservable;
struct input_file;
struct Molecule;

/**
 * @brief Utility class. It is used by observables to have access to SimBackend's private members.
 *
 * An instance to this class must be provided to initialise observables. In general, observables need
 * to have access to backends' private members, but since oxDNA is always in development, we are not sure
 * about what observables should or should not have access to. This is not a very clean way of
 * doing it, but having a class like this allows us to easily pass information to the observables without having
 * to change the overall design.
 */
class ConfigInfo {
private:
	ConfigInfo(std::vector<BaseParticle *> *ps, std::vector<std::shared_ptr<Molecule>> *mols);
	ConfigInfo() = delete;

	static std::shared_ptr<ConfigInfo> _config_info;

	/// A map associating list of callbacks to events
	std::map<std::string, std::vector<std::function<void()>>> _event_callbacks;

	number _temperature = 0.;

	FlattenedConfigInfo _flattened_conf;

public:
	virtual ~ConfigInfo();

	/**
	 * @brief Allows one to initialize everything at once.
	 *
	 * @param i the interaction object
	 * @param info
	 * @param l the list object
	 * @param abox the box object
	 */
	void set(BaseInteraction *i, std::string *info, BaseList *l, BaseBox *abox);

	void subscribe(std::string event, std::function<void()> callback);

	void notify(std::string event);

	int N() {
		return particles_pointer->size();
	}

	void update_temperature(number new_T);

	number temperature() {
		return _temperature;
	}

	void add_force_to_particles(std::shared_ptr<BaseForce> force, std::vector<int> particle_ids, std::string force_type);

	std::shared_ptr<BaseForce> get_force_by_id(std::string id);

	std::shared_ptr<BaseObservable> get_observable_by_id(std::string id);

	/**
	 * @brief Returns a pointer to the actual object. Static method to enforce the singleton pattern.
	 *
	 * @return Pointer to the ConfigInfo object
	 */
	static std::shared_ptr<ConfigInfo> instance();

	static void init(std::vector<BaseParticle *> *ps, std::vector<std::shared_ptr<Molecule>> *mols);

	static void clear();

	std::vector<BaseParticle *> &particles() {
		return *particles_pointer;
	}

	std::vector<std::shared_ptr<Molecule>> &molecules() {
		return *molecules_pointer;
	}

	const FlattenedConfigInfo &flattened_conf();

	/// Pointer to the array that stores all the particles' information.
	std::vector<BaseParticle *> *particles_pointer;

	/// Pointer to the vector that stores the molecules
	std::vector<std::shared_ptr<Molecule>> *molecules_pointer;

	/// Vector storing all the external forces
	std::vector<std::shared_ptr<BaseForce>> forces;

	/// Vector storing all the observables
	std::vector<std::shared_ptr<BaseObservable>> observables;

	/// Used to compute all different kinds of interaction energies (total, partial, between two particles, etc.).
	BaseInteraction *interaction = nullptr;

	/// Used by BackendInfo to print backend-related information such as Monte Carlo acceptance ratios.
	std::string *backend_info = nullptr;

	/// Pointer to lists
	BaseList *lists = nullptr;

	/// Pointer to box object
	BaseBox *box = nullptr;

	/// Current simulation step
	long long int curr_step = 0;

	input_file *sim_input = nullptr;
};

inline std::shared_ptr<ConfigInfo> ConfigInfo::instance() {
	if(_config_info == nullptr) throw oxDNAException("Trying to access an uninitialised ConfigInfo object");

	return _config_info;
}

#endif /* SRC_UTILITIES_CONFIGINFO_H_ */
