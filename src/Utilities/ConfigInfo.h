/*
 * ConfigInfo.h
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#ifndef SRC_UTILITIES_CONFIGINFO_H_
#define SRC_UTILITIES_CONFIGINFO_H_

#define CONFIG_INFO ConfigInfo::instance()

#include "oxDNAException.h"

#include <string>
#include <memory>
#include <vector>

class IBaseInteraction;
class BaseParticle;
class BaseList;
class BaseBox;
class input_file;

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
	static std::shared_ptr<ConfigInfo> _config_info;

	ConfigInfo(std::vector<BaseParticle *> &ps);
	ConfigInfo() = delete;
public:
	virtual ~ConfigInfo();

	/**
	 * @brief Allows one to initialize everything at once.
	 *
	 * @param p
	 * @param i
	 * @param Nn number of particles
	 * @param info
	 * @param l pointer to list object
	 */
	void set(IBaseInteraction *i, int *Nn, std::string *info, BaseList *l, BaseBox *abox);

	/**
	 * @brief Returns a pointer to the actual object. Static method to enforce the singleton pattern.
	 *
	 * @return Pointer to the ConfigInfo object
	 */
	static ConfigInfo *instance();

	static void init(std::vector<BaseParticle *> &ps);

	/// Pointer to the array which stores all the particles' information.
	std::vector<BaseParticle *> &particles;

	/// Used to compute all different kinds of interaction energies (total, partial, between two particles, etc.).
	IBaseInteraction *interaction = nullptr;

	/// Number of particles.
	int *N = nullptr;

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


inline ConfigInfo *ConfigInfo::instance() {
	if(_config_info == nullptr) throw oxDNAException("Trying to access an uninitialised ConfigInfo object");

	return _config_info.get();
}

#endif /* SRC_UTILITIES_CONFIGINFO_H_ */
