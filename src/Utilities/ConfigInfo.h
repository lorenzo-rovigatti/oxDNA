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

 class IBaseInteraction;
 class BaseParticle;
 class BaseList;
 class BaseBox;

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
	static ConfigInfo *_config_info;

	ConfigInfo();
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
	void set(BaseParticle **p, IBaseInteraction *i, int *Nn, std::string *info, BaseList *l, BaseBox *abox);

	/**
	 * @brief Returns a reference to the actual object. Static method to enforce the singleton pattern.
	 *
	 * @return Reference to the ConfigInfo object
	 */
	static ConfigInfo &ref_instance();

	/**
	 * @brief Returns a pointer to the actual object. Static method to enforce the singleton pattern.
	 *
	 * @return Pointer to the ConfigInfo object
	 */
	static ConfigInfo *instance();

	static void init();
	static void clear();

	///Pointer to the array which stores all the particles' information.
	BaseParticle **particles;


	/// Used to compute all different kinds of interaction energies (total, partial, between two particles, etc.).
	IBaseInteraction *interaction;

	/// Number of particles.
	int *N;

	/// Used by BackendInfo to print backend-related information such as Monte Carlo acceptance ratios.
	std::string *backend_info;

	/// Pointer to lists
	BaseList *lists;

	/// Pointer to box object
	BaseBox *box;

	/// Current simulation step
	long long int curr_step;
};


inline ConfigInfo &ConfigInfo::ref_instance() {
	return *instance();
}


inline ConfigInfo *ConfigInfo::instance() {
	if(_config_info == NULL) throw oxDNAException("Trying to access an uninitialised ConfigInfo object");

	return _config_info;
}

#endif /* SRC_UTILITIES_CONFIGINFO_H_ */
