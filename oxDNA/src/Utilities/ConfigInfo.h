/*
 * ConfigInfo.h
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#ifndef SRC_UTILITIES_CONFIGINFO_H_
#define SRC_UTILITIES_CONFIGINFO_H_

#define CONFIG_INFO ConfigInfo<number>::instance()

#include "oxDNAException.h"

#include <string>

template <typename number> class IBaseInteraction;
template <typename number> class BaseParticle;
template <typename number> class BaseList;
template <typename number> class BaseBox;

/**
 * @brief Utility class. It is used by observables to have access to SimBackend's private members.
 *
 * An instance to this class must be provided to initialise observables. In general, observables need
 * to have access to backends' private members, but since oxDNA is always in development, we are not sure
 * about what observables should or should not have access to. This is not a very clean way of
 * doing it, but having a class like this allows us to easily pass information to the observables without having
 * to change the overall design.
 */
template<typename number>
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
	void set(BaseParticle<number> **p, IBaseInteraction<number> *i, int *Nn, std::string *info, BaseList<number> *l, BaseBox<number> *abox);

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
	BaseParticle<number> **particles;


	/// Used to compute all different kinds of interaction energies (total, partial, between two particles, etc.).
	IBaseInteraction<number> *interaction;

	/// Number of particles.
	int *N;

	/// Used by BackendInfo to print backend-related information such as Monte Carlo acceptance ratios.
	std::string *backend_info;

	/// Pointer to lists
	BaseList<number> *lists;

	/// Pointer to box object
	BaseBox<number> *box;

	/// Current simulation step
	long long int curr_step;
};

template<typename number>
inline ConfigInfo<number> &ConfigInfo<number>::ref_instance() {
	return *instance();
}

template<typename number>
inline ConfigInfo<number> *ConfigInfo<number>::instance() {
	if(_config_info == NULL) throw oxDNAException("Trying to access an uninitialised ConfigInfo object");

	return _config_info;
}

#endif /* SRC_UTILITIES_CONFIGINFO_H_ */
