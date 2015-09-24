/*
 * BaseObservable.h
 *
 *  Created on: Feb 11, 2013
 *      Author: rovigatti
 */

#ifndef BASEOBSERVABLE_H_
#define BASEOBSERVABLE_H_

#include "../defs.h"
#include "../Interactions/BaseInteraction.h"
#include "../Lists/BaseList.h"

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
public:
	/**
	 * @brief Default constructor. It does nothing but zero initialization
	 */
	ConfigInfo() : particles(NULL), box_side(NULL), interaction(NULL), N(NULL), backend_info(NULL), lists(NULL), box(NULL) {}

	/**
	 * @brief Allows one to initialize everything at once.
	 *
	 * @param p
	 * @param b
	 * @param i
	 * @param Nn number of particles
	 * @param info
	 * @param l pointer to list object
	 */
	ConfigInfo(BaseParticle<number> **p, number *b, IBaseInteraction<number> *i, int *Nn, std::string *info, BaseList<number> *l) : particles(p), box_side(b), interaction(i), N(Nn), backend_info(info), lists(l), box(NULL) {}
	
	ConfigInfo(BaseParticle<number> **p, number *b, IBaseInteraction<number> *i, int *Nn, std::string *info, BaseList<number> *l, BaseBox<number> *abox) : particles(p), box_side(b), interaction(i), N(Nn), backend_info(info), lists(l), box(abox) {}

	virtual ~ConfigInfo() {}

	///Pointer to the array which stores all the particles' information.
	BaseParticle<number> **particles;

	/// Length of the box side.
	number *box_side;

	/// Used to compute all different kinds of interaction energies (total, partial, between two particles, etc.).
	IBaseInteraction<number> *interaction;

	/// Number of particles.
	int *N;

	/// Used by BackendInfo to print backend-related information such as Monte Carlo acceptance ratios.
	std::string *backend_info;
	
	/// Pointer to lists
	BaseList<number> *lists;

	/// Pointer to box object
	BaseBox<number> * box;
};

/**
 * @brief Abstract class defining the basic interface for observables.
 */
template <typename number>
class BaseObservable {
protected:
	/// Stores all the backend's information that may be needed by the observable
	ConfigInfo<number> _config_info;
public:
	BaseObservable() {};
	virtual ~BaseObservable() {};

	/**
	 * @brief Returns the string to be printed in the output stream.
	 *
	 * @param curr_step
	 * @return observable output
	 */
	virtual std::string get_output_string(llint curr_step) = 0;

	/**
	 * @brief Reads the options stored in the input file associated to this observable and in the general input file.
	 *
	 * @param my_inp
	 * @param sim_inp
	 */
	virtual void get_settings (input_file &my_inp, input_file &sim_inp) {}

	/**
	 * @brief Initializes the observable. This basic implementation copies the passed config_info.
	 *
	 * @param config_info
	 */
	virtual void init(ConfigInfo<number> &config_info) { _config_info = config_info; }
};

#endif /* BASEOBSERVABLE_H_ */
