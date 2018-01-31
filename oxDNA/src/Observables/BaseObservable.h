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
#include "../Utilities/ConfigInfo.h"

/**
 * @brief Abstract class defining the basic interface for observables.
 */
template <typename number>
class BaseObservable {
protected:
	/// Stores all the backend's information that may be needed by the observable
	ConfigInfo<number> &_config_info;
public:
	BaseObservable() : _config_info(ConfigInfo<number>::ref_instance()) {};
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
	virtual void init(ConfigInfo<number> &config_info) { }
};

#endif /* BASEOBSERVABLE_H_ */
