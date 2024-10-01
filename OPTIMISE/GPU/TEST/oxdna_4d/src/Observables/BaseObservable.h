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

class BaseObservable {
protected:
	/// Stores all the backend's information that may be needed by the observable
	ConfigInfo *_config_info;

	std::string _id;
	long long int _update_every = 0;
	long long int _times_updated = 0;
public:
	BaseObservable();

	virtual ~BaseObservable();

	bool is_update_every_set();

	virtual bool need_updating(llint curr_step);

	virtual bool require_data_on_CPU() {
		return true;
	}

	virtual void update_data(llint curr_step);

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
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	virtual void serialise() {

	}

	/**
	 * @brief Initializes the observable.
	 */
	virtual void init();

	void set_id(std::string id) {
		_id = id;
	}

	std::string get_id() {
		return _id;
	}
};

using ObservablePtr = std::shared_ptr<BaseObservable>;

#endif /* BASEOBSERVABLE_H_ */
