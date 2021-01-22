/*
 * contact_map.h
 *
 * Created on Mar 25, 2019
 *       Author: Poppleton
 */

#ifndef CONTACT_MAP_H_
#define CONTACT_MAP_H_

#include "BaseObservable.h"
#include "../Utilities/OrderParameters.h"

/**
 * @brief Outputs a contact map for the system
 */

class ContactMap: public BaseObservable {
protected:

public:
	ContactMap();
	virtual ~ContactMap();
	virtual void init();
	std::string get_output_string(llint curr_step);
};

#endif /*CONTACT_MAP_H_ */
