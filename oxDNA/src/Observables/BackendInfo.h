/*
 * BackendInfo.h
 *
 *  Created on: Feb 14, 2013
 *      Author: Flavio
 */

#ifndef BACKEND_INFO_H_
#define BACKEND_INFO_H_

#include "BaseObservable.h"

/**
 * @brief Outputs backend-related information (see ConfigInfo).
 */
template<typename number>
class BackendInfo: public BaseObservable<number> {
public:
	BackendInfo();
	virtual ~BackendInfo();

	std::string get_output_string(llint curr_step);
};

#endif /* HBENERGY_H_ */
