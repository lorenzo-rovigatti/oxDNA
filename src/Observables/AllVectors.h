/*
 * all_vectors.h
 *
 * Created on Mar 25, 2019
 *       Author: Poppleton
 */

#ifndef ALL_VECTORS_H_
#define ALL_VECTORS_H_

#include "BaseObservable.h"
#include "../Utilities/OrderParameters.h"

/**
 * @brief Outputs a contact map for the system
 */

class AllVectors: public BaseObservable {
protected:

public:
	AllVectors();
	virtual ~AllVectors();
	virtual void init();
	std::string get_output_string(llint curr_step);
};

#endif /*ALL_VECTORS_H_ */
